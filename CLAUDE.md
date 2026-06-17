# Spectroscopy Notebook — Project Context

## Overview

`spectroscopy.ipynb` is an interactive Jupyter notebook for fitting peaks in Raman, fluorescence, and other spectroscopic data. It is designed for use at Diamond Light Source beamline I18 and by students/collaborators. The notebook runs in JupyterLab with `ipywidgets` and `%matplotlib inline`.

## Notebook structure

The notebook has 4 cells (plus a trailing empty cell):

1. **Cell 0 — Markdown**: Title and high-level workflow summary.
2. **Cell 1 — Pure functions**: All algorithms and I/O — no widgets, no side effects, no `display()`. Baseline estimation, peak profiles, model evaluation, peak detection / fitting helpers, pandas DataFrame builders for the parameter tables, file loading, and results export.
3. **Cell 2 — Markdown**: Quick reference for the GUI controls (kept terse so it doesn't crowd the app).
4. **Cell 3 — Interactive app**: All `ipywidgets` construction, callbacks, layout, and `display(app)`. Depends on functions defined in cell 1.

This separation is deliberate: keep logic testable in cell 1, keep UI wiring in cell 3.

## Key architecture patterns

### State dictionary
All shared mutable state lives in a single dict called `S`. The notebook supports loading multiple spectra at once, so per-file data lives inside `S['files']`:
```python
S = dict(files=[], n_peaks=0, fig=None, fig2=None, _freeze=False)
# each file dict carries its own arrays + fit results:
#   name, x, y, x_trim, y_raw_trim, bl_trim, y_sub,
#   fit_curve, residual, popt, perr, peak_kept
```
`_freeze` is a transient flag used by "Tick all"/"Untick all" to suppress per-checkbox redraws during a bulk toggle. `fig2` holds the secondary "centre vs FWHM" summary plot. Callbacks read widget values directly but write results into `S`. Plotting and parameter display read from `S`.

### Peak detection cap
`detect_peaks()` honours a module-level `MAX_PEAKS` constant (default 30). When the threshold slider drops very low and `find_peaks` returns hundreds of candidates, the cap keeps `curve_fit` responsive — the tallest `MAX_PEAKS` are kept and the rest are discarded.

### Section headers
GUI sections use a lambda for consistent styling:
```python
hdr = lambda t: widgets.HTML(
    f'<h4 style="margin:8px 0 4px; border-bottom:1px solid #ccc">{t}</h4>')
```

### Callback wiring
All parameter widgets observe `update_plot()` via a loop:
```python
for _w in [lam_slider, poly_deg, trim_slider, profile_rb, ht_slider,
           range_slider, ref_dropdown]:
    _w.observe(lambda c: update_plot(), names='value')
```
Button clicks use `.on_click()`. The `on_load` callback handles both `FileUpload` and text path input.

### Range slider safety
`_update_range_slider()` batches the bound + value updates inside `range_slider.hold_sync()` so JupyterLab cannot display an intermediate state where the value has been clamped by stale bounds. Sequence: widen bounds → set value → narrow bounds to the exact target → step. Without `hold_sync`, the frontend would occasionally render `(xmin, xmin)` because the value-update sync arrived before the max-update sync.

The width gate in `peak_kept` also has a defensive fallback: if `range_slider.value` happens to be degenerate (`x_hi - x_lo < 1e-9`), it falls back to the per-spectrum `x_trim` extent so the gate doesn't silently disable itself.

## Algorithms

### Baseline estimation
- **arPLS** (`arpls()`): Asymmetrically Reweighted Penalized Least Squares. Single parameter λ controls smoothness. Default λ=1e4. Iterates until weight convergence. Exponent is clipped to [-500, 500] for numerical stability.
- **Polynomial** (`poly_baseline()`): Iterative polynomial fit. Each iteration removes points above the current fit, then refits. Default degree=5, up to 50 iterations.

The baseline is computed on the **edge-trimmed** spectrum (the *Edge trim* slice is applied to the raw data *before* `arpls`/`poly_baseline`), so instrumental edge ramps in extended-range scans can't drag the baseline below the real signal.

### Peak profiles
All profiles parameterised as (amplitude, centre, FWHM[, eta]):
- `gaussian(x, amp, cen, fwhm)` — σ derived from FWHM
- `lorentzian(x, amp, cen, fwhm)` — γ = FWHM/2
- `pseudo_voigt(x, amp, cen, fwhm, eta)` — linear combination: η·Lorentzian + (1−η)·Gaussian

### Multi-peak models
- `multi_gauss`, `multi_lor`: stride 3 (amp, cen, fwhm per peak)
- `multi_pv`: stride 4 (amp, cen, fwhm, eta per peak)
- `multi_pv_shared`: stride 3 per peak + 1 shared eta at the end of the parameter vector

All four models are vectorised across peaks via numpy broadcasting (single `(N_x, n_peaks)` grid, no Python loop). Each model has an analytical Jacobian counterpart (`jac_gauss`, `jac_lor`, `jac_pv`, `jac_pv_shared`) passed to `curve_fit` via `jac=...`. The model + analytical-Jacobian combination delivers a ~50× speedup over the original Python-loop + numerical-Jacobian implementation on 30-peak fits.

### Peak detection
Uses `scipy.signal.find_peaks` with height and prominence thresholds derived from the `Threshold %` slider (percentage of the post-baseline data range, computed by `masked_data_range` which excludes ±8 cm⁻¹ around any `exclude_centres` so a bright calibration lamp line can't inflate the range and push the threshold above genuine Raman peaks), with a hard `NOISE_K * sigma_noise` floor (default 5σ) so peakless / dark spectra return zero peaks. `estimate_noise_sigma` uses MAD of the first difference (robust to real peaks and slow baseline residuals). `detect_peaks` returns `(idx, info)` where `info` is a diagnostic dict consumed by the UI's live peak-count readout. Optional `exclude_centres` / `exclude_window` arguments drop any candidate within ±window of any listed centre — used to remove calibration-lamp emission lines from the Raman peak set after the wavenumber correction step. The detection exclusion uses the module constant `LAMP_MASK_HW` (8 cm⁻¹) — **the same half-width the fit later masks** (`fit_spectrum(fit_mask=...)`). The two must match: if detection excluded a *narrower* band than the fit mask, a feature picked up in the gap (outside detection, inside the fit mask) would enter the peak union and then be fitted against points that have all been masked away, producing a spurious "ghost" peak with a wildly batch-/window-dependent amplitude (this was a real bug — a detection at ~475 cm⁻¹, 3.5 cm⁻¹ from the Hg 479 lamp line, slipped past the old ±3 detection window but sat inside the ±8 fit mask).

### Wavenumber calibration (optional)
For the workflow where Hg(Ar) / Kr calibration lamps are on during acquisition so their sharp emission lines appear in every spectrum. Three helpers in cell 1:
- `parse_cal_lines(text)` — parses `label  wavenumber` lines (one per row, `#` for comments) from the textarea.
- `fit_calibration_lines(x, y, expected, window)` — for each expected (label, known) locates the local maximum within ±window of `known`, fits a single Gaussian for sub-sample precision, returns per-line dicts (`label, known, measured, fwhm, success, error`).
- `parse_cal_lines(text)` returns `(label, wavenumber, wavelength)` tuples — each input row is `label  wavenumber  [wavelength_nm]`; the optional wavelength is informational (shown in the diagnostic table). Leading non-numeric tokens are the label, first numeric is the operative wavenumber.
- `compute_wavenumber_correction(found, allow_offset=False, order='auto')` returns a **dict** (or None): `coeffs, order, rms, residuals, n_used`, the covariance machinery `x0, s, invVtV, sigma` (for the position-uncertainty band), per-order model-selection `stats`, and `knowns/measured/cen_err`. `order` may be `'auto'` (rule: 1→offset, 2–3→linear, 4+→quadratic), `'aicc'` (minimise AICc — the right call for few noisy lines; AICc is `inf` once params reach n−1, so it refuses unsupported orders), or an explicit int clamped to `min(order, 3, n−1)`. `stats[o]` holds `rms, chi2_red, aicc, loocv` — **don't pick order by RMS** (monotone); use χ²_red≈1, AICc/LOOCV minima.
- `cal_position_sigma(cal, x)` — 1σ calibration uncertainty at wavenumber(s) x, `sigma·sqrt(v(x)·invVtV·v(x))` in scaled coords. This systematic is propagated into the parameter tables (separate `Centre err (cal)` column alongside `Centre err (fit)`) and folded into the summary plot's centre-axis error bars in quadrature.
- `apply_wavenumber_correction(x_raw, coeffs)` — `x_raw - np.polyval(coeffs, x_raw)`.

The default lamp-line list spans the extended range: Hg(Ar) 479, Kr 813.2/839.3, Hg 1459.4/1522.4 cm⁻¹ (the Hg doublet from 576.96/579.07 nm, computed for the fitted 532.15 nm laser). The `cal_order_dd` dropdown (Auto/Offset/Linear/Quadratic/Cubic) drives the `order` argument.

`fit_calibration_lines` returns a per-line `cen_err` (Gaussian-fit centre uncertainty) used by the diagnostic plot. `_draw_cal_diagnostic()` renders, for the explored file, measured-vs-known (left) + residuals (right) with a **2σ confidence band propagated from `cen_err` through the fit** (not from residuals — so it survives exact fits and balloons under overfitting as the order rises). It's an `Output` (`cal_diag_out`) sitting to the right of the controls in `cal_body`, refreshed by the cal step and by the Explore dropdown.

In cell 3, the calibration step runs at the top of `update_plot` whenever `cal_use_cb.value` is True. Each file keeps an immutable `fd['x_raw']` set at load time; the working `fd['x']` is always re-derived from `x_raw` (either a copy when cal is off, or corrected when on). Status is reported per-file in `cal_status_html`. `fd['cal_info']` includes `coeffs`, `order`, `mode` (`offset`/`linear`/`quadratic`/`cubic`), `rms`, and `max_resid`. The status is green only when all expected lines are used AND `rms < 0.3` (the lamp lines lie cleanly on the fitted curve); an RMS that stays large as the order rises flags a genuine stitch step rather than smooth distortion. The status shows `N/M lines (<mode>), shift ±X at midpoint, RMS Y (max Z)`, except an exact fit (`n == order+1`, e.g. 3 lines + quadratic) hides the trivially-zero RMS. When **fewer lines than expected are used**, the status appends `— unused: <label> <known> [<reason>]` per dropped line, and the per-line diagnostic table (`cal_table_html`) carries a `status` column (`SNR <n>` for a used line; the rejection reason — `SNR <n> < <min>`, `no data in ±window`, or `fit failed` — in bold red for a dropped one). The reason text comes from `cal_line_note(d)` (cell 1). This makes a present-but-rejected lamp line (e.g. Hg 479 falling below the SNR gate in a particular spectrum) visible instead of silently absent. The known centres are also passed to `detect_peaks` as `exclude_centres` so they don't pollute the Raman peak union, and drawn as blue dashed vertical guides on the explorer / main / residual panels (separate style `_CAL_LINE` from the grey-dotted `_GUIDE_LINE` used for Raman peak centres).

### Multi-file fitting strategy (two-pass)
Designed for datasets where the spectra may be very different (e.g. a transect across a layered sample), so no single file can be treated as canonical.

1. **Pass 1 — per-file detection.** `detect_peaks` runs independently on every loaded file. `union_peak_centres` then merges the resulting centres across files via single-linkage clustering with `merge_tol = 5 * median_dx`; the cluster median is used as each union centre. If the union exceeds `MAX_PEAKS`, clusters are ranked first by the number of distinct files contributing (peaks that appear in many files survive), then by mean detected height.
2. **Pass 2 — per-file fit.** Each spectrum is fitted to the *union* peak set, with `snap_to_local_max` providing per-file initial centres/amplitudes from each spectrum's own data. `build_initial_guess` is rebuilt per file so amplitude bounds reflect what's actually in that spectrum. Peaks absent in a given file end up with tiny fitted amplitudes and are blanked by the per-spectrum `peak_kept` gates.

**Known consequence — results are mildly batch- and window-dependent by design.** Because the union peak *set* is shared across whichever files are loaded (and whichever x-range is active), a given spectrum's fit is not guaranteed bit-identical when the file set or x-range changes:
- *Isolated peaks are stable.* With plain (independent-η) Pseudo-Voigt, an isolated peak's centre and FWHM are essentially invariant to the batch and to the window (≈0.00 cm⁻¹ in tests on the WASH-WW-MYT transect).
- *Overlapping peaks shift when a neighbour enters/leaves the union.* If a peak is detected in only a few files (e.g. a shoulder present at one end of a transect), it still joins the union and is fitted into *every* file. Where it sits next to a real band, fitting it (even though it's later blanked) steals wing intensity and shifts the real band's centre/FWHM. Observed: a 1101 cm⁻¹ shoulder present in 2/15 files pulled the 1087 band by ~1.4 cm⁻¹ in the 13 files lacking it.
- *Shared-η couples everything.* `Pseudo-Voigt (shared η)` fits one global η for all peaks, so changing the peak set (more files, or a cropped window) changes η and hence *every* peak's FWHM — even isolated ones. Plain Pseudo-Voigt avoids this; prefer it when per-peak FWHM/position stability matters.

#### Per-file peak gate (`support_gate_cb`, opt-in, off by default)
The neighbour-stealing above is addressed by an opt-in checkbox in Section 4. When on, each spectrum is fitted in Pass 2 only to the union peaks that have **local support** in *that* spectrum — defined as a peak detected at the bare 5σ noise floor (`detect_peaks(..., height_pct=0)`) within `merge_tol` of the union centre (manual peaks are always kept). Unsupported peaks are dropped from that file's parameter vector *before* the fit, then `_scatter_popt` expands the subset result back to the standard full-length `popt`/`perr` layout with the absent peaks at zero amplitude — so every downstream gate / plot / table / export indexes it exactly as an ungated fit, and a real peak wrongly dropped still surfaces in the residual. `_fit_supported` does the support test + subset fit + scatter; `fd['support_n']` records how many peaks were actually fitted, and the live count readout appends "· per-file gate on".

Validated on the WASH-WW-MYT line_1 transect (15 files, calibration on): with the gate **on**, a given spectrum's affected bands fit *identically* to fitting that file alone (e.g. the 1087 band: 1088.4 cm⁻¹ ungated → 1086.92 cm⁻¹ gated = the single-file value, because the 1101 shoulder present in only 2/15 files is no longer force-fitted into the other 13). The support matrix on that data is sharply bimodal (files 01–04 carry 280/713/1101/1486; files 05–15 carry 205/703/1087), i.e. a genuine two-layer transect, so the gate's per-file peak switching reflects real structure. **Trade-off:** for a peak that fades in *gradually* along a transect, the gate switches it on at the file where it crosses 5σ, which can make a smooth trend render as a step — hence opt-in, not default.

### Manual peaks
`parse_manual_peaks` reads extra centres from the *Manual peaks* textarea; they're merged into `union_xs` (skipping any within `merge_tol` of an auto centre) for shoulders / blended peaks that aren't local maxima and so can't be found by `find_peaks`. `S['manual_mask']` flags them so pass-2 uses their exact centre as the initial guess instead of `snap_to_local_max` (which would collide a shoulder with its neighbouring main peak). They still pass through the `peak_kept` gates. An optional *Auto-add peaks from residual* mode (`residual_auto_cb` → `_augment_union_from_residuals`, default `max_iter=1`) adds peaks found in the fit residual that clear `k_sigma`·noise (default 5) and lie >`merge_tol` from every existing peak (the distance test rejects derivative-shaped residuals of mis-modelled peaks). A single residual pass is the default deliberately: the first residual exposes all genuinely-separated missing peaks at once, whereas iterating over-decomposes an asymmetric band into spurious extra components chasing its non-PV tail. Off by default.

### `peak_kept` gates (two per peak, per spectrum)
- **Amplitude:** `popt[amp] >= Threshold % × data range` of that spectrum (range via `masked_data_range`, lamp regions excluded).
- **Width range:** `width_slider.value` is a `(min%, max%)` tuple; the gate is `(width_min × x_range) <= popt[fwhm] <= (width_max × x_range)`. The **max** end catches `curve_fit` runaway-broadening fits where a pseudo-Voigt swallowed noise. The **min** end catches the symmetric collapse-onto-noise-spike failure where `curve_fit` drove the FWHM down to its lower bound to fit a single noisy sample, producing a tall narrow false peak that still passes the amplitude gate.

### Single-spectrum table
`build_single_file_df` honours `peak_kept` — it lists only peaks that passed the gates for that file. Without this, hiding all but one file would show the full union list (confusing because most rows correspond to peaks not present in that spectrum). `_draw_params` prints "No peaks survived the threshold / width filters for this spectrum." when the filtered DataFrame is empty.

### Blanked peaks in plots / residuals
Each spectrum carries both `popt` (the raw `curve_fit` output, used by the parameter table) and `popt_plot` (a copy with the amplitudes of `peak_kept`-False peaks zeroed). All plotting routines (`_draw_explorer`, `_draw_main_plot` single-spec and waterfall) use `popt_plot` so blanked peaks don't draw on screen, and `fd['fit_curve']` / `fd['residual']` are recomputed from `popt_plot` so the bogus signal surfaces in the residual strip as unfit instead of being silently captured by a rejected peak. The saved-data output uses the same `fit_curve`, so blanked peaks are excluded from exports for consistency.

The N peak columns of the results table line up across files so corresponding peaks can be compared directly. When calibration is on, the explorer / waterfall / single-spec panels and their residual strips set their y-limits from `_masked_ylim` (data with ±8 cm⁻¹ around each lamp centre removed) and the waterfall offset step uses `masked_data_range`, so tall lamp lines run off the top instead of squashing the sample peaks. The residual strips additionally NaN-blank the residual within ±6 cm⁻¹ of each lamp centre (`_blank_lamps`) so the huge unfit-lamp-line residual spikes aren't drawn as clipped vertical streaks. The `Explore` dropdown is purely cosmetic — it picks which spectrum's raw + baseline appears in the explorer panel; the fit is independent of that choice.

### Reload safety
`on_load` wraps the widget-mutation block in `S['_freeze'] = True` so that observer cascades fired by `_update_range_slider`, `trim_slider.max = ...`, and `ref_dropdown.options = ...` cannot run `update_plot` against transiently-mismatched state (`S['files']` vs stale `files_box.children`). `update_plot` and `update_summary_plot` additionally clamp `visible_idxs` to `< len(files)` as defensive belt-and-braces.

### Fitting
`scipy.optimize.curve_fit` with bounds, `jac=<analytical>`, and `x_scale='jac'`. When calibration is on, `fit_spectrum(..., fit_mask=...)` drops the lamp-line regions (±`LAMP_MASK_HW` = 8 cm⁻¹) from the fitted points so their huge unmodelled residuals can't dominate the least-squares cost or pull broad spurious peaks into the gaps between bright lamp lines (the model is still evaluated on the full axis for plotting/residuals). The same `LAMP_MASK_HW` is used by peak detection (above) and the residual-augment helper, so a region masked from the fit is never detected as a Raman peak. Initial guesses: amplitude from local peak height, centre from local peak position, FWHM from 5× median point spacing, eta=0.5. Bounds allow centres to shift ±15×dx and FWHM up to 150×dx. maxfev=30000.

## File I/O

### Loading (`auto_load_spectrum`)
Auto-detects header lines by scanning for the first row with two parseable floats. Tries tab, comma, and whitespace delimiters. Returns x, y arrays plus skip/nrows metadata.

### Saving (`save_fit_results`)
Writes a self-documenting text file with a comment header (source file, date, baseline method, profile, peak table) followed by tab-separated columns: x, raw, baseline, subtracted, fit, residual.

## GUI layout (cell 2)

The app is a `VBox` with seven numbered sections:

1. **Load spectra**: `FileUpload(multiple=True)` + multi-line `Textarea` for paths/globs + "Load paths" button; auto-detect checkbox; skip/max rows. Below: scrollable `VBox` of per-file visibility checkboxes, "Tick all"/"Untick all" buttons, and the *Explore* dropdown (formerly *Reference*) that picks which spectrum is shown in the **Raw spectrum explorer** panel. The dropdown has **no effect on the fit** — it only triggers `_draw_explorer()`.
2. **Wavenumber calibration** (optional, off by default): `cal_use_cb` toggle reveals/hides `cal_body` (a `Textarea` of `label  wavenumber` rows + a search-window slider + an "Allow single-peak offset" checkbox + a status `HTML`). When enabled, every spectrum's x-axis is corrected on the fly from its `x_raw` snapshot; the known centres are excluded from Raman peak detection and drawn as blue dashed guides on the plots.
3. **Background subtraction**: `RadioButtons` (arPLS/Polynomial) beside a `VBox` of *Smoothness λ* slider, *Poly degree* slider (hidden when arPLS), *Edge trim* slider. Each file gets its own baseline; the explorer panel renders the selected file's raw + baseline + model on top of one another.
4. **Peak fitting**: `RadioButtons` (4 profiles) beside *Threshold %* slider (with the live peak-count readout `peak_count_html` to its right), `Width %` range slider, x-range slider with reset button, and *Waterfall offset* slider (cosmetic y-spacing).
5. **Summary plot**: `summary_plot_out` (centre-vs-FWHM scatter, sized identically to the explorer above so peak x-positions line up vertically) followed by a `peak_select_box` `GridBox` of `Checkbox` toggles (one per useful union peak; 5-column grid wrapping to multiple rows above 5 peaks).
6. **Spectral decomposition** (optional, off by default): `decomp_use_cb` toggle reveals/hides `decomp_body`. NMF (via scikit-learn) is run on the baseline-subtracted spectra of all currently-visible files, on a common x-grid (intersection of x-ranges, interpolated). End-member spectra are rendered in `decomp_spectra_out` (stacked, with union Raman peak centre guides); per-spectrum weights are a stacked bar in `decomp_weights_out`. Re-runs automatically whenever upstream parameters change (via call from `update_plot`) and when the components slider / normalise checkbox change.
7. **Export**: Save plot (PNG) and save data buttons. *Save data* writes a single-spectrum text file when only one is visible, otherwise a multi-spectrum results-table TSV.

Output areas, ordered top-to-bottom on the page with shared width (~1100 px) so peak x-positions can be traced vertically:
- `explorer_out` — single-spectrum **Raw spectrum explorer**: raw + baseline + (model + baseline) + individual peaks (drawn on top of the baseline so peaks sit in their raw-data positions). Figure size `(FIG_W=11, PANEL_H=5.14)` matches the single-spec main plot.
- `plot_out` — main plot: waterfall (or single-spec data + fit) + residual strip. Two panels with `height_ratios=[main_h, RESIDUAL_H]` where `main_h = PANEL_H` for one spectrum or scales with N for the waterfall.
- `params_out` — parameter table via `itables.show`, searchable/sortable, `max_width='1080px'`.
- `summary_plot_out` — same `(FIG_W, PANEL_H)` as the explorer; `xlim` tied to `range_slider.value` so its x-axis matches the explorer's.

All sliders share `_SLIDER_LAYOUT = Layout(width='540px')` and `_SLIDER_STYLE = {'description_width': '110px'}` (range slider is slightly narrower so it can sit beside the Reset button) so labels and tracks line up.

Peak toggles: `peak_select_box` is a `GridBox` rebuilt by `_rebuild_peak_checkboxes(n_pk)` whenever the union peak count changes. The rebuild preserves ticked/unticked state for surviving peak indices and defaults new ones to ticked. `_selected_peaks()` returns the currently-ticked indices (replaces the old `peak_select.value`).

### Spectral decomposition
`build_decomp_matrix(spectra)` stacks visible files onto a common x-grid (intersection of x-ranges, interpolated linearly, finest sampling) and clips negatives so NMF accepts the input. `decompose_nmf(Y, n_components, normalize=True)` runs `sklearn.decomposition.NMF` with `init='nndsvda'`, then rescales each end-member row of H to unit area and absorbs the scale into W so each W[i,k] is the share of spectrum i accounted for by component k. Returns `(H, W, r2)` where r2 is the reconstruction quality. sklearn is imported lazily inside `decompose_nmf` so the rest of the notebook stays usable without it.

## Dependencies

- numpy, pandas, matplotlib, scipy (sparse, signal, optimize), ipywidgets, IPython.display, itables
- scikit-learn — only imported when spectral decomposition (Section 6) is enabled
- No external packages beyond standard scientific Python plus `itables` for the interactive parameter tables

## Conventions

- British English spelling in user-facing text (e.g. "analysed", "colour")
- Widget `style={'description_width': '...px'}` for alignment
- `widgets.Layout(width='...px')` for sizing — no percentage widths
- `clear_output(wait=True)` before redrawing to avoid flicker
- Tick marks: `ax.tick_params(which='both', top=True, right=True, direction='in')` on all axes
- No banner separators (no `=====` lines) in output text

## Planned features

The owner (Ollie) has mentioned wanting to add:
- A file browser GUI widget for searching and loading data files from a directory (the current Textarea + globs is a placeholder for this)

Already implemented (kept here for reference, since CLAUDE.md previously listed them as planned):
- Batch / multi-file fitting with summarised results table and centre-vs-FWHM scatter — done.

## Related notebook

`diffraction.ipynb` in the same folder is a sibling project — an interactive XRD analysis tool with similar widget patterns (same `hdr()` convention, same `S` dict pattern, same plot styling). It supports Cu and NaCl B1 pressure calibrants with BurnMan EoS integration, and includes a gas cylinder depletion calculator. Code structure is analogous but assembled from separate source files via a build script.
