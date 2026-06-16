# Interactive Spectroscopy Analysis

An interactive Jupyter notebook for fitting peaks in Raman, fluorescence, and
other spectroscopic data. Built for Diamond Light Source beamline I18 and
shared with students and collaborators. It runs in JupyterLab with `ipywidgets`
and provides a point-and-click workflow: load one or many spectra, subtract a
baseline, detect and fit peaks (Gaussian / Lorentzian / Pseudo-Voigt), and
export the results.

## Features

- **Multi-file fitting** — load many spectra at once; results are tabulated
  with a common set of peak columns so they line up across files, plus a
  centre-vs-FWHM summary scatter.
- **Background subtraction** — arPLS (default) or iterative polynomial, fitted
  per spectrum after edge-trimming.
- **Peak fitting** — Gaussian, Lorentzian, Pseudo-Voigt, or shared-η
  Pseudo-Voigt, with analytical Jacobians for speed. Add shoulder peaks by hand
  or let the residual auto-suggest them.
- **Wavenumber calibration** (optional) — for spectra recorded with Hg(Ar) / Kr
  calibration lamps switched on, correct each spectrum's x-axis on the fly with
  a polynomial fit, with a diagnostic plot and model-selection statistics to
  guard against over-fitting.
- **Spectral decomposition** (optional) — Non-negative Matrix Factorisation to
  pull end-member spectra and abundances out of a mixed dataset.
- **Export** — save the fit figure and a self-documenting results table.

## Getting started (no Python experience needed)

1. **Install Anaconda** (one time) from <https://www.anaconda.com/download>.
   Run the installer with the default options. This provides Python plus most
   of what the notebook needs.

2. **Install the extra packages.** Open *Anaconda Prompt* (Windows) or
   *Terminal* (macOS) and run:
   ```
   pip install itables
   ```
   (`itables` powers the interactive results tables; everything else ships with
   Anaconda. To install the full list explicitly, use
   `pip install -r requirements.txt`.)

3. **Launch JupyterLab** from Anaconda Navigator (click *Launch* under the
   JupyterLab tile), or run `jupyter lab` in the terminal.

4. **Open `spectroscopy.ipynb`** in the JupyterLab file browser.

5. **Run the cells** — click the first cell and press **Shift+Enter** to run
   each in turn (four cells). The last cell displays the interactive interface.

There is a terse *Quick reference* section just above the app explaining every
control. To get going: click *Upload*, choose one or more two-column
(x, y) spectrum files (`.txt` / `.csv` / `.dat`; header rows are fine), and
tweak *Threshold %* until the right peaks light up.

## Updating

This repository is the canonical, newest version. To grab the latest:

```
git pull
```

or download it fresh from the green **Code → Download ZIP** button on GitHub.

## Troubleshooting

- A red error after running a cell is usually a data-format issue — check the
  file is two numeric columns with any header rows above.
- If the sliders or plots don't appear, restart the kernel
  (*Kernel → Restart Kernel and Run All Cells*).

## Dependencies

numpy, pandas, matplotlib, scipy, ipywidgets, itables, and (only for the
optional spectral-decomposition section) scikit-learn. See `requirements.txt`.
