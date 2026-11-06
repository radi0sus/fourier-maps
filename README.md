# Fourier maps

Browser tool for Fourier contour maps (Fo-Fc, Fo, Fc) from a CIF with reflection data.
No server, no build step: unzip, keep the folder structure, open `index.html`.

```
fourier-maps/
  index.html
  static/css/style.css
  static/js/{cif,hkl-parser,scattering-factors,map,app}.js
```

Reads `_shelx_hkl_file`, `_refln` and `_diffrn_refln` data. Fc is computed from the model in the CIF
(IT92 form factors, anisotropic ADPs, f', f'', special positions); weights follow `_refine_ls_weighting_details`.

Live view: after loading a CIF the Fo and Fc densities are computed once as a 3D grid (FFT, a few tenths of a second).
Sections (plane atoms, offset, width, grid, map type, contour step, colours) are then cut out of that grid by tricubic
interpolation in a few milliseconds, so the map follows every control without a calculate button. Only a change of the
resolution limit recomputes the 3D grid.

Mouse: drag the map to pan, wheel to zoom around the pointer, double-click to reset the view.

Performance: while dragging or moving a slider a cheap preview is drawn (coarser grid, trilinear); the full-quality map follows 150 ms after the last change. Atoms are generated for the visible region only, so panning never runs out of atoms.
