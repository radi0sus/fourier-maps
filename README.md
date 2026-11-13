# Fourier maps

Browser tool for Fourier contour maps (Fo-Fc, Fo, Fc) and Patterson maps from a CIF with reflection data.
No server, no build step: unzip, keep the folder structure, open `index.html`.

```
fourier-maps/
  index.html
  static/css/style.css
  static/js/{cif,hkl-parser,scattering-factors,map,app,simulate,patterson}.js
```

Reads `_shelx_hkl_file`, `_refln` and `_diffrn_refln` data. Fc is computed from the model in the CIF
(IT92 form factors, anisotropic ADPs, f', f'', special positions); weights follow `_refine_ls_weighting_details`.

Live view: after loading a CIF the Fo and Fc densities are computed once as a 3D grid (FFT, a few tenths of a second).
Sections (plane atoms, offset, width, grid, map type, contour step, colours) are then cut out of that grid by tricubic
interpolation in a few milliseconds, so the map follows every control without a calculate button. Only a change of the
resolution limit recomputes the 3D grid.

Mouse: drag the map to pan, wheel to zoom around the pointer, double-click to reset the view.

Performance: while dragging or moving a slider a cheap preview is drawn (coarser grid, trilinear); the full-quality map follows 150 ms after the last change. Atoms are generated for the visible region only, so panning never runs out of atoms.

## Simulated refinement (Model card, Fourier tab)

Teaching mode, no real refinement. Switch on *Simulate refinement*: only the atoms ticked in the model list contribute to Fc
(start: the heaviest element). The other atoms are still in Fo, so they show up as peaks in Fo-Fc. Add atoms step by step
(click a peak in the map, click an atom in the list, or *+ strongest peak*; Shift+click removes an atom) and watch the peaks
disappear and R1 / wR2 fall. Positions and ADPs always come from the CIF.
- Fo keeps the scale of the full model, so the difference map is easy to read. *Refit scale* fits the scale factor again after every step.
- The phases of the Fo map come from the partial model (model bias, as in a real structure solution). The H atoms stay as weak peaks after *Non-H*.
- *Auto contour step* (on by default) sets the contour step after every change of the model, so that the strongest missing peak gets about 12 lines; moving the contour slider by hand switches it off.
- Missing atoms can be drawn as dashed rings (*Show missing atoms*) or hidden. The log lists the strongest Fo-Fc densities at the missing
  atoms and the R1 of the last steps. *Undo* goes back one step.
- The Patterson tab is not affected by the partial model.

## Patterson tab

Second tab next to the Fourier map. P(u,v,w) = 1/V · Σ |F|² cos 2π(hu+kv+lw) (F000 left out), computed once as a 3D grid like the Fourier density.
- Coefficients: |F|² (origin peak = 100 %) or sharpened E²−1 (Wilson-normalised, units of σ); maps from Fo², Fc² or Fo²−Fc².
- Sections u, v or w = const at any position (0 … 1); the cursor shows u v w and the map value; the log lists the highest peaks of the section.
- Harker elements are derived from the symmetry operators of the CIF (rotation/screw axes → Harker planes, mirror/glide planes → Harker lines)
  and drawn into the section (orange dashed = plane, green = line / point; thick orange frame = the section itself is a Harker section).
  The element buttons jump to a matching section (e.g. P2₁/c: v = ½ plane and [0, v, ½] line).
