# Fourier maps

Browser tool for Fourier contour maps (Fo-Fc, Fo, Fc) and Patterson maps from a CIF with reflection data.
No server, no build step: unzip, keep the folder structure, open `index.html`.

```
fourier-maps/
  index.html
  static/css/style.css
  static/js/{cif,hkl-parser,scattering-factors,map,app,simulate,orient,profile,patterson}.js
```

Reads `_shelx_hkl_file`, `_refln` and `_diffrn_refln` data. Fc is computed from the model in the CIF
(IT92 form factors, anisotropic ADPs, f', f'', special positions); weights follow `_refine_ls_weighting_details`.

Live view: after loading a CIF the Fo and Fc densities are computed once as a 3D grid (FFT, a few tenths of a second).
Sections (plane atoms, offset, width, grid, map type, contour step, colours) are then cut out of that grid by tricubic
interpolation in a few milliseconds, so the map follows every control without a calculate button. Only a change of the
resolution limit recomputes the 3D grid.

Mouse: drag the map to pan, wheel to zoom around the pointer, double-click to reset the view.

Performance: while dragging or moving a slider a cheap preview is drawn (coarser grid, trilinear); the full-quality map follows 150 ms after the last change. Atoms are generated for the visible region only, so panning never runs out of atoms.

The map type *2Fo − Fc* (relative) is the Fo map plus the Fo − Fc map: atoms of the model keep their height, missing atoms come out about twice as strong as in the Fo map. With a partial model it shows the framework and the gaps in one picture.

## Turning the section

The plane atoms of the list are only a starting point. **Cmd / Ctrl / Alt + drag** turns the section around the midpoint of its three atoms (dragging over the whole width = 180°); the map is cut from the 3D grid, so any plane works. While the key is held (or with *Show molecule*) the bonded molecule around the centre is drawn as a pale skeleton with atom names, also the parts outside the atom slab, so you can see which atoms you are looking at. **Cmd / Ctrl / Alt + click** picks atoms (the clicked symmetry images); the third pick defines a new plane and brings it into the picture plane. *Reset tilt* turns back; choosing atoms from the lists goes back to a plane through those atoms. The log shows the angle against the atom plane. Plain click / Shift+click still add / remove atoms of the model in the simulated refinement.

## Profile

**Shift + drag** on the map from atom to atom draws a line; the ends snap to atoms (also to atoms of the molecule beside the section, which is shown while Shift is held). A strip below the map shows Fo, Fc and Fo-Fc along the line (e/Å³ over Å, with the atoms near the line marked). The profile is cut from the 3D grid, so a bond that leaves the plane works too. Extension continues the line beyond both atoms. Hovering the strip moves a marker on the line in the map and vice versa; the × at the top right (or *Clear line*) removes it. The log lists Fo, Fc and Fo-Fc at the atom centres. The curves follow the map type (absolute types include F000) and the partial model of the simulated refinement.

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
- The Model card moves with the tab. In the Patterson tab the Fc² and Fo²−Fc² maps are computed from the partial model, so Fo²−Fc² shows the vectors that are still missing.

## Layout

Like the other tools: the controls are in a sidebar on the left (full window height, scrolls on its own), the map uses the rest of the window. The cards in the sidebar can be folded (click the title; the state is remembered in the browser). The Model card moves with the tab and is the last card of the sidebar. Drop a CIF anywhere in the window (or use *Open CIF*).

## Patterson tab

Second tab next to the Fourier map. P(u,v,w) = 1/V · Σ |F|² cos 2π(hu+kv+lw) (F000 left out), computed once as a 3D grid like the Fourier density.
- Coefficients: |F|² (origin peak = 100 %) or sharpened E²−1 (Wilson-normalised, units of σ); maps from Fo², Fc² or Fo²−Fc².
- Sections u, v or w = const at any position (0 … 1); the cursor shows u v w and the map value; the log lists the highest peaks of the section.
- Harker elements are derived from the symmetry operators of the CIF (rotation/screw axes → Harker planes, mirror/glide planes → Harker lines)
  and drawn into the section (orange dashed = plane, green = line / point; thick orange frame = the section itself is a Harker section).
  The element buttons jump to a matching section (e.g. P2₁/c: v = ½ plane and [0, v, ½] line).
- Interatomic vectors: for every pair of atoms (all symmetry images) the vector is drawn as a ring when it lies within 0.4 Å of the section; ring size ~ Z·Z. Purple = atom with its own symmetry image (these sit on the Harker elements), blue-green = other vectors, dashed = involves an atom that is not in the model. *Min. weight* hides weak vectors; the atom pairs can be limited to the model or to vectors with a missing atom. Hovering a ring shows the atom pair.
