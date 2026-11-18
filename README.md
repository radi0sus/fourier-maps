# Fourier maps

Browser tool for Fourier contour maps (Fo-Fc, Fo, Fc) and Patterson maps from a CIF with reflection data.
No server, no build step: unzip, keep the folder structure, open `index.html`.

```
fourier-maps/
  index.html
  static/css/style.css
  static/js/{cif,hkl-parser,scattering-factors,map,app,simulate,orient,profile,readout,info,patterson,touch,fullscreen}.js
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

## Structure information and read-out

Under the file field: a short line (formula, space group, Z) and the **Structure** card with cell, space group (with number, Hall symbol, crystal system, number of operations, centrosymmetric or not), composition (formula, formula weight, Z, F(000), calculated density, content of the asymmetric unit), the measurement (temperature, radiation, reflections, resolution, Rint) and the refinement (R1, wR2, GooF, parameters; values of this tool next to those of the CIF). *Copy* puts all lines on the clipboard. Missing items are left out.

Under the map: a live read-out with the fractional coordinates at the pointer and the density of the shown map (Fourier tab), or the Patterson value (Patterson tab).

## Turning the section

The plane atoms of the list are only a starting point. **Cmd / Ctrl / Alt + drag** turns the section around the midpoint of its three atoms (dragging over the whole width = 180°); the map is cut from the 3D grid, so any plane works. While the key is held (or with *Show molecule*) the bonded molecule around the centre is drawn as a pale skeleton with atom names, also the parts outside the atom slab, so you can see which atoms you are looking at. **Cmd / Ctrl / Alt + click** picks atoms (the clicked symmetry images); the third pick defines a new plane and brings it into the picture plane. *Reset tilt* turns back; choosing atoms from the lists goes back to a plane through those atoms. The log shows the angle against the atom plane. Plain click / Shift+click still add / remove atoms of the model in the simulated refinement.

## Profile

**Shift + drag** on the map from atom to atom draws a line; the ends snap to atoms (also to atoms of the molecule beside the section, which is shown while Shift is held). A strip below the map shows the electron density along the line (ρ(Fo), ρ(Fc) and Δρ = ρ(Fo) − ρ(Fc), in e/Å³ over Å, with the atoms near the line marked; relative to the cell mean unless an absolute map type with F000 is chosen). Peak heights depend on Z, ADP and resolution, they are not electron counts. The profile is cut from the 3D grid, so a bond that leaves the plane works too. Extension continues the line beyond both atoms. Hovering the strip moves a marker on the line in the map and vice versa; the × at the top right (or *Clear line*) removes it. The log lists ρ(Fo), ρ(Fc) and Δρ at the atom centres. The line goes with the section: if you change *Offset* (Fourier) or *Position* (Patterson), the line moves by the same amount, so the profile shows the same line in the neighbouring plane (the strip says how far it has moved). The curves follow the map type (absolute types include F000) and the partial model of the simulated refinement.

The Patterson tab has the same profile. The line snaps to the vector rings and the origin peaks; the strip shows P(Fo²), P(Fc²) and ΔP = P(Fo²) − P(Fc²) in the unit chosen in the Map card, and along a cell axis a second axis in fractional coordinates (read a coordinate off a Harker line). *Clip origin peak* leaves the origin peak out of the height axis.

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

## Resolution limit

The resolution sliders (Map card of both tabs) start at the d<sub>min</sub> of the data: the left end means all reflections ("all data (d ≥ 0.80 Å)"), moving to the right cuts reflections with a smaller d-spacing. The slider can only make the resolution worse, never better than the data, and every position has an effect.

## Touch / tablet

On touch and pen input a small bar appears at the bottom (also with `?touch=1`; `?touch=0` hides it). It replaces the modifier keys: *Pan* (drag pans, tap adds a model atom), *Turn* (drag turns the section, tap picks atoms; Fourier tab only), *Profile* (drag from atom to atom draws the line, tap removes a model atom). Two fingers always pinch-zoom and pan. A tap shows the read-out at that point and stays until the next tap; tapping a control with a hint shows the hint. *Reset view* is the double-click. Mouse and keyboard work exactly as before.

## Full screen

The button at the top right of a map (Fourier and Patterson tab) lets the map fill the whole browser window (not the browser's own full-screen mode, so it works the same on desktop and tablet). The sidebar and the log are hidden; the section, zoom, tilt and profile line stay as they are, and clicking, panning and the simulated refinement keep working. A profile line is shown as an overlay in the lower left corner of the map; *Hide profile* at the top right switches it off. **Esc** or the button at the top right goes back. On touch input the button is in the touch bar (*Full screen*), the mode bar stays available, and the exit button is at the top right.
