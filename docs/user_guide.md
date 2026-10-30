# Potential Flow Lab: User Guide

Potential Flow Lab builds two-dimensional, incompressible, irrotational flows by adding up
elementary solutions, and redraws the streamlines every time the configuration changes.

## Getting started

Run `PotentialFlowLab.exe`. It needs no installation; Windows 10/11 with an OpenGL 3.3 graphics
driver is enough. To start from an example, pick one from the **Examples** menu (hover an entry
for a short description); it replaces the current scene, view and View-menu settings in one
click. The same scenes are in `examples/` as files, and any scene file can also be opened with
**File → Open...** or dragged onto the exe.

The window has four parts:

| Part | What it does |
|---|---|
| Toolbar (left) | One icon per element type. Hover for its name, click to add one. Below them, the dye rake tool. |
| Plot (centre) | The flow in the (x, y) plane: streamlines, arrows, element markers. |
| Menu bar (top) | **File** (including exports), **Examples**, **Elements**, **View** and **Verify** menus. |
| Status bar (bottom) | Velocity under the mouse, streamline count, wall check, warnings. |

## The flow view

- **Zoom** with the mouse wheel and **pan** by dragging. Both axes always use the same scale, so
  circles look like circles. **View → Reset view** returns to −4 < x < 4, −2.5 < y < 2.5.
- Streamlines are drawn evenly spaced, with arrows in the flow direction. They are recomputed
  after every change and whenever the view moves.
- Hover the plot to read x, y, the velocity (u, v), the speed |V| and, when a uniform flow is
  present, the pressure coefficient Cp = 1 − |V|²/U∞² in the status bar.
- Hover an element marker to see its parameters; **double-click** it to edit it.
- **Drag** a marker to move the element. Everything follows live while you drag: streamlines,
  equipotentials, contours, stagnation points, body surfaces and the Cp panel. The tracer particles
  keep flowing through the changing field. Positions are rounded to about one pixel; hold **Shift**
  to snap to a coarser grid (0.1 at the default zoom). **Esc** during a drag puts the element back.
  A wall is moved parallel to itself (grab it anywhere along its line). A uniform flow has no
  position, so it cannot be dragged; double-click its arrow to change it.
- Stagnation points (V = 0) are marked with red diamonds; hover one to read its coordinates.
- Regions outside a sector flow's wedge are darkened: the flow is not defined there.

## Adding flow elements

Click a toolbar icon. A dialog opens with default values; change them and press **Apply** (or
**Cancel** / Esc). New elements are placed at the centre of the view. Angles are in degrees,
measured counter-clockwise from +x. If a value is invalid, Apply is disabled and the dialog says
why; if the result is legal but noteworthy (for example two singularities at the same point) a
yellow warning appears before you apply.

| Element | Parameters | Sign convention |
|---|---|---|
| Uniform flow | speed U∞ ≥ 0, direction α | u = U∞ cos α, v = U∞ sin α |
| Source / sink | strength m = Λ, position (x₀, y₀) | Λ > 0 source, Λ < 0 sink (Λ = flow rate per unit depth) |
| Point vortex | circulation Γ, position | Γ > 0 turns **clockwise** (Anderson) |
| Doublet | strength μ, position, axis angle β | F = μ e^(iβ)/(z − z₀); β = α with a uniform flow at α gives a cylinder of radius √(μ/U∞) |
| Sector (wedge) flow | strength V, vertex, wall angles θ₁, θ₂ | 0 < θ₂ − θ₁ ≤ 360°, n = 180°/(θ₂ − θ₁) |
| Symmetry wall | a point on the wall, angle φ | infinite straight wall with zero normal velocity |

To change or remove an element later, open the **Elements** menu: every element is listed with
its parameters; hover one and choose **Edit...** or **Delete**. **Delete all** (or **File → Clear
All**) removes every element but keeps the view.

## Streamlines, equipotentials, contours and tracer particles

The **View** menu holds the plot settings, which are saved with the scene:

- **Streamline density**: number of streamlines across the height of the view.
- **Direction arrows** and **Arrow spacing** (screen pixels between arrows on a line).
- **Show wall images**: draw the mirrored elements that enforce the symmetry walls (faded).
- **Equipotential lines** (orange): lines of constant velocity potential φ. They cross every
  streamline at right angles and use the same spacing, so in regions of typical speed the two
  families form roughly square cells.
- **Background**: filled bands under the streamlines, with a colour scale in the top-right corner.
  - *Speed contours |V|*: velocity magnitude from 0 to twice the typical (median) speed in the
    view. Faster regions, such as those next to a source or vortex, use the top band.
  - *Pressure coefficient Cp*: Cp = 1 − |V|²/V∞² in bands of 0.25 (wider for very low Cp). Red
    is above and blue below the free-stream pressure; Cp = 0 is neutral grey. The scale runs from 1
    (stagnation) down to the Cp of twice the typical speed. Lower values, as next to a singular
    point, use the bottom band. Cp needs a uniform flow; without one, the scale box says so.
- **Surface pressure**: see the next section.
- **Stagnation points**: the diamond markers. They are found by Newton's method on the analytic
  complex velocity, so they sit exactly where W = 0.
- **Tracer particles** and **Particle trails**: animated particles coloured by speed (blue slow,
  red fast). They are a visual aid only; the streamlines are the result.

## Dye rakes

A dye rake is a line of dye emitters, as in a water or smoke tunnel. Click the rake tool at the
bottom of the toolbar, then drag across the plot to draw the rake (a plain click places a vertical
one). Each emitter releases dye continuously, and the dye is carried by the flow. The resulting
filaments are streaklines. In a steady flow these coincide with the streamlines through the
emitters, so the dye shows the same curves as the streamlines, but you can watch it travel.

The dye is drawn in pulses of half a second, alternately bright and dim. Each pulse holds fluid
that left the emitter during the same time span, so pulses stretch where the flow is fast and
bunch up where it is slow, for example ahead of a stagnation point. Dye fades out after 15 seconds
of animation and disappears where it leaves the view or runs into a singular point.

- Drag a rake's line to move it, or drag one of its square end handles to turn or stretch it.
- Right-click a rake to change its number of emitters (1–64) or delete it. **Delete**
  removes the rake under the mouse.
- **Elements → Delete all dye rakes** removes them all; **File → Clear All** removes rakes along
  with the elements.
- Rakes change nothing in the flow; they are only a visualization. They are saved with the scene.

## Surface pressure and forces

**View → Surface pressure** finds the bodies in the flow: the closed bodies (cylinder, lifting
cylinder, Rankine oval) and half-bodies bounded by the dividing streamlines through stagnation
points. It needs a uniform flow.

- In the plot, each body surface is coloured by Cp on the same scale as the Cp background. Arrows
  normal to the surface point towards it where Cp > 0 (pressure) and away from it where Cp < 0
  (suction). Arrow length is proportional to |Cp|. With several bodies, each is numbered.
- A panel on the right plots Cp along the surface. The x axis can be the angle θ around the body
  (from the free-stream direction, counter-clockwise), x/c along the free stream, or the arc
  length s from the front stagnation point. Half-bodies are open, so they are only shown against
  s. The Cp axis is inverted (suction up), as usual in aerodynamics. For a circular cylinder the
  exact curve is drawn as a light band under the computed one. Hover the plot to read a point;
  the plot on the left marks it on the body.
- Below the plot, for every closed body: the chord, the circulation enclosed, the lowest Cp, and the
  drag and lift per unit span from integrating Cp around the surface (D′/q∞, L′/q∞, c_l). When the
  body is alone in a uniform stream, the exact values are shown next to them: zero drag
  (d'Alembert) and L′/q∞ = 2Γ/V∞ (Kutta–Joukowski).

A body is only found if its dividing streamlines close. Next to other bodies the flow is no longer
symmetric and they miss each other slightly. Gaps up to 10 % of the distance between the
stagnation points are closed with a short segment, and the panel says so. Some flows have no closed
body at all: for a source and sink at incidence, part of the source's fluid escapes downstream.

## Live parameters and sweeps

**View → Parameters** (Ctrl+P) opens a panel on the right with a slider for every parameter of
every element: strengths, speeds, positions and angles (in degrees). With the Cp panel open too,
the parameters take the top part of the column.

- Drag a slider and the flow follows live, as when dragging an element on the plot: streamlines,
  contours, stagnation points, bodies and the Cp panel update while the particles keep flowing.
- **Ctrl+click** a slider to type an exact value. A slider's range is chosen when it first shows:
  ±2|value| rounded up to 1, 2 or 5 × 10ⁿ (at least ±2) for strengths, 0 to that for speeds, ±180°
  for angles and the current view for positions. Typing a value outside the range widens it.
- Values that would make an element invalid (a negative speed, a sector with θ₂ ≤ θ₁) are not
  applied.

The **Sweep** section steps one parameter from **From** to **To**. Choose it in the list, or
right-click its slider and pick *Sweep this parameter*; its label turns blue.

- **Play** animates it back and forth over **Duration** seconds (or From → To and start again,
  with *Back and forth* off). Stop leaves the parameter where it is. A classic one: sweep Γ of the
  lifting cylinder from 0 to 4π and watch the two stagnation points slide down and merge.
- **Export CSV...** evaluates the flow at **Steps** evenly spaced values (both ends included) over
  the current view, a few per frame with a progress bar, and writes one file:

  | Column | Meaning |
  |---|---|
  | `step`, `value` | step number (from 1) and the parameter value (angles in radians) |
  | `kind`, `index` | `stagnation` (one row per stagnation point) or `body` (one row per body found) |
  | `x`, `y` | the stagnation point, or the body's centroid (a half-body: its stagnation point) |
  | `closed` | 1 for a closed body, 0 for a half-body |
  | `chord`, `circulation` | of a closed body |
  | `lift_q`, `drag_q`, `cl`, `cd` | L′/q∞, D′/q∞ and the coefficients from the surface Cp |
  | `cp_min` | the lowest Cp on the surface |
  | `lift_q_exact` | 2Γ/V∞ (Kutta–Joukowski) for an isolated closed body |

  Columns that do not apply are `NaN`. Steps where the value makes the element invalid are left
  out. To plot lift against Γ in Python:
  `d = pandas.read_csv("flow_sweep.csv"); b = d[d.kind == "body"]; plot(b.value, b.lift_q)`.

## Recording animations

**File → Record animation (MP4 / GIF)...** records the plot, axes included, as a video:

| Option | |
|---|---|
| Format | **MP4** (H.264): small, smooth colours, plays everywhere (PowerPoint, browsers, phones). **Animated GIF**: loops forever and goes wherever an image can (web pages, chat, Markdown), with at most 255 colours per frame and larger files. |
| Length | seconds of video (Ctrl+click to type up to 10 minutes) |
| Frame rate | MP4 24–60 fps, GIF 10–50 fps |
| Size | 100 %, 75 % or 50 % of the plot on screen. Half-size GIFs are about a quarter of the file size. |
| Start from fresh particles and dye | respawn the particles and clear the dye first, so the video shows the dye leaving the rakes |
| Sweep ... | play the parameter sweep once, From to To, over the length of the video |

After choosing the file name, recording starts; the menu bar shows a red dot, the time recorded and
a **Stop** button (File → Stop recording does the same and keeps what was recorded). You can keep
working: drag elements, move sliders or pan, and the video shows it.

The animation advances exactly one video frame (1/fps) for every frame recorded. The video
therefore plays at the right speed even when the screen is slower than the video's frame rate, or
when encoding makes it stutter. On a fast screen the frames drawn between two video frames simply
hold the animation still. Resizing the window stops the recording, and the video so far is kept.

MP4 export uses Windows Media Foundation. Windows "N" editions do not include it unless the
Media Feature Pack is installed; there, record a GIF instead.

## Symmetry walls

A symmetry wall is an infinite straight line that the flow cannot cross. The program enforces it
with the method of images rather than just drawing a line:

- sources, vortices and doublets are mirrored across the wall (vortices with −Γ, doublets with
  axis angle 2φ − β);
- a uniform flow keeps only its component along the wall (with walls that cross, none of it);
- sector flows are not mirrored, so a wall is a streamline of a sector only if it lies along one
  of the sector's own walls (the status bar warns otherwise).

The status bar shows the numerical check: the largest |v · n| sampled along each wall in view,
divided by the largest speed there. It should be around 10⁻¹⁶. Two parallel walls need an
infinite row of images, so the program truncates it and warns that the wall condition is then
approximate.

## Saving and loading configurations

**File → Save** (Ctrl+S) and **Save As...** (Ctrl+Shift+S) write the scene as a JSON file: every
element with its parameters, the walls, the dye rakes, the current view and the View-menu settings. **File →
Open...** (Ctrl+O) restores it exactly. **File → New** (Ctrl+N) starts an empty scene with
default settings. The title bar shows `*` when there are unsaved changes, and New, Open, the
Examples menu and Quit ask before discarding them. An example opens as an unsaved copy, so **Save**
asks for a file name and never changes the built-in example.

Files are plain text and can be edited by hand; angles are stored in degrees. A file with an
invalid element (for example θ₂ ≤ θ₁) is rejected with a message naming the element.

## Exporting images and data

The **File** menu has four exports, besides **Record animation** (see *Recording animations*);
the sweep CSV is exported from the Parameters panel (see *Live parameters and sweeps*). Each one
asks where to save the file.

| Export | What it writes |
|---|---|
| **Export image (PNG)...** | The plot, axes included, exactly as shown (hover highlights are left out). |
| **Export field data (CSV)...** | `x, y, u, v, speed, phi, psi, cp` at the nodes of a 201-column grid over the current view. |
| **Export lines (CSV)...** | `kind, line, x, y`: every point of every streamline and equipotential line shown. |
| **Export surface Cp (CSV)...** | `body, closed, x, y, s, streamwise, theta_deg, cp, cp_exact`: every point of every body surface found (see *Surface pressure and forces*). `streamwise` is x/c for a closed body and the distance downstream of the stagnation point for a half-body; `cp_exact` is the closed-form Cp for a circular cylinder, otherwise `NaN`. |

In the field file, `NaN` marks nodes outside the flow domain or on a singular point, and `cp` is
`NaN` when there is no uniform flow. φ is multivalued around a vortex and ψ around a source. The
file uses the principal branch, so φ jumps by Γ (and ψ by Λ) across the line running in −x from
the element. The files open directly in Excel, MATLAB (`readtable`) or Python
(`numpy.genfromtxt(..., names=True)`).

### From the command line

The same exports can be run without the interface, for example to regenerate figures:

```powershell
PotentialFlowLab.exe examples\cylinder.json --export-image cylinder.png --export-field cylinder.csv --export-lines lines.csv --export-surface surface.csv
```

The window opens briefly, the scene is drawn, the files are written and the program exits. The
exit code is 0 on success and 1 if a file could not be written. `tools/report_figures.py` uses this
to produce every figure and validation number in the report.

Videos too: `--export-video` records an MP4, or a GIF when the name ends in `.gif`, after the other
exports. `--video-seconds` (default 6) and `--video-fps` (default 30 for MP4, 20 for GIF) set its
length and frame rate:

```powershell
PotentialFlowLab.exe examples\lifting_cylinder.json --export-video lifting_cylinder.mp4 --video-seconds 10
```

## Verifying the physics

**Verify → Physics checks...** (F5) tests the program against potential-flow theory, in front of
you. Each of the 22 checks builds a flow with the program's own code and compares it with the exact
result, written out separately in polar form, without using the program's equations. The checks
are in six groups:

| Group | What is checked |
|---|---|
| Foundations | div V = 0 and curl V = 0 (Laplace), W = dF/dz for every element, and the flux and circulation around closed curves |
| Elementary flows | Uniform flow, source and sink, vortex, doublet, corner flow and plate edge, at every node of a 201 × 117 grid |
| Superposition and bodies | Cylinder with and without circulation, its stagnation points, the Rankine oval and half-body (ψ along the traced surface), a source near a wall |
| Symmetry walls | No flow through an inclined wall or through two perpendicular walls |
| Pressure and forces | Cp on the cylinder, d'Alembert's paradox, Kutta–Joukowski lift |
| Streamlines and equipotentials | Vortex streamlines are circles, and equipotentials cross the streamlines at 90° |

The checks run one after another, and each turns into a green tick when it passes (red cross if
not). The bar beside it shows −log₁₀ of the error, about how many digits agree with theory, and
the white tick on the bar marks the tolerance. Where both sides are exact formulas, only rounding
separates them, so the tolerance is 10⁻¹². Checks that trace a surface, integrate a force or
differentiate numerically are only as accurate as that step, and their tolerances are looser.

Click a check to open it. You see what theory says, how the program was tested, the program's
numbers next to the exact ones, the error and the tolerance. **Show in plot** loads that flow into
the plot, so you can look at it and change it. **Copy report** copies every result as plain text,
for example for a lab report. The same checks run in the unit tests (`tests/test_verification.cpp`).

## Troubleshooting

- **"OpenGL 3.3 required" at start-up**: update the graphics driver, and run the program on the
  computer itself, not through Remote Desktop or a virtual machine.
- **Greek letters show as "?"**: the program uses the Segoe UI font from Windows; if it is
  missing, it falls back to a basic font without Greek letters.
- **(!) warnings in the status bar**: hover them for the full text. They describe legal but easily
  misread configurations (coincident singularities, elements on a wall, truncated images).
