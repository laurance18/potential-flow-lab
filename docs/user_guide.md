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
| Toolbar (left) | At the top, the **click-to-place** switch (crosshairs). Below it one icon per element type. Hover for its name, click to add one. Below the source-panel body, the pencil that draws one, then the vortex-sheet airfoil. At the bottom, the dye rake tool. |
| Plot (centre) | The flow in the (x, y) plane: streamlines, arrows, element markers. |
| Menu bar (top) | **File** (including exports), **Edit** (undo / redo), **Examples**, **Elements**, **View**, **Verify** and **Help** menus. |
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
  **Click** one to show its dividing streamlines in bold yellow, and click it again to hide them.
  Every streamline through the point is drawn (four at a simple stagnation point, six where two
  have merged), inside a body as well as outside it, up to the next stagnation point, a singular
  point or the edge of the view, with an arrow on each. Several points can be selected at once.
  The lines follow the point through drags, sliders, sweeps and pans, and go away with it if it
  disappears. The selection is not saved with the scene and is not part of the undo history.
- Regions outside a sector flow's wedge are darkened: the flow is not defined there.
- The right end of the menu bar shows the **frame rate** (green at 55 FPS and above, yellow from
  30, red below), the frame time, the **GPU load** of the whole graphics card as Task Manager shows
  it, with this program's share in brackets, and the video memory the program uses. VSync caps the
  frame rate at the monitor's refresh rate. **View → Performance stats** (Ctrl+Shift+F) hides it.

## Adding flow elements

Click a toolbar icon. A dialog opens with default values; change them and press **Apply** (or
**Cancel** / Esc). New elements are placed at the centre of the view. Angles are in degrees,
measured counter-clockwise from +x. If a value is invalid, Apply is disabled and the dialog says
why; if the result is legal but noteworthy (for example two singularities at the same point) a
yellow warning appears before you apply.

**Click to place.** The crosshairs button at the top of the toolbar switches it on or off (off at
start). While it is on, choosing an element shows its marker under the mouse with its coordinates;
click on the plot where it should go (hold Shift to snap to a grid) and the dialog opens with that
position filled in. The position can still be typed over there, and every other parameter is set
as usual. Esc, or clicking the same icon again, cancels. A uniform flow has no position, so its
dialog always opens at once.

| Element | Parameters | Sign convention |
|---|---|---|
| Uniform flow | speed U∞ ≥ 0, direction α | u = U∞ cos α, v = U∞ sin α |
| Source / sink | strength m = Λ, position (x₀, y₀) | Λ > 0 source, Λ < 0 sink (Λ = flow rate per unit depth) |
| Point vortex | circulation Γ, position | Γ > 0 turns **clockwise** (Anderson) |
| Doublet | strength μ, position, axis angle β | F = μ e^(iβ)/(z − z₀); β = α with a uniform flow at α gives a cylinder of radius √(μ/U∞) |
| Sector (wedge) flow | strength V, vertex, wall angles θ₁, θ₂ | 0 < θ₂ − θ₁ ≤ 360°, n = 180°/(θ₂ − θ₁) |
| Symmetry wall | a point on the wall, angle φ | infinite straight wall with zero normal velocity |
| Source-panel body | shape, semi-axes a and b, centre, rotation θ, panels N | a solid ellipse or rectangle, solved by the source panel method (see *Source-panel bodies*) |
| Vortex-sheet airfoil | camber line, chord c, middle of the chord, rotation θ (α = −θ), panels N | a thin airfoil: a vortex sheet γ > 0 clockwise on the camber line, with the Kutta condition (see *Vortex-sheet airfoils*) |

To change or remove an element later, open the **Elements** menu: every element is listed with
its parameters; hover one and choose **Edit...** or **Delete**. **Delete all** (or **File → Clear
All**) removes every element but keeps the view.

### Undo and redo

**Edit → Undo** (Ctrl+Z) takes back the last change to the elements or dye rakes: an element
added, edited, moved or deleted, a parameter slider moved, a sweep played, a rake drawn, moved or
removed, Clear All. **Edit → Redo** (Ctrl+Y or Ctrl+Shift+Z) puts it back. A whole drag or slider
move is one step, however long it lasts, and so is a sweep from Play to Stop (undo stops a playing
sweep). Making a new change after an undo ends the redo list. The view and the View-menu settings
are not part of the history, so undo never moves the plot. The last 200 steps are kept; opening a
file or example, or File → New, starts a new history.

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
  complex velocity, so they sit exactly where W = 0. Switching them off also hides any dividing
  streamlines.
- **Tracer particles** and **Particle trails**: animated particles coloured by speed (blue slow,
  red fast). They are a visual aid only; the streamlines are the result.

**Speed.** The work runs on every core of the processor:
- the tracer particles and the dye;
- the |V| and Cp contours, and the search for stagnation points;
- while you drag, slide or sweep, the streamlines through the previous frame's seeds.

Each result is exactly what one core would compute, only sooner. Near source-panel bodies and
thin airfoils, each group of 16 panels sums its velocity as a series of its own wherever a point is
far enough from it. That gives the same values to rounding as summing the panels one by one, about
twice as fast near a 100-panel body. Animations stay smooth with several panel bodies, and a drag
redraws the streamlines several times a second even around an airfoil next to a body.

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
points, and the source-panel bodies (see *Source-panel bodies*), which come first. It needs a
uniform flow.

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
  | `kind`, `index` | `stagnation` (one row per stagnation point), `body` (one row per body found) or `airfoil` (one row per vortex-sheet airfoil) |
  | `x`, `y` | the stagnation point, the body's centroid (a half-body: its stagnation point), or the middle of the airfoil's chord |
  | `closed` | 1 for a closed body, 0 for a half-body |
  | `chord`, `circulation` | of a closed body or an airfoil |
  | `lift_q`, `drag_q`, `cl`, `cd` | L′/q∞, D′/q∞ and the coefficients from the surface Cp (an airfoil: from ΔCp along the chord, no drag) |
  | `cp_min` | the lowest Cp on the surface (an airfoil: on either side) |
  | `lift_q_exact` | 2Γ/V∞ (Kutta–Joukowski) for an isolated closed body, and for every airfoil |
  | `alpha_deg`, `cm_c4` | an airfoil's angle of attack in degrees and its c_m about the quarter chord |
  | `cl_thin`, `cm_c4_thin` | thin-airfoil theory's c_l and c_m,c/4 at that angle |

  Columns that do not apply are `NaN`. Steps where the value makes the element invalid are left
  out. To plot lift against Γ in Python:
  `d = pandas.read_csv("flow_sweep.csv"); b = d[d.kind == "body"]; plot(b.value, b.lift_q)`.
  For an airfoil's c_l–α line, sweep its **α** and plot
  `a = d[d.kind == "airfoil"]; plot(a.alpha_deg, a.lift_q_exact / a.chord, a.alpha_deg, a.cl_thin)`.

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
divided by the largest speed there. It should be around 10⁻¹⁶.

Two parallel walls make a channel. It needs infinite rows of images, which the program sums in
closed form, so a channel is exact too (see the *Source in a channel* example). Half of a source's
flow leaves each way, at speed Λ/2h far down a channel of width h. Put elements between the walls:
one placed outside has an image inside, which acts as a real element there, and the status bar warns
about it. Other wall layouts with no closed image system, such as three parallel walls or walls not
meeting at 180°/k, are truncated at 32 images per element, and the status bar warns that the wall
condition is then approximate.

## Source-panel bodies

Superposition only makes the bodies that its elements happen to produce: a cylinder, a Rankine
oval, a half-body. The **source panel method** (Anderson, §3.17) works the other way round: you
give the shape, and the program finds the flow around it. Add one with the **Source-panel body**
tool at the bottom of the element icons (below a separator: it is solved numerically, not
written down in closed form).

| Parameter | Meaning |
|---|---|
| Shape | **Ellipse** (a circle when a = b), **Rectangle**, **NACA 4-digit airfoil** or an **outline** you drew or imported (see below) |
| a, b | Half the size along and across the body's own x axis: the semi-axes of the ellipse, the half-width and half-height of the rectangle or the outline, half the chord and half the thickness of the airfoil |
| x₀, y₀ | Centre of the body (drag it there, or grab it anywhere inside it) |
| θ | Turns the body counter-clockwise about its centre |
| N | Number of panels, 4 to 400 (a rectangle or an airfoil needs an even number) |

How it works:

- The surface is cut into N straight panels, numbered clockwise. Each one is a source sheet of
  its own constant strength λ (volume flow per unit length).
- The strengths are not given: they are solved for, so that at the middle of every panel (its
  **control point**) no fluid crosses the surface. The flow there is the free stream, every other
  element of the scene (and its wall images) and all the panels, of every panel body at once, so
  bodies feel each other.
- The surface speed is the tangential velocity at the control points, and Cp follows from
  Bernoulli. The surface-pressure panel plots it at the control points (one marker per panel);
  for an ellipse or circle alone in a uniform stream it also draws the exact curve, from the
  conformal map of the cylinder flow, and prints the largest difference.
- Inside the body the panels make a flow of their own that means nothing physically, so the
  inside is drawn as a solid and left out of the flow: streamlines stop at the surface and
  particles never go in. The status bar says so when the mouse is inside.

What to expect:

- **A circle is exact at the control points**, for any N: Anderson's 8-panel cylinder already
  gives Cp = 1 − 4 sin²θ there. Between the control points and away from the body it is not.
- **An ellipse converges**: the surface Cp error falls about 4 times each time N doubles (second
  order), the velocity away from the body only about 2 times (first order). Sweep N in the
  **Parameters** panel to watch it.
- **No lift, no drag.** Source panels carry no circulation, so a body alone in a stream feels no
  force (d'Alembert), at any angle. Put a **vortex inside** the body to give it circulation: inside
  a circle that is the lifting cylinder, and the lift from the panels approaches L′ = ρV∞Γ. Only
  the circulation of what lies inside a body shows outside it; the panels cancel a source or
  doublet inside it, and a warning says so. Nothing in this method fixes how much circulation an
  airfoil gets: that is the Kutta condition, which needs vortex panels.
- Sharp corners (the rectangle) have infinite speed in exact potential flow. The panels are made
  shorter near the corners, but Cp right at a corner stays only an estimate.

**Seeing the source strengths.** **View → Panel strengths (λ, γ)** draws the solution itself on
every panel body: a stick from each control point, normal to its panel, as long as the panel's
strength λ (the longest stick is the largest |λ| of all the panels; a caption gives its value).
It points outwards, in green, for a source (λ > 0, fluid pushed out of the surface) and inwards,
in red, for a sink. The line through the tips is the strength distribution as a polar plot around
the body. Hover a stick for the panel's number, λ, λ/V∞, λ/2πV∞, its Cp and the volume flow λs
it lets out. With **Surface pressure** on, the side panel also plots λ/V∞ along the surface, on the
same axis as Cp, over the exact continuous sheet where it is known (an ellipse or circle alone in a
uniform stream). The exact sheet leaves a uniform flow inside the body, so λ = −v_inside · n; on a
circle λ = −2V∞ cos(θ − α), 2V∞ on the point facing the stream.

The strengths converge more slowly than Cp, at first order: Anderson's 8 panels give 2.366 V∞ on
the panel facing the stream, 18 % above the sheet's 2V∞, while their Cp at the control points is
already exact. Sweep N and watch the sticks settle onto the sheet. On the rectangle the short
corner panels carry the largest strengths, where the flow has to turn round the corner.

The status bar shows the boundary condition as it is solved: the largest |v · n| at the control
points over the speed there (about 10⁻¹⁵), and in its tooltip Σλs, the net source strength of each
body, which is zero for the exact solution of a closed body (and zero by symmetry for an ellipse or
rectangle alone in a stream).

### Shapes of your own: drawing, airfoils and imported coordinates

Besides the ellipse and the rectangle, a source-panel body can be:

- **A NACA 4-digit airfoil.** Choose *NACA 4-digit airfoil* as the shape and give the camber m
  (% of the chord), its position p (tenths of the chord), the chord c and the thickness t/c; the
  dialog shows the designation (NACA 2412, ...). The airfoil lies along the body's x axis, leading
  edge in front; with a stream along +x its angle of attack is α = −θ. The N panels (an even number)
  are cosine-spaced along the chord, short at both edges.
- **A drawn outline.** Pick the **pencil tool** (just below the source-panel body icon) and click
  points around the shape you want, or hold the button and drag to sketch it. Click the first
  point again (or press Enter, or double-click) to close it; Backspace takes back the last point
  and Esc gives up. The dialog then opens with the new body where you drew it.
- **Imported airfoil coordinates.** **File → Import airfoil coordinates (.dat)...** reads an x y
  coordinate file in the Selig or Lednicer layout, the two used by the UIUC airfoil database. The
  airfoil gets a chord of 2 and its own thickness, 120 panels, and the name from the file.

For a drawn or imported outline the points only give the shape; the N panels are laid along it
independently of them:

- **Smooth outline** (the default): a smooth curve (a cubic spline) through the points, except at
  sharp corners, turning by more than 100° (a trailing edge), which stay sharp. A dozen points
  already give a smooth blob; points on an ellipse are joined into the ellipse to within 3·10⁻⁴
  with 24 points and 2·10⁻⁵ with 48.
- **Straight outline** (untick *Smooth outline*): straight lines from point to point, with a corner,
  and so a panel vertex, at every point where the outline turns. That needs at least as many panels
  as corners.
- Every corner is a panel vertex. Between corners the panels are shared out by length, a little
  shorter where the outline curves, and cosine-spaced towards corners sharper than 30°.
- a and b are the half-width and half-height of the outline: it is stretched to fit 2a × 2b.
  *Keep the outline's proportions* sets b back so that it is not distorted.
- An outline must not cross itself (the dialog says so if it does). Concave shapes are fine: a body
  that no point inside it sees whole (a crescent, a C) sums the panels of its complex potential one
  by one instead of by vertex, with the same result. Its φ and ψ exports then jump across lines
  running in −x from the vertices; the velocity does not.

What to expect from these shapes: Cp converges at first order on them (the error halves when N
doubles), as constant-strength panels do on a general spacing. The ellipse preset's equal steps
are a special case, with second order. An airfoil at incidence shows the limit of source panels
most clearly. Without circulation it has no lift, and the flow has to turn round the sharp trailing
edge, where its exact speed is infinite: the panels there need huge strengths (the sticks of
**Panel strengths** are cut short), and the force from Cp comes to d'Alembert's zero only
slowly. Putting the Kutta condition on that trailing edge is what a vortex sheet does (see *Vortex-sheet airfoils*).

Limits: panel bodies have no images in symmetry walls (the walls are then not streamlines near the
body, and a warning says so), and overlapping bodies have no meaningful solution (also warned).
Stagnation points on a panel body are found where the tangential speed changes sign between two
control points.

## Vortex-sheet airfoils (thin-airfoil theory)

A source-panel airfoil gets no lift: nothing in it sets the circulation. Classical **thin-airfoil
theory** (Anderson, §4.7 to 4.8) sets it with the Kutta condition. It replaces the airfoil by its
**camber line** and puts a vortex sheet γ(s) on it, with γ = 0 at the trailing edge so the flow
leaves it smoothly. The **Vortex-sheet airfoil** tool (the last of the panel-method icons, a curved
line with sticks on it) does exactly that, numerically.

| Parameter | Meaning |
|---|---|
| Camber line | **Flat plate**, **NACA 4-digit mean line** (camber m % of the chord at p tenths of the chord, as for the NACA body) or **circular arc** (camber h/c, its height at mid-chord, up to 25 %) |
| Chord c | from the leading to the trailing edge |
| x₀, y₀ | the **middle of the chord** (drag it there, or grab the airfoil anywhere along its line) |
| θ | turns the airfoil counter-clockwise about the middle of its chord; in a stream along +x the angle of attack is **α = −θ**, as for the NACA source-panel body |
| N | number of vortex panels along the camber line, 4 to 400 |

The camber line has no thickness: the flow goes past both sides of it. In the Parameters panel the
airfoil has a slider for **α** as well as θ (the same rotation, read the other way), so a sweep of α
gives its c_l–α line (see *Live parameters and sweeps*).

How it works:

- The camber line is cut into N straight panels, cosine-spaced along the chord
  (x = (c/2)(1 − cos θ₀) at equal steps of θ₀), so they are short at both edges.
- γ varies linearly along each panel between its values at the vertices, so the sheet strength is
  continuous. The trailing-edge value is γ = 0: the **Kutta condition**. The other N values are
  solved so that no fluid crosses the camber line at the middle of each panel.
- The airfoils are solved **in one system with the source-panel bodies**, against everything else
  in the scene (wall images included), so airfoils and bodies feel each other.
- Why these panels: constant-strength vortex panels collocated at their midpoints have no Kutta
  condition (on a flat plate they converge to the flow without circulation). Point vortices, as in
  the lumped-vortex 1/4–3/4 rule, make the streamlines and particles near the line jump. Linear
  panels keep the velocity finite everywhere along the line except at the leading edge, where the
  exact γ is infinite.
- The panels apply the boundary condition on the curved line itself, not on its chord, and at the
  real angle. So they converge to the exact potential flow past the infinitely thin line: for a
  flat plate c_l = 2π sin α, for a circular arc the conformal-map result. Thin-airfoil theory is
  the linearization of that (2πα for the plate). The difference is of order α² and (dz/dx)²: about
  0.1 % for a NACA 2412 at 4°.

What you see:

- **The plot**: the camber line as a thin solid line. The flow passes smoothly off the trailing
  edge, and the front stagnation point sits just under the leading edge at positive α. The
  stagnation points on the sheet are where the speed along one side changes sign between two
  control points. Streamlines stop rather than cut through the line.
- **View → Panel strengths (λ, γ)** draws γ at every vertex as a stick normal to the camber line:
  above it (amber) for γ > 0, which is clockwise and makes the flow over the top faster, and below
  it (cyan) for γ < 0. The line through the tips is the loading. The leading-edge stick is cut
  short (the exact γ is infinite there), and the trailing-edge one is zero. Hover a stick for γ,
  ΔCp = 2γ/V∞ and thin-airfoil theory's value there.
- **View → Surface pressure** opens a section for each airfoil in the side panel:
  - ΔCp(x/c) = 2γ/V∞ along the chord, the panels' values at their vertices over the
    thin-airfoil curve 4[A₀(1 + cos θ)/sin θ + Σ Aₙ sin nθ]. Hover a point to mark it on the
    airfoil.
  - The theory's Fourier coefficients A₀ to A₃: A₀ = α − (1/π)∫ dz/dx dθ₀ and
    Aₙ = (2/π)∫ dz/dx cos nθ₀ dθ₀. They are exact for the flat plate and the NACA mean line, and
    summed to full precision for the arc.
  - The panels beside the theory:
    - Γ;
    - c_l from Kutta–Joukowski, 2Γ/(V∞c), and from ∫ΔCp d(x/c) along the chord, both against
      2π(A₀ + A₁/2);
    - c_m about the quarter chord against (π/4)(A₂ − A₁);
    - α_L0 against −(1/π)∫ dz/dx (cos θ₀ − 1) dθ₀;
    - the lift slope against 2π;
    - γ at the trailing edge.
  - The two lift values agree exactly on a flat plate. With camber the ΔCp one is slightly lower,
    because γ is per unit length of the curved line while ΔCp is integrated along the chord
    (dx = ds cos ε). The theory does not tell them apart.
  - The panels' α_L0 and lift slope come from solving for unit streams along x and y as well. They
    are shown only when the rest of the flow is uniform (other panel bodies and airfoils are
    allowed).
- **Status bar**: the boundary condition at the airfoil's control points (about 10⁻¹⁵), and in its
  tooltip Γ.

Signs: γ and Γ are positive clockwise, like a vortex's Γ, so Γ > 0 lifts in a stream along +x
(L′ = ρV∞Γ). ΔCp = Cp_lower − Cp_upper, c_m > 0 nose up. φ is multivalued around an airfoil
with circulation. The field export puts its jump of Γ on the ray running on from the trailing edge
along the chord. Across the sheet itself φ jumps as well, as it does across any vortex sheet.

Warnings: no uniform flow (so no stream past the trailing edge), an angle of attack beyond 15° (a
real airfoil would have stalled), an airfoil running into a source-panel body or crossing another
airfoil, an element lying on the camber line, and walls. Like panel bodies, airfoils have no
images in walls.

Not included: the airfoils' surfaces are not in the surface-Cp CSV export, and thick airfoils
with vortex panels on their surface are not modelled.

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
you. Each of the 61 checks runs the program's own code and compares the result with the exact one,
written out separately, mostly in polar form, without using the program's equations. The checks are
in ten groups:

| Group | What is checked |
|---|---|
| Foundations | div V = 0 and curl V = 0 (Laplace), W = dF/dz for every element, flux and circulation around closed curves, the sum of every element, finite velocity at singular points, turned and mirrored scenes, φ and ψ as integrals of V, and Bernoulli against Euler's equation |
| Elementary flows | Uniform flow, source and sink, vortex, doublet, corner flow and plate edge, at every node of a 201 × 117 grid |
| Superposition and bodies | Cylinder with and without circulation, its stagnation points on and off the surface, dividing streamlines, the Rankine oval and half-body (ψ along the traced surface), a source near a wall |
| Symmetry walls | No flow through an inclined wall, two perpendicular walls or a channel, the channel's closed form against its series of images, a 60° corner against its conformal map, and a vortex in a channel |
| Pressure and forces | Cp on the cylinder and the Rankine bodies, d'Alembert's paradox, Kutta–Joukowski lift, forces from the far field (Blasius), and forces along a parameter sweep |
| Streamlines and equipotentials | Straight lines in a uniform flow, rays from a source or sink, circles around a vortex in the right direction, ψ constant along streamlines, even spacing, arrows along the flow, and equipotentials at 90° to the streamlines |
| Particles and dye | Tracer particles against the exact orbit around a vortex and against ψ past a cylinder; dye against ψ and the spacing V∞Δt |
| Panel methods | Source panels on a cylinder (exact at the control points, and Anderson's 8-panel strengths), on an ellipse (second-order convergence of Cp), the panel strengths against the exact source sheet (first order), the panel sums near and far, the field around an ellipse against its conformal map, a lifting cylinder from a vortex inside panels, no flow through the panels, drawn outlines (the spline through points on an ellipse, and its flow), NACA 4-digit geometry and Σλs, and a concave body against the panels summed directly |
| Thin airfoils | A flat plate's c_l → 2π sin α (second order) and its γ(x) = 2V∞ sin α √((c − x)/x) with γ(TE) = 0; the NACA 2412 mean line's α_L0 and c_m,c/4 against thin-airfoil theory (itself checked by quadrature); a circular arc against the exact Joukowski solution, c_l and the whole field; ∮V·dl = Γ and the Blasius lift ρV∞Γ with no drag; the panel sums against the Biot–Savart integral; and no flow through the camber line beside a source-panel body |
| Saving and reopening | A scene with every kind of element saved to a file and reopened: every parameter, the rakes, the plot settings and the streamlines come back |

The checks run one after another, and each turns into a green tick when it passes (red cross if
not). The bar beside it shows −log₁₀ of the error, about how many digits agree with theory, and
the white tick on the bar marks the tolerance. Where both sides are exact formulas, only rounding
separates them, so the tolerance is 10⁻¹². Checks that trace a surface, integrate a force,
differentiate numerically or move particles are only as accurate as that step, and their tolerances
are looser.

Click a check to open it. You see what theory says, how the program was tested, the program's
numbers next to the exact ones, the error and the tolerance. **Show in plot** loads that flow into
the plot, so you can look at it and change it. **Copy report** copies every result as plain text,
for example for a lab report. The same checks run in the unit tests (`tests/test_verification.cpp`).

## Help menu

- **User guide** (F1): this guide, built into the exe. The chapters are listed on the left; the
  search box keeps only the chapters that mention a word. Esc closes it.
- **Sign conventions...**: the coordinate and angle conventions, the complex potential and
  velocity of every element with the meaning of each sign (Λ, Γ, μ, V, θ₁, θ₂, φ), how the wall
  images are formed, and Cp, lift and drag.
- **Keyboard shortcuts...**: every key and mouse gesture, grouped by what it acts on.
- **About Potential Flow Lab...**: the version, the libraries built in, the OpenGL version and
  the graphics card in use, with a Copy button for problem reports.

## Troubleshooting

- **"OpenGL 3.3 required" at start-up**: update the graphics driver, and run the program on the
  computer itself, not through Remote Desktop or a virtual machine.
- **Greek letters show as "?"**: the program uses the Segoe UI font from Windows; if it is
  missing, it falls back to a basic font without Greek letters.
- **(!) warnings in the status bar**: hover them for the full text. They describe legal but easily
  misread configurations (coincident singularities, elements on a wall, truncated images).
