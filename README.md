# Potential Flow Lab

An interactive 2D potential-flow simulator in C++20, built as a senior mechanical-engineering
semester project. It superposes elementary flows (uniform flow, source/sink, vortex, doublet,
sector flow) and straight symmetry walls, draws evenly spaced streamlines with direction arrows,
equipotential lines and velocity-magnitude contours, and animates tracer particles on the GPU.

- Left toolbar: one icon per element type (hover for the name); a click opens a parameter dialog
  with Apply / Cancel. Below them, the dye rake tool.
- Main plot: labelled (x, y) axes, scroll to zoom, drag to pan; streamlines are recomputed after
  every add, edit or delete. Hover an element to name it, double-click it to edit it, or **drag it**
  to move it: the streamlines, contours, stagnation points and body surfaces follow live while the
  particles keep flowing.
- **Dye rakes**: drag a line across the plot and its emitters release continuous, pulsed dye
  streaks that the flow carries away (streaklines; in steady flow they trace the streamlines and
  the pulse spacing shows the speed). Rakes can be moved, turned and edited, and are saved with
  the scene.
- **Parameters** panel (View → Parameters, Ctrl+P): a live slider for every element parameter;
  the flow follows as you drag. A **sweep** steps one parameter through a range: Play animates it,
  Export CSV writes the stagnation points and each body's lift, drag, c_l, c_d and lowest Cp at
  every step (with the Kutta–Joukowski lift for comparison).
- **Record animation** (File menu): the animated plot as an MP4 video (H.264, Windows Media
  Foundation) or a looping animated GIF, at a chosen length, frame rate and size. The animation
  advances exactly one video frame per captured frame, so the video plays at the right speed. You
  can drag elements and move sliders while recording, or play the sweep once over the video.
- **Elements** menu: every element with its parameters, each with Edit and Delete.
- **Edit** menu: Undo (Ctrl+Z) and Redo (Ctrl+Y / Ctrl+Shift+Z) of every change to the elements
  and dye rakes. A whole drag, slider move or sweep is one step.
- **Examples** menu: the scenes in `examples/` (cylinder, lifting cylinder, Rankine bodies, corner
  flow, ...) built into the exe; one click loads the elements, view and plot settings.
- **File** menu: New, Open, Save, Save As (JSON scenes, see `examples/`), Clear All, and exports:
  PNG image of the plot, CSV field data (u, v, |V|, φ, ψ, Cp on a grid), CSV streamline and
  equipotential polylines, CSV surface Cp of the bodies, and MP4/GIF recordings. The exports also
  run from the command line.
- **View** menu: streamline density, arrows, wall images, equipotential lines, a |V| or Cp
  background, surface pressure, stagnation point markers, tracer particles.
- **Verify** menu (F5): 22 physics checks run in front of you. Each compares the program with an
  exact result: elementary flows over a whole grid, Laplace's equation, the cylinder and Rankine
  bodies, wall images, surface Cp, d'Alembert and Kutta–Joukowski. Each check shows its error
  against its tolerance and the program's numbers beside theory, and can load its flow into the plot.
- **Surface pressure**: the bodies (cylinder, Rankine oval, half-body, ...) are found from their
  stagnation points; their surfaces are coloured by Cp with pressure arrows, and a side panel plots
  Cp along the surface (exact curve for a cylinder) with the integrated lift and drag compared with
  Kutta–Joukowski and d'Alembert.
- Status bar: velocity and Cp under the mouse, the numerical wall check, and warnings.

See [docs/user_guide.md](docs/user_guide.md) for usage, [docs/report.md](docs/report.md) for the
project report (equations, architecture, validation, examples) and
[docs/report_notes.md](docs/report_notes.md) for the working notes behind it.

## Report figures and validation

```powershell
python tools/report_figures.py     # needs numpy and matplotlib, and a Release build
```

This renders the report figures into `docs/figures/` with the exe's command-line exports and
compares the exported fields with the closed-form solutions (`docs/figures/validation.md`).

The release build is one self-contained `PotentialFlowLab.exe`. It needs no installer, DLLs,
Visual C++ Redistributable or extra folders.

## Prerequisites (Windows 10/11 x64)

| Tool | Notes |
|---|---|
| Visual Studio 2022 with the **Desktop development with C++** workload | Provides MSVC and the Windows SDK |
| CMake ≥ 3.25 | Put it on `PATH` (cmake.org installer or `winget install Kitware.CMake`) |
| Git | Used by CMake to fetch dependencies |
| VS Code + **CMake Tools** + **C/C++** extensions | Optional; recommended in `.vscode/extensions.json` |

All libraries are downloaded and built from source on the first configure, so the first configure
needs internet. Nothing else has to be installed.

## Build, run, test (PowerShell, from the repository root)

```powershell
cmake --preset msvc                     # configure (Visual Studio 2022 generator, x64)
cmake --build --preset debug            # or: cmake --build --preset release
.\build\msvc\bin\Debug\PotentialFlowLab.exe
ctest --preset debug                    # unit tests (doctest), also: ctest --preset release
```

A plain PowerShell works; you don't need a Developer Prompt. Outputs land in
`build\msvc\bin\<Config>\`.

### VS Code

1. Open the repository folder. CMake Tools detects `CMakePresets.json` and configures it.
2. In the status bar, pick configure preset **msvc** and build preset **debug**.
3. Press **F5** ("Debug PotentialFlowLab"). It builds first, then starts the debugger.
4. Tests: the CMake Tools *Testing* panel, or the "Debug unit tests" launch configuration.

### Shader hot reload (Debug only)

Debug builds read the GLSL files straight from `shaders/`. Edit a shader, then press **Ctrl+R**
(or *File → Reload shaders*). A compile error is printed to the console, and the previous working
shader stays active. Release builds use only the copies embedded in the exe. To make Debug use
the embedded copies too, configure with `-DPFL_SHADERS_FROM_DISK=OFF`.

## Packaging for distribution

```powershell
cmake --build --preset release --target dist
```

This produces `dist\PotentialFlowLab.zip`, containing `PotentialFlowLab.exe` and `README.txt`.

### Verifying the exe is standalone

Open **Developer PowerShell for VS 2022**, which puts `dumpbin` on `PATH`:

```powershell
dumpbin /dependents build\msvc\bin\Release\PotentialFlowLab.exe
dumpbin /headers    build\msvc\bin\Release\PotentialFlowLab.exe | findstr subsystem
```

- The dependency list should contain **only Windows system DLLs**. Currently those are:
  `COMDLG32.dll ole32.dll KERNEL32.dll USER32.dll GDI32.dll SHELL32.dll IMM32.dll` (COMDLG32
  provides the Open / Save dialogs; ole32 starts the Windows Imaging Component that writes PNGs).
- `OPENGL32.dll` doesn't appear because GLFW loads it at runtime.
- `MFPlat.DLL` and `MFReadWrite.dll` (Media Foundation, for MP4 recording) are listed separately
  as *delay load* dependencies: they are loaded only when an MP4 is recorded. Windows "N" editions
  without the Media Feature Pack lack them, and the program still starts there; only MP4 export
  then reports that it is unavailable (GIF export still works).
- Any `VCRUNTIME140*.dll`, `MSVCP140*.dll`, `api-ms-win-crt-*.dll` or `glfw3.dll` means the
  static linking broke.
- The subsystem should read `2 subsystem (Windows GUI)`, which means no console window.

As a final check, copy the exe alone into an empty folder on another PC and double-click it.

## How the standalone build works

- **Static runtime:** `CMAKE_MSVC_RUNTIME_LIBRARY = MultiThreaded$<$<CONFIG:Debug>:Debug>` (`/MT`,
  `/MTd`). It is set before any dependency is fetched, so they all use it.
- **Static dependencies:** `BUILD_SHARED_LIBS=OFF`. ImGui/ImPlot are compiled into `pfl_imgui`.
- **Embedded assets:** `cmake/EmbedResources.cmake` turns `shaders/*` and `examples/*.json` into
  byte arrays inside the exe at build time, using a generated `.cpp` that is rebuilt when a shader changes. Byte arrays
  (not string literals) avoid MSVC's string-literal size limit and will also work for fonts and
  icons.
- **Windows resources:** `resources/app.rc.in` adds the app icon (`resources/app.ico`, regenerate
  with `python tools/make_icon.py`) and version info taken from `project(VERSION ...)`.
- **No console in Release:** the exe links with `/SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup`, so
  `main()` stays a normal `main()`.
- **No OpenGL 3.3:** a Windows message box explains the requirement (and that VMs / Remote
  Desktop often lack it), then the program exits. To test this message in a Debug build, set
  `$env:PFL_FORCE_GL_FAIL=1` before running.
- **Laptops:** the exe exports `NvOptimusEnablement` / `AmdPowerXpressRequestHighPerformance` so
  hybrid-graphics laptops use the discrete GPU.

## Repository layout

```
src/core/     pure math/physics (no OpenGL/ImGui): flow elements, superposition, wall images,
              scene model and warnings, undo history, RK4 streamlines, tracer particles
src/render/   OpenGL wrappers: shaders (with error logs), buffers, framebuffer, particle layer
src/ui/       ImGui/ImPlot: menus, toolbar, parameter dialog, plot view, status bar
src/io/       JSON scene save/load (nlohmann/json), built-in example scenes, CSV exports
src/app/      main(), window/context setup, main loop, Windows helpers
shaders/      GLSL sources (embedded at build time)
tests/        doctest unit tests (link only pfl_core and pfl_io)
resources/    icon, Windows .rc template, README.txt for the zip
cmake/        dependency pins, warnings, resource embedding, packaging
third_party/  generated glad loader (see third_party/glad/README.md)
tools/        helper scripts
docs/         user guide and report notes
examples/     example scenes (cylinder, lifting cylinder, Rankine bodies, corner flow, ...), also
              built into the exe for the Examples menu (list them in src/io/CMakeLists.txt and
              src/io/examples.cpp)
```

## Dependencies (pinned in `cmake/Dependencies.cmake`)

| Library | Version |
|---|---|
| GLFW | 3.5.1 |
| Dear ImGui (docking branch) | v1.93.0-docking |
| ImPlot | master @ `851cdbd` (2026-10-10). v1.0 plus fixes needed for ImGui ≥ 1.92.9 |
| nlohmann/json | v3.12.0 |
| doctest | v2.5.3 |
| glad | glad2 2.0.8, GL 4.6 core + KHR_debug (generated, committed) |

