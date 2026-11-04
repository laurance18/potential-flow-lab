// main_menu.hpp: the application's main menu bar (File, Edit, Examples, Elements, View, Verify) and its shortcuts.
//
// File: New, Open, Save, Save As, the image and CSV exports, Record animation (or Stop recording),
// Clear All, Quit (and Reload shaders in Debug builds). Edit: Undo (Ctrl+Z) and Redo (Ctrl+Y or
// Ctrl+Shift+Z) of element and dye-rake edits (core/history.hpp). Examples: the built-in example scenes
// (io/examples.hpp), each loaded with one click. Elements: see elements_menu.hpp. View: the
// Parameters panel, reset the view and the plot settings that are saved with the scene (streamline
// density, arrows, wall images, equipotentials, |V| or Cp background, surface pressure, stagnation
// points, particles), edited in place. Verify: the physics checks (verification_window.hpp), F5.
// At the right end of the bar: the recording in progress with
// a Stop button, or a short notice (e.g. the file just saved).
#pragma once

#include "core/scene.hpp"
#include "ui/elements_menu.hpp"

#include <cstddef>
#include <optional>
#include <string>

namespace pfl {

struct RecordingProgress {
    double seconds = 0.0; // of video recorded so far
    double total = 0.0;
};

// What the menus show besides the scene.
struct MainMenuState {
    bool canReloadShaders = false; // Debug builds reading shaders from disk (Ctrl+R)
    bool parameterPanel = false;   // View → Parameters is ticked
    bool canUndo = false;
    bool canRedo = false;
    std::optional<RecordingProgress> recording;
    std::string notice;
};

struct MainMenuActions {
    bool newScene = false;
    bool open = false;
    bool save = false;
    bool saveAs = false;
    bool undo = false;
    bool redo = false;
    bool exportImage = false;
    bool exportField = false;
    bool exportLines = false;
    bool exportSurface = false;
    bool record = false;
    bool stopRecording = false;
    bool toggleParameterPanel = false;
    bool verify = false; // open the physics checks
    bool clearAll = false;
    bool quit = false;
    bool reloadShaders = false;
    bool resetView = false;
    bool plotChanged = false; // a View setting in `plot` was edited
    std::optional<std::size_t> example; // index into exampleScenes()
    ElementsMenuAction elements;
};

// Draws the menu bar and handles its keyboard shortcuts.
MainMenuActions drawMainMenu(const Scene& scene, PlotSettings& plot, const MainMenuState& state);

// "Discard unsaved changes?" confirmation. Open it, then draw it every frame; it returns true
// (discard) or false (cancel) once the user has chosen.
void openDiscardChangesDialog();
std::optional<bool> drawDiscardChangesDialog();

} // namespace pfl
