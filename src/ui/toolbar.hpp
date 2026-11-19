// toolbar.hpp: left-hand icon toolbar. At the top, the click-to-place switch; below it one button per
// element type (uniform flow, source/sink,
// vortex, doublet, sector flow, symmetry wall; below a separator the panel methods: the source-panel
// body, its outline pencil and the vortex-sheet airfoil) and, below them, the dye rake tool. Hovering a button
// shows its name; clicking an element asks the caller to open the parameter dialog for a new element
// of that type (or, with click-to-place on, to let the user click its position on the plot first),
// clicking the rake tool asks it to start placing a rake on the plot.
#pragma once

#include "core/scene.hpp"

#include <imgui.h>

#include <optional>

namespace pfl {

struct ToolbarAction {
    std::optional<ElementKind> element;
    bool placeRake = false;
    bool drawBody = false; // start (or cancel) drawing the outline of a source-panel body
    bool toggleClickToPlace = false;
};

// What the toolbar highlights.
struct ToolbarState {
    bool clickToPlace = false;          // new elements are placed with a click on the plot
    std::optional<ElementKind> placing; // the element being placed that way
    bool placingRake = false;
    bool drawingBody = false; // the outline tool, just below the source-panel body
};

float toolbarWidth();

ToolbarAction drawToolbar(ImVec2 pos, ImVec2 size, const ToolbarState& state);

} // namespace pfl
