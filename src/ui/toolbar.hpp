// toolbar.hpp: left-hand icon toolbar with one button per element type (uniform flow, source/sink,
// vortex, doublet, sector flow, symmetry wall) and, below them, the dye rake tool. Hovering a button
// shows its name; clicking an element asks the caller to open the parameter dialog for a new element
// of that type, clicking the rake tool asks it to start placing a rake on the plot.
#pragma once

#include "core/scene.hpp"

#include <imgui.h>

#include <optional>

namespace pfl {

struct ToolbarAction {
    std::optional<ElementKind> element;
    bool placeRake = false;
    bool drawBody = false; // start (or cancel) drawing the outline of a source-panel body
};

float toolbarWidth();

// `placingRake` highlights the rake tool while a rake is being placed, `drawingBody` the outline tool
// (just below the source-panel body) while an outline is being drawn.
ToolbarAction drawToolbar(ImVec2 pos, ImVec2 size, bool placingRake, bool drawingBody);

} // namespace pfl
