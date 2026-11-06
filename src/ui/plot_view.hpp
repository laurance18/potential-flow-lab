// plot_view.hpp: the main (x, y) plot.
//
// An ImPlot plot with labelled axes and equal x/y scales (scroll to zoom, drag to pan). Drawn
// back to front: the optional |V| or Cp bands, the region outside the flow domain (darkened), the
// tracer-particle layer (a texture the app renders), the optional equipotential lines, streamlines
// with direction arrows, the dye streaks, the body surfaces coloured by Cp with pressure arrows, the
// dye rakes, sector and symmetry walls, one marker per element (wall images faded), the stagnation
// points and the colour scale.
//
// Interaction: hovering a marker names the element; double-clicking it asks the caller to edit it,
// and dragging it moves it (every element but a uniform flow; a wall moves parallel to itself).
// Dragged positions are rounded to about a pixel, or to a coarse grid with Shift held, and Esc puts
// the element back. Dye rakes are placed by dragging across the plot after the toolbar's rake tool
// was chosen, moved by dragging their line, turned or stretched by dragging an end, and edited or
// deleted from a right-click menu. While anything is grabbed, dragging no longer pans the plot.
// Clicking a stagnation point shows or hides its dividing streamlines, which are drawn above the dye.
#pragma once

#include "core/flow_field.hpp"
#include "core/scene.hpp"
#include "core/streamline.hpp"
#include "core/surface_pressure.hpp"
#include "core/vec2.hpp"

#include <imgui.h>

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace pfl {

struct ParticleLayer {
    ImTextureID texture = 0;             // 0: nothing to show
    ImDrawCallback beginBlend = nullptr; // sets premultiplied-alpha blending for the texture
};

// Filled contours of |V| or Cp: the quantity sampled on a grid over the view, drawn in `bands`
// equal bands from low to high. Values beyond either end use the end band (|V| is +infinity and Cp
// -infinity at singular points). |V| bands use Viridis, Cp bands the diverging pressureColor().
struct ContourMap {
    Background quantity = Background::None;
    ScalarGrid values;
    double low = 0.0;
    double high = 0.0;
    int bands = 0; // 0 with quantity Cp: no free stream, so Cp is undefined
};

// Body surfaces drawn coloured by Cp (on the scale lowestCp..1) with pressure arrows.
struct SurfaceLayer {
    const std::vector<BodySurface>& bodies; // empty when switched off
    double lowestCp = -3.0;
    std::optional<Vec2> marker; // the surface point hovered in the Cp panel
};

// Dye streaks, drawn in the colour of the rake that released them, in pulses of `pulseLength`
// releases (alternately bright and dim) and fading out over the last part of their life.
struct DyeLayer {
    std::span<const DyeStreak> streaks;
    double maxAge = 1.0;
    int pulseLength = 15;
};

// What the app computes from the flow for the plot to draw.
struct PlotLayers {
    const StreamlineSet& streamlines;
    const StreamlineSet& equipotentials; // empty when they are switched off
    const ContourMap& contours;
    SurfaceLayer surfaces;
    ParticleLayer particles;
    const std::vector<Vec2>& stagnationPoints; // empty when they are switched off
    const std::vector<Vec2>& dividingPoints;   // stagnation points showing their dividing streamlines
    const StreamlineSet& dividingLines;        // those streamlines
    DyeLayer dye;
};

struct ElementMove {
    std::size_t index;
    Vec2 position; // new anchor point (see anchorPoint())
};

struct RakeChange {
    std::optional<std::size_t> index; // rake to replace; nullopt adds a new one
    DyeRake rake;
};

struct PlotViewResult {
    Rect limits;                   // world rectangle currently shown
    ImVec2 pixelSize{0.0f, 0.0f};  // size of the plot area in pixels
    std::optional<Vec2> mouse;     // world position under the mouse, if it is over the plot
    std::optional<std::size_t> editRequest;
    std::optional<ElementMove> move;      // an element dragged to a new position this frame
    std::optional<RakeChange> rakeChange; // a rake placed, dragged or edited this frame
    std::optional<std::size_t> rakeRemove;
    std::optional<Vec2> dividingAdd;           // a stagnation point clicked to show its dividing streamlines
    std::optional<std::size_t> dividingRemove; // index in PlotLayers::dividingPoints, clicked to hide them
    bool dragging = false;                // an element or a rake is being dragged
};

class PlotView {
public:
    // Shows exactly `view` from the next frame on (after loading a file or resetting the view).
    void setView(const Rect& view) { m_pendingView = view; }

    // Starts (or, if already started, cancels) placing a dye rake: the next drag on the plot draws it.
    void toggleRakePlacement();
    bool placingRake() const { return m_placingRake; }

    PlotViewResult draw(ImVec2 pos, ImVec2 size, const Scene& scene, const FlowModel& model, const PlotLayers& layers);

private:
    enum class Target { Element, Rake, RakeStart, RakeEnd };
    struct Grab {
        Target target = Target::Element;
        std::size_t index = 0;
        Vec2 pressWorld;     // mouse position when the button went down
        Vec2 anchor;         // element: its anchor point then
        DyeRake rake;        // rake: the rake then
        bool moved = false;  // the mouse has moved past the drag threshold
    };

    void updateGrab(const Scene& scene, Vec2 mouseWorld, const Rect& limits, float heightPx, PlotViewResult& result);
    void updateRakePlacement(Vec2 mouseWorld, bool hovered, const Rect& limits, float heightPx, PlotViewResult& result);

    std::optional<Rect> m_pendingView;
    std::optional<Grab> m_grab;
    bool m_grabbableHovered = false; // last frame: something draggable under the mouse (blocks panning)
    bool m_placingRake = false;
    std::optional<Vec2> m_rakeStart;       // while a new rake is being drawn
    std::optional<std::size_t> m_rakeMenu; // rake whose right-click menu is open
};

} // namespace pfl
