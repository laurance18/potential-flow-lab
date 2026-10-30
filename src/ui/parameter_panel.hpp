// parameter_panel.hpp: live parameter sliders and the parameter sweep (View → Parameters).
//
// Sliders: one per numeric parameter of every element (core/parameters.hpp), angles in degrees.
// Moving one changes the flow live, like dragging an element on the plot: streamlines, contours,
// stagnation points, bodies and the Cp panel follow while the particles keep flowing. Ctrl+click a
// slider to type an exact value. A slider's range is set when it first shows its value: ±2|value|
// rounded up to 1, 2 or 5 × 10ⁿ (at least ±2) for strengths, 0 to that for speeds, ±180° for
// angles (±360° beyond that) and the view for positions. A value entered outside widens it.
//
// Sweep: one parameter stepped from From to To. Play moves it back and forth (or From to To and
// over again) in Duration seconds of animation time; a recording can play it once over its length.
// Export writes the sweep CSV (io/data_export.hpp): the stagnation points and the bodies with
// their forces at each of Steps values. Right-click any slider to sweep its parameter.
#pragma once

#include "core/parameters.hpp"
#include "core/scene.hpp"
#include "core/vec2.hpp"
#include "ui/parameter_dialog.hpp"

#include <imgui.h>

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace pfl {

struct SweepSetup {
    ParameterRef parameter;
    double from = 0.0;
    double to = 1.0;
    int steps = 41;           // values in the CSV export
    double seconds = 5.0;     // one pass From to To when playing
    bool backAndForth = true; // else From to To and jump back
};

struct ParameterPanelResult {
    std::vector<ElementEdit> edits; // slider changes, at most one per element
    bool exportSweep = false;
    bool cancelExport = false;
};

// Width of the panel when the window's plot area is `available` pixels wide.
float parameterPanelWidth(float available);

class ParameterPanel {
public:
    // Advances a playing sweep by `dt` seconds of animation time; returns the element with the swept
    // parameter at its new value. Call every frame, shown or not.
    std::optional<ElementEdit> advance(const Scene& scene, double dt);

    // `exportProgress` is the fraction of the sweep export done while one runs.
    ParameterPanelResult draw(ImVec2 pos, ImVec2 size, const Scene& scene, const Rect& view,
                              std::optional<float> exportProgress);

    // The sweep set up in the panel, if its element is still there.
    std::optional<SweepSetup> sweep(const Scene& scene) const;
    // "#3 Point vortex: Γ from 0 to 12.6" (angles in degrees), empty without a sweep.
    std::string sweepDescription(const Scene& scene) const;

    // Plays the sweep once from From to To over `seconds`, then stops (for a recording).
    void playOnce(double seconds);
    void stop() { m_playing = false; }
    bool playing() const { return m_playing; }

    // Forgets the ranges and the sweep (a new scene was loaded).
    void reset();

private:
    struct Range {
        double low = 0.0;
        double high = 1.0;
    };

    // Forgets the ranges and the sweep when the scene's elements are not the ones they were set for.
    void sync(const Scene& scene);
    Range& range(ParameterRef ref, const ElementParameter& parameter, double value, const Rect& view);
    void chooseSweep(ParameterRef ref, double value, const Range& range);
    bool drawSlider(ParameterRef ref, const ElementParameter& parameter, SceneElement& element, const Rect& view);
    void drawSweep(const Scene& scene, const Rect& view, std::optional<float> exportProgress,
                   ParameterPanelResult& result);

    std::vector<ElementKind> m_kinds; // the scene's element kinds the ranges and the sweep belong to
    std::map<std::pair<std::size_t, std::size_t>, Range> m_ranges;
    std::optional<SweepSetup> m_sweep;
    bool m_playing = false;
    bool m_once = false;      // stop at To (a recording)
    double m_onceSeconds = 0.0;
    double m_time = 0.0;      // animation time since Play
};

} // namespace pfl
