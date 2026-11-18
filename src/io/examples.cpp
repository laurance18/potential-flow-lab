#include "io/examples.hpp"

#include "embedded_resources.hpp"

#include <string>

namespace pfl {
namespace {

// Keep in step with the files listed in src/io/CMakeLists.txt.
constexpr ExampleScene k_examples[] = {
    {"uniform_45deg.json", "Uniform flow at 45°", "Uniform flow U∞ = 1 at α = 45°: straight, parallel streamlines."},
    {"rankine_half_body.json", "Rankine half-body",
     "Uniform flow U∞ = 1 and a source Λ = 2π at (-1, 0).\n"
     "The dividing streamline is a half-body, with its stagnation point at (-2, 0)\n"
     "and a half-width of Λ/2U∞ = π far downstream."},
    {"rankine_oval.json", "Rankine oval",
     "Uniform flow U∞ = 1, a source Λ = 2π at (-1, 0) and an equal sink at (1, 0).\n"
     "The dividing streamline closes into an oval."},
    {"cylinder.json", "Cylinder",
     "Uniform flow U∞ = 1 and a doublet μ = 1: a cylinder of radius a = √(μ/U∞) = 1.\n"
     "Stagnation points at 0° and 180°, |V| = 2U∞ at the top and bottom, no lift or drag."},
    {"lifting_cylinder.json", "Lifting cylinder",
     "The cylinder plus a clockwise vortex Γ = 2π at its centre.\n"
     "The stagnation points move down to 210° and 330°; the lift is L′ = ρ∞U∞Γ (Kutta–Joukowski)."},
    {"source_near_wall.json", "Source near a wall",
     "A source Λ = 2π at (0, 1) above the symmetry wall y = 0.\n"
     "An image source below the wall keeps the flow from crossing it."},
    {"channel_source.json", "Source in a channel",
     "A source Λ = 2π at (0, 0.4) between symmetry walls at y = ±1.\n"
     "Two infinite rows of images keep the flow in the channel; half of it leaves each way,\n"
     "at speed Λ/2h = π/2 far from the source."},
    {"source_sink_vortex.json", "Source, sink and vortex",
     "A source Λ = 2π at (-2, 0), a sink of equal strength at (0, 0) and a clockwise vortex Γ = 2π at (2, 0)."},
    {"corner_flow_90deg.json", "Flow in a 90° corner", "Sector flow between walls at 0° and 90° (n = 2)."},
    {"flat_plate_edge.json", "Flow around a flat plate edge",
     "Sector flow with walls at 0° and 360° (n = 1/2):\n"
     "flow around the edge of a semi-infinite flat plate along +x."},
    {"panel_cylinder.json", "Source panels: cylinder (8 panels)",
     "Anderson's example: a circle of radius 1 cut into 8 straight source panels, in U∞ = 1.\n"
     "Cp at the 8 control points is exactly 1 − 4 sin²θ; between them the octagon shows.\n"
     "Raise N in the Parameters panel to round it off."},
    {"panel_ellipse.json", "Source panels: ellipse at incidence",
     "A 4:1 ellipse with 60 source panels, at 15° to the stream U∞ = 1,\n"
     "against the exact Cp from the conformal map. No lift: source panels carry no circulation."},
    {"panel_lifting_cylinder.json", "Source panels: lifting cylinder",
     "A vortex Γ = 2π inside a source-panel circle (48 panels): the panels keep the flow off the\n"
     "surface and the vortex gives the circulation. Compare with the Lifting cylinder example."},
    {"panel_bluff_bodies.json", "Source panels: square and circle",
     "A square and a circle, solved together so each feels the other, in U∞ = 1.\n"
     "The square's corners are where the speed peaks."},
    {"panel_naca2412.json", "Source panels: NACA 2412 airfoil",
     "A NACA 2412 of chord 3 at 5° angle of attack, 100 source panels. Without circulation it gets no\n"
     "lift, and the flow has to turn round the sharp trailing edge, where the panels need huge strengths.\n"
     "That is what the Kutta condition of a vortex sheet fixes."},
    {"panel_drawn_shape.json", "Source panels: a drawn shape",
     "A body drawn as 14 points, joined by a smooth curve and cut into 80 panels.\n"
     "Draw your own with the pencil tool below the source-panel body."},
};

} // namespace

std::span<const ExampleScene> exampleScenes() { return k_examples; }

SceneLoadResult loadExampleScene(const ExampleScene& example) {
    const auto text = examples::findText(std::string("examples/") + example.file);
    if (!text) return {.scene = std::nullopt, .error = std::string("Example ") + example.file + " is not built in."};
    return sceneFromJson(*text);
}

} // namespace pfl
