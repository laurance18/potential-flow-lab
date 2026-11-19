// surface_panel.hpp: the surface-pressure panel beside the plot.
//
// Plots Cp along every body surface found in the flow (core/surface_pressure.hpp) against the
// angle around the body, x/c along the free stream or the arc length from the stagnation point, on
// an inverted Cp axis (suction up, as usual in aerodynamics). For a circular cylinder the exact
// curve is drawn under the computed one. Below, for each closed body, the drag and lift per unit
// span from integrating Cp around the surface, next to the exact values for an isolated body
// (d'Alembert: no drag; Kutta-Joukowski: L' = rho V_inf Gamma). Hovering the plot picks the
// nearest surface point, which the main plot marks.
//
// Vortex-sheet airfoils come first: for each, Delta Cp = 2 gamma / V_inf along the chord against
// the thin-airfoil curve, the theory's Fourier coefficients A0 .. A3, and the circulation, c_l (from
// Kutta-Joukowski and from Delta Cp), c_m about the quarter chord, alpha_L0 and the lift slope from
// the panels beside the theory's values.
#pragma once

#include "core/scene.hpp"
#include "core/surface_pressure.hpp"
#include "core/vec2.hpp"

#include <imgui.h>

#include <optional>
#include <vector>

namespace pfl {

// Width of the panel when the window's plot area is `available` pixels wide.
float surfacePanelWidth(float available);

class SurfacePanel {
public:
    // Returns the surface point picked in the Cp plot, if any.
    std::optional<Vec2> draw(ImVec2 pos, ImVec2 size, const std::vector<BodySurface>& bodies,
                             const std::vector<AirfoilLoading>& airfoils, const FlowModel& model);

private:
    enum class Axis { Angle, Chord, ArcLength };
    Axis m_axis = Axis::Angle;
};

} // namespace pfl
