// sweep.hpp: a parameter sweep. One element parameter is stepped through a range, and at each
// value the flow is rebuilt and analysed like the plot does it: the stagnation points in the view,
// the bodies with their surface-pressure forces (core/surface_pressure.hpp), and the thin airfoils
// with their lift and moment (core/vortex_sheet.hpp). Sweeping an airfoil's α gives its c_l-α curve.
#pragma once

#include "core/parameters.hpp"
#include "core/scene.hpp"
#include "core/vec2.hpp"

#include <optional>
#include <vector>

namespace pfl {

// `steps` evenly spaced values from `from` to `to`, both included (just `from` if steps < 2).
std::vector<double> sweepValues(double from, double to, int steps);

struct SweepBody {
    bool closed = false;
    Vec2 centroid;            // closed: area centroid; open: the stagnation point
    double chord = 0.0;       // closed: extent along the free stream
    double circulation = 0.0; // closed: clockwise vortex strength inside
    double lift = 0.0;        // closed: L'/q_inf and D'/q_inf from the surface Cp (NaN for a half-body)
    double drag = 0.0;
    double exactLift = 0.0;   // 2 Gamma / V_inf (Kutta-Joukowski) for an isolated closed body, else NaN
    double cpMin = 0.0;
};

// A vortex-sheet airfoil (vortex_sheet.hpp) at one step: what the panels give and what thin-airfoil
// theory says for an airfoil alone at the same angle of attack.
struct SweepAirfoil {
    Vec2 position;            // the middle of the chord
    double chord = 0.0;
    double alpha = 0.0;       // angle of attack (rad)
    double circulation = 0.0; // Gamma
    double lift = 0.0;        // L'/q_inf from Delta Cp along the chord
    double exactLift = 0.0;   // 2 Gamma / V_inf (Kutta-Joukowski)
    double cpMin = 0.0;       // lowest Cp on either side at the control points
    double cmQuarter = 0.0;
    double clTheory = 0.0;    // 2 pi (A0 + A1/2)
    double cmTheory = 0.0;    // (pi/4)(A2 - A1)
};

struct SweepSample {
    double value = 0.0;
    std::vector<Vec2> stagnationPoints;
    std::vector<SweepBody> bodies;       // empty without a free stream
    std::vector<SweepAirfoil> airfoils;  // likewise
};

// The scene with `ref` set to `value`, analysed over `view`. nullopt if the value makes the element
// invalid (a sector with θ₂ ≤ θ₁, a negative uniform speed, ...).
std::optional<SweepSample> sampleSweep(const Scene& scene, ParameterRef ref, double value, const Rect& view);

} // namespace pfl
