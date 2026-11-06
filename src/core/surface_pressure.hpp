// surface_pressure.hpp: pressure coefficient on the bodies of a flow, and the forces on them.
//
// Cp = 1 - |V|^2 / V_inf^2 (Bernoulli with the free stream as reference). The flow is steady,
// incompressible and irrotational, so Bernoulli holds between any two points and Cp is defined
// everywhere, but only when a free stream gives the reference speed (FlowModel::freeStreamSpeed > 0).
//
// Bodies. Any streamline can be replaced by a solid wall; the bodies built by superposition are
// bounded by the dividing streamlines through the stagnation points. findBodySurfaces()
//   1. finds the stagnation points (W = 0) by Newton's method on the analytic W(z), started from a
//      grid of seeds,
//   2. finds the streamlines through each one by sampling v on a small circle around it: they
//      leave or arrive where v is radial (four directions at a simple stagnation point, six at a
//      double one),
//   3. traces those streamlines (RK4 on v/|v|, outgoing ones forward and arriving ones backward)
//      until they reach another stagnation point, a singular point or the edge of the region.
// Two streamlines that leave one stagnation point and meet again at another (or one that returns
// to its start) enclose a closed body: cylinder, lifting cylinder, Rankine oval. With other bodies
// nearby they only nearly meet (see BodySurface::closureGap). Two that leave
// one stagnation point and both run out of the region bound an open body (Rankine half-body), and
// likewise two that arrive at one (the half-body of a sink).
//
// Limitations: only bodies with stagnation points on them are found. A cylinder with
// Gamma > 4 pi a V_inf has its single stagnation point off the surface; what is found then is the
// closed dividing streamline around it, which encloses the fluid circulating with the cylinder
// (a valid body too, with the same forces). Bodies with more than two stagnation points on the
// surface are not assembled.
#pragma once

#include "core/scene.hpp"
#include "core/streamline.hpp"
#include "core/vec2.hpp"

#include <optional>
#include <vector>

namespace pfl {

inline double pressureCoefficient(double speed, double freeStreamSpeed) {
    const double ratio = speed / freeStreamSpeed;
    return 1.0 - ratio * ratio;
}

struct BodySurface {
    bool closed = false;
    // Always ordered with the body on the right, so perp(next - point) points out of it.
    // Closed: clockwise (over the upper side first when the free stream points right) from the
    // front stagnation point back to it; the last point repeats the first.
    // Open: from the far end of one branch through the stagnation point to the far end of the
    // other (for a source's half-body: lower branch first, "upper" = left of the free stream).
    std::vector<Vec2> points;
    std::vector<double> cp;
    std::vector<double> arcLength;  // from the (front) stagnation point; open: negative before it
    std::vector<double> streamwise; // closed: x/c, 0 at the most upstream point, 1 at the most downstream;
                                    // open: distance downstream of the stagnation point
    std::vector<double> angle;      // closed: polar angle about the centroid, from the free-stream
                                    // direction, counter-clockwise, degrees in [0, 360); open: empty
    std::vector<Vec2> stagnationPoints;
    Vec2 centroid;                  // closed: area centroid; open: the stagnation point
    // Closed: how far the two surfaces miss the rear stagnation point, where a short straight segment
    // closes the body. 0 for an exactly closed body (a lone body in a uniform stream, symmetric about
    // the line through its stagnation points). Other bodies nearby make the flow asymmetric, so the
    // dividing streamlines no longer meet; the body is accepted while the gap stays below 10% of the
    // distance between its stagnation points. (A Rankine oval at incidence, for instance, is not
    // closed at all: part of the source's fluid escapes downstream.)
    double closureGap = 0.0;
};

struct SurfaceSettings {
    double stepFraction = 1.0 / 2500.0; // tracing step as a fraction of the body region's size
    int maxPointsPerBranch = 40'000;
    int seedGrid = 32;                  // Newton seeds per side of the search region
};

// The stagnation points of `field` inside `region` (excluding singular points), deduplicated.
// `speedScale` sets the |W| that counts as zero (1e-9 of it).
std::vector<Vec2> findStagnationPoints(const FlowField& field, const Rect& region, double speedScale, int seedGrid = 32);

// The stagnation point Newton's method reaches from `from`, if it is within `reach`: where a known
// stagnation point has moved to after a small change of the flow. Nullopt if it has gone.
std::optional<Vec2> followStagnationPoint(const FlowField& field, Vec2 from, double reach, double speedScale);

// The dividing streamlines through the stagnation point `point`: the lines of step 2 above, traced
// (step 3) until they reach one of the stagnation points in `stagnation`, a singular point or the
// edge of `bounds`, which sets the step and need not contain `point`. Each line runs in the flow
// direction with one arrow halfway along it, so a line arriving at `point` ends there. No seeds.
StreamlineSet dividingStreamlines(const FlowField& field, Vec2 point, std::vector<Vec2> stagnation, const Rect& bounds,
                                  const SurfaceSettings& settings = {});

// The bodies near the singular points of the flow, traced far enough to cross `view`. Empty
// without a free stream.
std::vector<BodySurface> findBodySurfaces(const FlowModel& model, const Rect& view,
                                          const SurfaceSettings& settings = {});

// Force per unit span on a closed body from integrating the surface pressure,
// F / q_inf = -sum Cp n ds with the outward normal n, split along and across the free stream.
struct BodyForces {
    double drag = 0.0;  // D' / q_inf, along the free stream (a length)
    double lift = 0.0;  // L' / q_inf, perpendicular to it, positive 90 degrees counter-clockwise
    double chord = 0.0; // extent along the free stream
    double circulation = 0.0; // sum of the clockwise vortex strengths inside the body
    // Every singular point of the flow is inside the body (no other bodies, walls or sector flows),
    // so the exact values are D' = 0 (d'Alembert) and L' / q_inf = 2 Gamma / V_inf (Kutta-Joukowski).
    bool isolated = false;
};
BodyForces bodyForces(const BodySurface& body, const FlowModel& model);

// A circular cylinder: uniform flow plus one doublet aligned with it (and optionally a vortex at the
// same point), with the surface found at r = sqrt(mu / V_inf).
struct CylinderReference {
    Vec2 centre;
    double radius = 1.0;
    double circulation = 0.0;    // Gamma, > 0 clockwise
    double freeStreamAngle = 0.0;
    double freeStreamSpeed = 1.0;

    // Exact Cp at polar angle theta about the centre (radians):
    // Cp = 1 - (2 sin(theta - alpha) + Gamma / (2 pi a V_inf))^2.
    double cp(double theta) const;
};
std::optional<CylinderReference> cylinderReference(const BodySurface& body, const FlowModel& model);

} // namespace pfl
