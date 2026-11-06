// flow_element.hpp: elementary potential-flow solutions.
//
// Each element gives its complex potential F = phi + i psi and complex velocity W = dF/dz = u - i v
// at a point, using the implementation form of the course reference equations (section 4), with
// zeta = z - z0. Sign conventions follow Anderson, Ch. 3:
//     Lambda > 0 is a source, Gamma > 0 turns clockwise, mu = kappa / 2 pi.
// All angles are radians measured counter-clockwise from +x. Pure math, no graphics.
#pragma once

#include "core/vec2.hpp"

#include <complex>
#include <numbers>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace pfl {

using Complex = std::complex<double>;

constexpr Complex toComplex(Vec2 p) { return {p.x, p.y}; }

// u = Re W, v = -Im W.
constexpr Vec2 velocityFromW(Complex w) { return {w.real(), -w.imag()}; }

// W = V e^{-i alpha}.
struct UniformFlow {
    double speed = 1.0; // V_inf, >= 0
    double angle = 0.0; // alpha

    friend bool operator==(const UniformFlow&, const UniformFlow&) = default;
};

// W = Lambda / (2 pi zeta).
struct SourceSink {
    double strength = 1.0; // Lambda: volume flow rate per unit depth; > 0 source, < 0 sink
    Vec2 position;

    friend bool operator==(const SourceSink&, const SourceSink&) = default;
};

// W = i Gamma / (2 pi zeta).
struct Vortex {
    double circulation = 1.0; // Gamma; > 0 clockwise
    Vec2 position;

    friend bool operator==(const Vortex&, const Vortex&) = default;
};

// W = -mu e^{i beta} / zeta^2. With beta = 0 the flow leaves the centre towards -x and returns
// from +x; uniform flow at angle alpha plus a doublet with beta = alpha is a cylinder of radius
// sqrt(mu / V_inf).
struct Doublet {
    double strength = 1.0; // mu
    Vec2 position;
    double angle = 0.0; // beta: rotates the doublet axis counter-clockwise

    friend bool operator==(const Doublet&, const Doublet&) = default;
};

// F = V (zeta e^{-i theta1})^n with n = pi / (theta2 - theta1): flow inside the wedge between the
// walls at theta1 and theta2 (both streamlines). Defined only inside the wedge.
struct SectorFlow {
    double strength = 1.0; // V, in units of speed / length^(n-1)
    Vec2 vertex;
    double wall1 = 0.0;                     // theta1
    double wall2 = std::numbers::pi / 2.0; // theta2; 0 < theta2 - theta1 <= 2 pi

    double opening() const { return wall2 - wall1; }
    double exponent() const { return std::numbers::pi / opening(); }

    friend bool operator==(const SectorFlow&, const SectorFlow&) = default;
};

using FlowElement = std::variant<UniformFlow, SourceSink, Vortex, Doublet, SectorFlow>;

// Complex velocity and potential of one element at p. Infinite or NaN at a singular point; the
// potentials of sources and vortices use the principal branch of log (cut along -x from the centre).
Complex complexVelocity(const FlowElement& element, Vec2 p);
Complex complexPotential(const FlowElement& element, Vec2 p);
inline Vec2 velocity(const FlowElement& element, Vec2 p) { return velocityFromW(complexVelocity(element, p)); }

// The point where the element's velocity is infinite: the centre of a source, vortex or doublet,
// or the vertex of a sector flow opening wider than pi (n < 1). Uniform flow has none.
std::optional<Vec2> singularPoint(const FlowElement& element);

// `element` moved by `offset`. A uniform flow is unchanged.
FlowElement translated(FlowElement element, Vec2 offset);

// Angle of p around the sector vertex measured from the first wall, in [0, 2 pi).
double sectorLocalAngle(const SectorFlow& sector, Vec2 p);

// False only for points outside a sector flow's wedge, where that element is undefined.
bool inElementDomain(const FlowElement& element, Vec2 p);

// A message describing why the parameters are not physically valid, or nullopt if they are.
std::optional<std::string> validationError(const FlowElement& element);

// An infinite row of copies of a source, vortex or doublet, at z0 + n T for every integer n. The sum
// over the row has a closed form, with k = i pi / T and w = k (z - z0):
//     sum 1 / (z - z0 - n T)     = k coth w         (symmetric partial sums)
//     sum 1 / (z - z0 - n T)^2   = k^2 csch^2 w
// so a row is the element's own W and F with 1 / zeta -> k coth w, 1 / zeta^2 -> k^2 csch^2 w and
// log zeta -> log sinh w (up to a constant).
struct PeriodicRow {
    FlowElement element;   // the copy n = 0; a source, vortex or doublet, never uniform or sector
    Vec2 period;           // T, not zero
    bool original = false; // the copy n = 0 is an element of the scene rather than an image

    friend bool operator==(const PeriodicRow&, const PeriodicRow&) = default;
};

// Complex velocity and potential of the whole row at p. Infinite or NaN at a copy. Between the two
// walls of a channel F has one branch cut, running along the channel from the element.
Complex complexVelocity(const PeriodicRow& row, Vec2 p);
Complex complexPotential(const PeriodicRow& row, Vec2 p);

// The singular point of the copy closest to p along the row.
Vec2 nearestCopy(const PeriodicRow& row, Vec2 p);

// The copies whose singular points lie in `region`, at most a few hundred (a stretch from the middle
// of the region if there are more than that, which only happens far beyond what can be resolved).
std::vector<FlowElement> rowCopies(const PeriodicRow& row, const Rect& region);

} // namespace pfl
