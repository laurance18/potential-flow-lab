#include "core/flow_element.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace pfl {
namespace {

using std::numbers::pi;
constexpr double k_twoPi = 2.0 * pi;
constexpr Complex k_i{0.0, 1.0};

// Points within this angle of a sector wall count as on the wall, so rounding in atan2 does not
// push wall points out of the wedge.
constexpr double k_angleTolerance = 1e-9;

template <typename... Ts>
struct Overloaded : Ts... {
    using Ts::operator()...;
};

bool isFinite(Vec2 p) { return std::isfinite(p.x) && std::isfinite(p.y); }

// Beyond this |Re w|, coth w = +-1 and csch^2 w = 0 to machine precision, while cosh and sinh would
// overflow a little further on (far down a narrow channel).
constexpr double k_farRow = 20.0;
constexpr int k_maxRowCopies = 256;

Complex coth(Complex w) {
    if (std::abs(w.real()) > k_farRow) return w.real() > 0.0 ? 1.0 : -1.0;
    return std::cosh(w) / std::sinh(w);
}

Complex cschSquared(Complex w) {
    if (std::abs(w.real()) > k_farRow) return 0.0;
    const Complex s = std::sinh(w);
    return 1.0 / (s * s);
}

// Principal branch of log(sinh w).
Complex logSinh(Complex w) {
    if (std::abs(w.real()) <= k_farRow) return std::log(std::sinh(w));
    // sinh(a + i b) -> e^|a| / 2 (+-cos b + i sin b), the sign being that of a.
    const double cosine = std::cos(w.imag());
    return {std::abs(w.real()) - std::numbers::ln2, std::atan2(std::sin(w.imag()), w.real() > 0.0 ? cosine : -cosine)};
}

// k = i pi / T.
Complex rowWavenumber(const PeriodicRow& row) { return k_i * pi / toComplex(row.period); }

} // namespace

FlowElement translated(FlowElement element, Vec2 offset) {
    std::visit(Overloaded{
                   [](UniformFlow&) {},
                   [offset](SectorFlow& e) { e.vertex += offset; },
                   [offset](auto& e) { e.position += offset; },
               },
               element);
    return element;
}

double sectorLocalAngle(const SectorFlow& sector, Vec2 p) {
    const Vec2 d = p - sector.vertex;
    double angle = std::fmod(std::atan2(d.y, d.x) - sector.wall1, k_twoPi);
    if (angle < 0.0) angle += k_twoPi;
    // Just below the first wall wraps to ~2 pi; snap it back unless the wedge really is 2 pi wide.
    if (angle > k_twoPi - k_angleTolerance && sector.opening() < k_twoPi - k_angleTolerance) angle = 0.0;
    return angle;
}

Complex complexVelocity(const FlowElement& element, Vec2 p) {
    return std::visit(
        Overloaded{
            [](const UniformFlow& e) -> Complex { return e.speed * std::polar(1.0, -e.angle); },
            [p](const SourceSink& e) -> Complex { return e.strength / (k_twoPi * toComplex(p - e.position)); },
            [p](const Vortex& e) -> Complex { return k_i * e.circulation / (k_twoPi * toComplex(p - e.position)); },
            [p](const Doublet& e) -> Complex {
                const Complex zeta = toComplex(p - e.position);
                return -e.strength * std::polar(1.0, e.angle) / (zeta * zeta);
            },
            [p](const SectorFlow& e) -> Complex {
                // V n e^{-i theta1} |zeta|^{n-1} e^{i (n-1) theta_l}
                const double n = e.exponent();
                const double r = length(p - e.vertex);
                return e.strength * n * std::polar(1.0, -e.wall1) * std::pow(r, n - 1.0) *
                       std::polar(1.0, (n - 1.0) * sectorLocalAngle(e, p));
            },
        },
        element);
}

Complex complexPotential(const FlowElement& element, Vec2 p) {
    return std::visit(
        Overloaded{
            [p](const UniformFlow& e) -> Complex { return e.speed * std::polar(1.0, -e.angle) * toComplex(p); },
            [p](const SourceSink& e) -> Complex { return e.strength / k_twoPi * std::log(toComplex(p - e.position)); },
            [p](const Vortex& e) -> Complex {
                return k_i * e.circulation / k_twoPi * std::log(toComplex(p - e.position));
            },
            [p](const Doublet& e) -> Complex {
                return e.strength * std::polar(1.0, e.angle) / toComplex(p - e.position);
            },
            [p](const SectorFlow& e) -> Complex {
                // V |zeta|^n e^{i n theta_l}
                const double n = e.exponent();
                return e.strength * std::pow(length(p - e.vertex), n) * std::polar(1.0, n * sectorLocalAngle(e, p));
            },
        },
        element);
}

std::optional<Vec2> singularPoint(const FlowElement& element) {
    return std::visit(Overloaded{
                          [](const UniformFlow&) -> std::optional<Vec2> { return std::nullopt; },
                          [](const SourceSink& e) -> std::optional<Vec2> { return e.position; },
                          [](const Vortex& e) -> std::optional<Vec2> { return e.position; },
                          [](const Doublet& e) -> std::optional<Vec2> { return e.position; },
                          [](const SectorFlow& e) -> std::optional<Vec2> {
                              if (e.exponent() < 1.0) return e.vertex;
                              return std::nullopt;
                          },
                      },
                      element);
}

bool inElementDomain(const FlowElement& element, Vec2 p) {
    const auto* sector = std::get_if<SectorFlow>(&element);
    return sector == nullptr || sectorLocalAngle(*sector, p) <= sector->opening() + k_angleTolerance;
}

Complex complexVelocity(const PeriodicRow& row, Vec2 p) {
    const Complex k = rowWavenumber(row);
    return std::visit(
        Overloaded{
            [&](const SourceSink& e) -> Complex {
                return e.strength / k_twoPi * k * coth(k * toComplex(p - e.position));
            },
            [&](const Vortex& e) -> Complex {
                return k_i * e.circulation / k_twoPi * k * coth(k * toComplex(p - e.position));
            },
            [&](const Doublet& e) -> Complex {
                return -e.strength * std::polar(1.0, e.angle) * k * k * cschSquared(k * toComplex(p - e.position));
            },
            [](const auto&) -> Complex { return {}; },
        },
        row.element);
}

Complex complexPotential(const PeriodicRow& row, Vec2 p) {
    const Complex k = rowWavenumber(row);
    return std::visit(
        Overloaded{
            [&](const SourceSink& e) -> Complex {
                return e.strength / k_twoPi * logSinh(k * toComplex(p - e.position));
            },
            [&](const Vortex& e) -> Complex {
                return k_i * e.circulation / k_twoPi * logSinh(k * toComplex(p - e.position));
            },
            [&](const Doublet& e) -> Complex {
                return e.strength * std::polar(1.0, e.angle) * k * coth(k * toComplex(p - e.position));
            },
            [](const auto&) -> Complex { return {}; },
        },
        row.element);
}

Vec2 nearestCopy(const PeriodicRow& row, Vec2 p) {
    const Vec2 centre = *singularPoint(row.element);
    return centre + row.period * std::round(dot(p - centre, row.period) / lengthSquared(row.period));
}

std::vector<FlowElement> rowCopies(const PeriodicRow& row, const Rect& region) {
    const Vec2 centre = *singularPoint(row.element);
    double lo = std::numeric_limits<double>::infinity();
    double hi = -lo;
    const std::array<Vec2, 4> corners = {region.min, region.max, Vec2{region.min.x, region.max.y},
                                         Vec2{region.max.x, region.min.y}};
    for (const Vec2 corner : corners) {
        const double n = dot(corner - centre, row.period) / lengthSquared(row.period);
        lo = std::min(lo, n);
        hi = std::max(hi, n);
    }
    double first = std::ceil(lo);
    double last = std::floor(hi);
    if (last - first >= k_maxRowCopies) {
        first = std::floor(0.5 * (lo + hi)) - k_maxRowCopies / 2;
        last = first + k_maxRowCopies - 1.0;
    }
    std::vector<FlowElement> copies;
    for (double n = first; n <= last; ++n) {
        FlowElement copy = translated(row.element, row.period * n);
        if (region.contains(*singularPoint(copy))) copies.push_back(std::move(copy));
    }
    return copies;
}

std::optional<std::string> validationError(const FlowElement& element) {
    using Result = std::optional<std::string>;
    const Result notFinite = "All parameters must be finite numbers.";
    return std::visit(
        Overloaded{
            [&](const UniformFlow& e) -> Result {
                if (!std::isfinite(e.speed) || !std::isfinite(e.angle)) return notFinite;
                if (e.speed < 0.0) return "Speed must not be negative; reverse the flow with the angle instead.";
                return std::nullopt;
            },
            [&](const SourceSink& e) -> Result {
                if (!std::isfinite(e.strength) || !isFinite(e.position)) return notFinite;
                return std::nullopt;
            },
            [&](const Vortex& e) -> Result {
                if (!std::isfinite(e.circulation) || !isFinite(e.position)) return notFinite;
                return std::nullopt;
            },
            [&](const Doublet& e) -> Result {
                if (!std::isfinite(e.strength) || !isFinite(e.position) || !std::isfinite(e.angle)) return notFinite;
                return std::nullopt;
            },
            [&](const SectorFlow& e) -> Result {
                if (!std::isfinite(e.strength) || !isFinite(e.vertex) || !std::isfinite(e.wall1) ||
                    !std::isfinite(e.wall2)) {
                    return notFinite;
                }
                if (e.opening() <= 0.0 || e.opening() > k_twoPi + k_angleTolerance) {
                    return "Wall angles must satisfy 0 < θ₂ − θ₁ ≤ 360°.";
                }
                return std::nullopt;
            },
        },
        element);
}

} // namespace pfl
