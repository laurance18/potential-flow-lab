#include "core/verification.hpp"

#include "core/dye.hpp"
#include "core/flow_element.hpp"
#include "core/flow_field.hpp"
#include "core/image_system.hpp"
#include "core/particles.hpp"
#include "core/streamline.hpp"
#include "core/surface_pressure.hpp"
#include "core/sweep.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <format>
#include <functional>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <type_traits>
#include <utility>

namespace pfl {
namespace {

constexpr double pi = std::numbers::pi;
constexpr double deg = pi / 180.0;
constexpr double inf = std::numeric_limits<double>::infinity();

// Field comparisons use the default view sampled like the field CSV export (201 x 117 nodes).
const Rect k_view = PlotSettings{}.view;
constexpr int k_gridColumns = 201;
constexpr int k_gridRows = 117;

// Closed form against closed form: only rounding separates them.
constexpr double k_exactTolerance = 1e-12;

// --- Formatting ---------------------------------------------------------------------------

std::string num(double x) {
    if (x == 0.0) return "0";
    if (!std::isfinite(x)) return "—";
    if (std::abs(x) < 1e-4) return std::format("{:.1e}", x);
    return std::format("{:.6g}", x);
}

// Components below 1e-12 of `scale` are rounding and shown as 0.
std::string vec(Vec2 v, double scale = 1.0) {
    const auto snap = [scale](double x) { return std::abs(x) < 1e-12 * scale ? 0.0 : x; };
    return std::format("({}, {})", num(snap(v.x)), num(snap(v.y)));
}

std::string point(Vec2 p) { return std::format("({:.4g}, {:.4g})", p.x + 0.0, p.y + 0.0); }

// 23517 -> "23,517".
std::string count(int n) {
    std::string digits = std::to_string(n);
    for (int i = static_cast<int>(digits.size()) - 3; i > 0; i -= 3) digits.insert(static_cast<std::size_t>(i), ",");
    return digits;
}

VerificationResult finish(VerificationResult r) {
    r.passed = std::isfinite(r.error) && r.error <= r.tolerance && r.samples > 0;
    return r;
}

VerificationResult unable(std::string why) {
    VerificationResult r;
    r.error = inf;
    r.tolerance = 0.0;
    r.note = std::move(why);
    return r;
}

FlowModel modelOf(const Scene& scene) { return buildFlowModel(scene); }

// --- Exact solutions, in polar form about each element ------------------------------------

struct Polar {
    double r;
    double theta;
};

Polar polarAbout(Vec2 p, Vec2 centre) {
    const Vec2 d = p - centre;
    return {length(d), std::atan2(d.y, d.x)};
}

// Cartesian velocity from its radial and tangential components at polar angle theta.
Vec2 fromPolar(double vr, double vt, double theta) {
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    return {vr * c - vt * s, vr * s + vt * c};
}

Vec2 exactUniform(double speed, double alpha) { return {speed * std::cos(alpha), speed * std::sin(alpha)}; }

// V_r = Λ / 2πr, V_θ = 0.
Vec2 exactSource(double strength, Vec2 centre, Vec2 p) {
    const auto [r, theta] = polarAbout(p, centre);
    return fromPolar(strength / (2.0 * pi * r), 0.0, theta);
}

// V_r = 0, V_θ = -Γ / 2πr: clockwise for Γ > 0.
Vec2 exactVortex(double circulation, Vec2 centre, Vec2 p) {
    const auto [r, theta] = polarAbout(p, centre);
    return fromPolar(0.0, -circulation / (2.0 * pi * r), theta);
}

// V_r = -μ cos(θ - β) / r², V_θ = -μ sin(θ - β) / r² (Anderson's doublet with μ = κ / 2π).
Vec2 exactDoublet(double strength, double beta, Vec2 centre, Vec2 p) {
    const auto [r, theta] = polarAbout(p, centre);
    return fromPolar(-strength * std::cos(theta - beta) / (r * r), -strength * std::sin(theta - beta) / (r * r), theta);
}

// Cylinder of radius a in a stream V∞ at angle α, with circulation Γ:
// V_r = V∞ (1 - a²/r²) cos(θ - α), V_θ = -V∞ (1 + a²/r²) sin(θ - α) - Γ / 2πr.
Vec2 exactCylinder(double speed, double alpha, double radius, double circulation, Vec2 centre, Vec2 p) {
    const auto [r, theta] = polarAbout(p, centre);
    const double a2 = radius * radius / (r * r);
    return fromPolar(speed * (1.0 - a2) * std::cos(theta - alpha),
                     -speed * (1.0 + a2) * std::sin(theta - alpha) - circulation / (2.0 * pi * r), theta);
}

// Flow in the wedge between the walls θ1 and θ2, n = π / (θ2 - θ1), with θ measured from θ1:
// V_r = n V r^(n-1) cos(nθ), V_θ = -n V r^(n-1) sin(nθ). nullopt outside the wedge.
std::optional<Vec2> exactSector(double strength, double wall1, double wall2, Vec2 vertex, Vec2 p) {
    constexpr double edge = 1e-9; // points this close to a wall count as on it
    const double opening = wall2 - wall1;
    const auto [r, theta] = polarAbout(p, vertex);
    double local = std::fmod(theta - wall1, 2.0 * pi);
    if (local < 0.0) local += 2.0 * pi;
    if (local > 2.0 * pi - edge && opening < 2.0 * pi - edge) local = 0.0;
    if (local > opening + edge) return std::nullopt;
    local = std::min(local, opening);
    const double n = pi / opening;
    const double magnitude = n * strength * std::pow(r, n - 1.0);
    return fromPolar(magnitude * std::cos(n * local), -magnitude * std::sin(n * local), theta);
}

// Cylinder of radius a at the origin in a stream V∞ along +x, with circulation Γ:
// ψ = V∞ (r − a²/r) sin θ + (Γ / 2π) ln(r / a), which is 0 on the cylinder.
double cylinderPsi(double speed, double radius, double circulation, Vec2 p) {
    const auto [r, theta] = polarAbout(p, {});
    return speed * (r - radius * radius / r) * std::sin(theta) + circulation / (2.0 * pi) * std::log(r / radius);
}

// Half-height h of the Rankine oval with V∞ = 1, b = 1, Λ = 2π: V∞ h = (Λ / 2π)(π − 2 atan(h / b)).
double ovalHalfHeight() {
    double lo = 0.0;
    double hi = 5.0;
    for (int k = 0; k < 200; ++k) {
        const double mid = 0.5 * (lo + hi);
        (mid - (pi - 2.0 * std::atan(mid)) > 0.0 ? hi : lo) = mid;
    }
    return 0.5 * (lo + hi);
}

bool at(Vec2 p, Vec2 singularity) { return length(p - singularity) < 1e-9; }

// --- Comparing a whole field with theory --------------------------------------------------

struct Exact {
    enum Kind { Value, Singular, Outside } kind = Value;
    Vec2 v;
};
using ExactField = std::function<Exact(Vec2)>;

Exact value(Vec2 v) { return {Exact::Value, v}; }
constexpr Exact k_singular{Exact::Singular, {}};
constexpr Exact k_outside{Exact::Outside, {}};

struct GridStats {
    double error = 0.0;
    int compared = 0;
    int singular = 0;   // nodes on a singular point, skipped
    int masked = 0;     // nodes outside the flow domain that the program masks, as it should
    int maskErrors = 0; // nodes masked when they should not be, or the other way round

    GridStats& operator+=(const GridStats& o) {
        error = std::max(error, o.error);
        compared += o.compared;
        singular += o.singular;
        masked += o.masked;
        maskErrors += o.maskErrors;
        return *this;
    }
};

// Every node of the grid over `view`: |V_program - V_exact| relative to |V_exact|, but to no less
// than the median exact speed, so that nodes where the flow nearly stops (stagnation points) do not
// turn rounding into large relative errors. The same measure as tools/report_figures.py.
GridStats compareGrid(const FlowField& field, const ExactField& exact, const Rect& view = k_view) {
    GridStats stats;
    std::vector<std::pair<Vec2, Vec2>> pairs; // program, exact
    pairs.reserve(static_cast<std::size_t>(k_gridColumns * k_gridRows));
    for (int j = 0; j < k_gridRows; ++j) {
        for (int i = 0; i < k_gridColumns; ++i) {
            const Vec2 p{view.min.x + view.size().x * i / (k_gridColumns - 1),
                         view.min.y + view.size().y * j / (k_gridRows - 1)};
            const Exact e = exact(p);
            if (e.kind == Exact::Singular) {
                ++stats.singular;
                continue;
            }
            const bool inside = field.inDomain(p);
            if (e.kind == Exact::Outside) {
                ++(inside ? stats.maskErrors : stats.masked);
                continue;
            }
            if (!inside) {
                ++stats.maskErrors;
                continue;
            }
            pairs.emplace_back(field.velocity(p), e.v);
        }
    }
    if (pairs.empty()) {
        stats.error = inf;
        return stats;
    }
    std::vector<double> speeds;
    speeds.reserve(pairs.size());
    for (const auto& [program, ex] : pairs) speeds.push_back(length(ex));
    const auto middle = speeds.begin() + static_cast<std::ptrdiff_t>(speeds.size() / 2);
    std::nth_element(speeds.begin(), middle, speeds.end());
    const double median = *middle;
    for (const auto& [program, ex] : pairs) {
        const double e = length(program - ex) / std::max(length(ex), median);
        stats.error = std::isfinite(e) ? std::max(stats.error, e) : inf;
    }
    stats.compared = static_cast<int>(pairs.size());
    if (stats.maskErrors > 0) stats.error = inf;
    return stats;
}

VerificationResult gridResult(const GridStats& g, std::vector<Comparison> comparisons, int fields = 1) {
    VerificationResult r;
    r.error = g.error;
    r.tolerance = k_exactTolerance;
    r.errorLabel = "max |V − V_exact| / |V_exact| over every grid node";
    r.samples = g.compared;
    r.comparisons = std::move(comparisons);
    r.note = fields == 1 ? std::format("{} grid nodes compared over the view", count(g.compared))
                         : std::format("{} grid nodes compared over the view, in {} flows", count(g.compared), fields);
    if (g.singular > 0) r.note += std::format("; {} on a singular point skipped", g.singular);
    if (g.masked > 0) r.note += std::format("; {} outside the flow domain, correctly masked", count(g.masked));
    if (g.maskErrors > 0) r.note += std::format("; {} masked wrongly", count(g.maskErrors));
    r.note += ".";
    return finish(std::move(r));
}

Comparison velocityAt(const FlowField& field, Vec2 p, Vec2 exact) {
    return {"V at " + point(p), vec(field.velocity(p)), vec(exact)};
}

// --- Scenes -------------------------------------------------------------------------------

Scene withPlot(Scene scene, void (*setup)(PlotSettings&)) {
    setup(scene.plot);
    return scene;
}

// A bit of everything, for the checks that hold for any flow.
Scene mixedScene() {
    return {.elements = {UniformFlow{.speed = 1.0, .angle = 15.0 * deg},
                         SourceSink{.strength = 2.5, .position = {-1.5, 0.6}},
                         SourceSink{.strength = -1.8, .position = {1.6, -0.4}},
                         Vortex{.circulation = 3.0, .position = {0.2, 1.1}},
                         Doublet{.strength = 0.3, .position = {0.4, -1.2}, .angle = 35.0 * deg}}};
}

Scene uniformScene() { return {.elements = {UniformFlow{.speed = 1.0, .angle = 45.0 * deg}}}; }
Scene sourceSinkScene() {
    return {.elements = {SourceSink{.strength = 2.0 * pi, .position = {-1.5, 0.0}},
                         SourceSink{.strength = -2.0 * pi, .position = {1.5, 0.0}}}};
}
Scene vortexScene() { return {.elements = {Vortex{.circulation = 2.0 * pi}}}; }
Scene doubletScene() { return {.elements = {Doublet{.strength = 1.0}}}; }
Scene cornerScene() {
    return withPlot({.elements = {SectorFlow{.strength = 1.0, .wall1 = 0.0, .wall2 = 90.0 * deg}}},
                    [](PlotSettings& p) { p.view = {{-1.0, -0.6}, {4.0, 2.5}}; });
}
Scene plateScene() { return {.elements = {SectorFlow{.strength = 1.0, .wall1 = 0.0, .wall2 = 2.0 * pi}}}; }
Scene cylinderScene() { return {.elements = {UniformFlow{}, Doublet{.strength = 1.0}}}; }
Scene liftingCylinderScene() {
    return {.elements = {UniformFlow{}, Doublet{.strength = 1.0}, Vortex{.circulation = 2.0 * pi}}};
}
Scene rankineOvalScene() {
    return {.elements = {UniformFlow{}, SourceSink{.strength = 2.0 * pi, .position = {-1.0, 0.0}},
                         SourceSink{.strength = -2.0 * pi, .position = {1.0, 0.0}}}};
}
Scene halfBodyScene() {
    return {.elements = {UniformFlow{}, SourceSink{.strength = 2.0 * pi, .position = {-1.0, 0.0}}}};
}
Scene sourceNearWallScene() {
    return {.elements = {SourceSink{.strength = 2.0 * pi, .position = {0.0, 1.0}}, SymmetryWall{}}};
}
Scene inclinedWallScene() {
    return {.elements = {UniformFlow{.speed = 1.0, .angle = 10.0 * deg},
                         SourceSink{.strength = 4.0, .position = {1.0, 1.0}},
                         Vortex{.circulation = -3.0, .position = {-1.5, 0.5}},
                         Doublet{.strength = 0.7, .position = {0.5, 2.0}, .angle = 70.0 * deg},
                         SymmetryWall{.point = {0.3, -0.8}, .angle = 25.0 * deg}}};
}
Scene cornerWallsScene() {
    return {.elements = {Vortex{.circulation = 1.0, .position = {1.0, 2.0}}, SymmetryWall{},
                         SymmetryWall{.angle = 90.0 * deg}}};
}
Scene channelSourceScene() {
    return {.elements = {SourceSink{.strength = 2.0 * pi, .position = {0.0, 0.4}}, SymmetryWall{.point = {0.0, -1.0}},
                         SymmetryWall{.point = {0.0, 1.0}}}};
}
Scene inclinedChannelScene() {
    const SymmetryWall lower{.point = {0.3, -0.8}, .angle = 25.0 * deg};
    return {.elements = {UniformFlow{.speed = 1.0, .angle = 10.0 * deg},
                         SourceSink{.strength = 4.0, .position = {1.0, 0.0}},
                         Vortex{.circulation = -3.0, .position = {-1.5, -0.5}},
                         Doublet{.strength = 0.7, .position = {0.5, 0.9}, .angle = 70.0 * deg}, lower,
                         SymmetryWall{.point = lower.point + lower.normal() * 2.0, .angle = lower.angle}}};
}
Scene surfacePressureScene() {
    return withPlot(liftingCylinderScene(), [](PlotSettings& p) { p.showSurfacePressure = true; });
}
Scene dragScene() { return withPlot(rankineOvalScene(), [](PlotSettings& p) { p.showSurfacePressure = true; }); }
Scene liftScene() {
    return withPlot({.elements = {UniformFlow{.speed = 2.0, .angle = 20.0 * deg},
                                  Doublet{.strength = 2.0, .position = {0.5, -0.3}, .angle = 20.0 * deg},
                                  Vortex{.circulation = 2.0 * pi, .position = {0.5, -0.3}}}},
                    [](PlotSettings& p) { p.showSurfacePressure = true; });
}
Scene equipotentialScene() {
    return withPlot({.elements = {UniformFlow{.speed = 1.0, .angle = 20.0 * deg},
                                  Doublet{.strength = 1.0, .angle = 20.0 * deg}}},
                    [](PlotSettings& p) { p.showEquipotentials = true; });
}

Scene offBodyStagnationScene() {
    return {.elements = {UniformFlow{}, Doublet{.strength = 1.0}, Vortex{.circulation = 6.0 * pi}}};
}
Scene twoSourcesScene() {
    return {.elements = {SourceSink{.strength = 2.0 * pi, .position = {-1.0, 0.0}},
                         SourceSink{.strength = 2.0 * pi, .position = {1.0, 0.0}}}};
}
const Vec2 k_wedgeElement{1.5 * std::cos(20.0 * deg), 1.5 * std::sin(20.0 * deg)};
Scene wedgeSourceScene() {
    return {.elements = {SourceSink{.strength = 2.0 * pi, .position = k_wedgeElement}, SymmetryWall{},
                         SymmetryWall{.angle = 60.0 * deg}}};
}
Scene channelVortexScene() {
    return {.elements = {Vortex{.circulation = 2.0 * pi}, SymmetryWall{.point = {0.0, -1.0}},
                         SymmetryWall{.point = {0.0, 1.0}}}};
}
Scene halfBodyPressureScene() {
    return withPlot(halfBodyScene(), [](PlotSettings& p) { p.showSurfacePressure = true; });
}
Scene obliqueSectorScene() {
    return {.elements = {SectorFlow{.strength = 1.3, .vertex = {0.31, -0.97}, .wall1 = 30.0 * deg, .wall2 = 120.0 * deg}}};
}
Scene radialScene() { return {.elements = {SourceSink{.strength = 2.0 * pi, .position = {0.3, -0.2}}}}; }
Scene dyeScene() {
    Scene scene = liftingCylinderScene();
    scene.rakes.push_back({.start = {-3.0, -1.5}, .end = {-3.0, 1.5}, .emitters = 5});
    return scene;
}

// Exact sum of the elements of mixedScene().
Vec2 exactMixed(Vec2 p) {
    return exactUniform(1.0, 15.0 * deg) + exactSource(2.5, {-1.5, 0.6}, p) + exactSource(-1.8, {1.6, -0.4}, p) +
           exactVortex(3.0, {0.2, 1.1}, p) + exactDoublet(0.3, 35.0 * deg, {0.4, -1.2}, p);
}

// --- Foundations ------------------------------------------------------------------------

// Points spread over the view, at least 0.25 from the singular points of mixedScene(). A fixed
// sequence, so every run checks the same points.
std::vector<Vec2> scatteredPoints(int count) {
    const Vec2 singular[] = {{-1.5, 0.6}, {1.6, -0.4}, {0.2, 1.1}, {0.4, -1.2}};
    std::uint32_t state = 20'240'601u;
    const auto uniform = [&state] {
        state = state * 1'664'525u + 1'013'904'223u;
        return static_cast<double>(state >> 8) / static_cast<double>(1u << 24);
    };
    std::vector<Vec2> points;
    while (static_cast<int>(points.size()) < count) {
        const Vec2 p{k_view.min.x + k_view.size().x * uniform(), k_view.min.y + k_view.size().y * uniform()};
        if (std::ranges::none_of(singular, [p](Vec2 s) { return length(p - s) < 0.25; })) points.push_back(p);
    }
    return points;
}

VerificationResult laplace() {
    const FlowField field = modelOf(mixedScene()).field;
    constexpr double h = 1e-5;
    double divergence = 0.0;
    double curl = 0.0;
    int samples = 0;
    for (const Vec2 p : scatteredPoints(2000)) {
        const Vec2 dx = (field.velocity(p + Vec2{h, 0.0}) - field.velocity(p - Vec2{h, 0.0})) / (2.0 * h);
        const Vec2 dy = (field.velocity(p + Vec2{0.0, h}) - field.velocity(p - Vec2{0.0, h})) / (2.0 * h);
        const double gradient = std::sqrt(lengthSquared(dx) + lengthSquared(dy));
        divergence = std::max(divergence, std::abs(dx.x + dy.y) / gradient);
        curl = std::max(curl, std::abs(dx.y - dy.x) / gradient);
        ++samples;
    }
    VerificationResult r;
    r.error = std::max(divergence, curl);
    r.tolerance = 1e-6;
    r.errorLabel = "max |div V| and |curl V| relative to the velocity gradient, by central differences";
    r.samples = samples;
    r.comparisons = {{"div V = ∂u/∂x + ∂v/∂y (incompressible)", num(divergence), "0"},
                     {"curl V = ∂v/∂x − ∂u/∂y (irrotational)", num(curl), "0"}};
    r.note = std::format("{} points spread over the view, at least 0.25 from every singular point. "
                         "The tolerance allows for the finite differences (step 10⁻⁵).",
                         count(samples));
    return finish(std::move(r));
}

VerificationResult cauchyRiemann() {
    const Vec2 c{0.3, -0.2};
    const FlowElement elements[] = {
        UniformFlow{.speed = 1.7, .angle = 25.0 * deg},
        SourceSink{.strength = -2.4, .position = c},
        Vortex{.circulation = 3.1, .position = c},
        Doublet{.strength = 0.8, .position = c, .angle = 40.0 * deg},
        SectorFlow{.strength = 1.2, .vertex = c, .wall1 = 30.0 * deg, .wall2 = 160.0 * deg},
        SectorFlow{.strength = 0.9, .vertex = c, .wall1 = -45.0 * deg, .wall2 = 225.0 * deg},
    };
    const char* names[] = {"Uniform flow", "Sink", "Vortex", "Doublet", "Sector, n ≈ 1.38", "Sector, n = 2/3"};
    // Inside every wedge and away from the branch cut of log along -x from c.
    const Vec2 points[] = {{1.1, 0.9}, {0.5, 1.4}, {-0.4, 0.8}, {1.6, 0.1}};
    constexpr double h = 1e-6;
    VerificationResult r;
    for (std::size_t k = 0; k < std::size(elements); ++k) {
        double worst = 0.0;
        for (const Vec2 p : points) {
            const Complex w = complexVelocity(elements[k], p);
            const Complex dFdx =
                (complexPotential(elements[k], p + Vec2{h, 0.0}) - complexPotential(elements[k], p - Vec2{h, 0.0})) / (2.0 * h);
            const Complex dFdy =
                (complexPotential(elements[k], p + Vec2{0.0, h}) - complexPotential(elements[k], p - Vec2{0.0, h})) / (2.0 * h);
            // dF/dz = ∂F/∂x = -i ∂F/∂y, i.e. u = ∂φ/∂x = ∂ψ/∂y and v = ∂φ/∂y = -∂ψ/∂x.
            const double scale = 1.0 + std::abs(w);
            worst = std::max({worst, std::abs(dFdx - w) / scale, std::abs(dFdy - Complex{0.0, 1.0} * w) / scale});
            r.samples += 2;
        }
        r.error = std::max(r.error, worst);
        r.comparisons.push_back({std::string(names[k]) + ": |dF/dz − W|", num(worst), "0"});
    }
    r.tolerance = 1e-6;
    r.errorLabel = "max |dF/dz − W| / (1 + |W|), dF/dz by central differences";
    r.note = "Every element type, at four points each, in x and in y. The tolerance allows for the finite "
             "differences (step 10⁻⁶).";
    return finish(std::move(r));
}

VerificationResult fluxAndCirculation() {
    const FlowField field = modelOf(mixedScene()).field;
    struct Contour {
        const char* name;
        Vec2 centre;
        double radius;
        double flux;        // the source strengths inside
        double circulation; // counter-clockwise, so minus the (clockwise) vortex strengths inside
    };
    const Contour contours[] = {
        {"Circle around the source", {-1.5, 0.6}, 0.4, 2.5, 0.0},
        {"Circle around the sink", {1.6, -0.4}, 0.4, -1.8, 0.0},
        {"Circle around the vortex", {0.2, 1.1}, 0.4, 0.0, -3.0},
        {"Circle around the doublet", {0.4, -1.2}, 0.4, 0.0, 0.0},
        {"Circle around everything", {0.0, 0.0}, 3.0, 0.7, -3.0},
        {"Circle around nothing", {2.6, 1.6}, 0.5, 0.0, 0.0},
    };
    constexpr int n = 1024; // the trapezoidal rule on a circle converges geometrically
    VerificationResult r;
    for (const Contour& c : contours) {
        double flux = 0.0;
        double circulation = 0.0;
        const double ds = 2.0 * pi * c.radius / n;
        for (int k = 0; k < n; ++k) {
            const double t = 2.0 * pi * k / n;
            const Vec2 normal{std::cos(t), std::sin(t)};
            const Vec2 v = field.velocity(c.centre + normal * c.radius);
            flux += dot(v, normal) * ds;
            circulation += dot(v, perp(normal)) * ds;
        }
        r.samples += n;
        r.error = std::max({r.error, std::abs(flux - c.flux) / 3.0, std::abs(circulation - c.circulation) / 3.0});
        const auto zero = [](double x) { return std::abs(x) < 1e-12 ? 0.0 : x; };
        r.comparisons.push_back({c.name, std::format("Q = {}, Γ = {}", num(zero(flux)), num(zero(circulation))),
                                 std::format("Q = {}, Γ = {}", num(c.flux), num(c.circulation))});
    }
    r.tolerance = 1e-9;
    r.errorLabel = "max |∫ v·n ds − ΣΛ| and |∫ v·t ds + ΣΓ| around each circle, relative to the largest strength";
    r.note = "Q = ∫ v·n ds is the volume flow out of the circle, Γ = ∫ v·t ds the circulation around it "
             "(counter-clockwise, so a clockwise vortex Γ = 3 gives −3). 1024 points per circle.";
    return finish(std::move(r));
}

VerificationResult superposition() {
    const FlowField field = modelOf(mixedScene()).field;
    const GridStats g = compareGrid(field, [](Vec2 p) {
        const bool singular = at(p, {-1.5, 0.6}) || at(p, {1.6, -0.4}) || at(p, {0.2, 1.1}) || at(p, {0.4, -1.2});
        return singular ? k_singular : value(exactMixed(p));
    });
    return gridResult(g, {velocityAt(field, {0.0, 0.0}, exactMixed({0.0, 0.0})),
                          velocityAt(field, {-3.0, -2.0}, exactMixed({-3.0, -2.0}))});
}

VerificationResult singularPoints() {
    using Term = std::function<Vec2(Vec2)>;
    const Term terms[] = {[](Vec2) { return exactUniform(1.0, 15.0 * deg); },
                          [](Vec2 p) { return exactSource(2.5, {-1.5, 0.6}, p); },
                          [](Vec2 p) { return exactSource(-1.8, {1.6, -0.4}, p); },
                          [](Vec2 p) { return exactVortex(3.0, {0.2, 1.1}, p); },
                          [](Vec2 p) { return exactDoublet(0.3, 35.0 * deg, {0.4, -1.2}, p); }};
    struct Centre {
        const char* name;
        Vec2 position;
        std::size_t term; // the element centred here
    };
    const Centre centres[] = {{"Source", {-1.5, 0.6}, 1}, {"Sink", {1.6, -0.4}, 2}, {"Vortex", {0.2, 1.1}, 3}, {"Doublet", {0.4, -1.2}, 4}};
    const FlowField field = modelOf(mixedScene()).field;
    VerificationResult r;
    for (const Centre& c : centres) {
        Vec2 others{};
        for (std::size_t k = 0; k < std::size(terms); ++k) {
            if (k != c.term) others += terms[k](c.position);
        }
        const Vec2 v = field.velocity(c.position);
        const double e = length(v - others) / length(others);
        r.error = std::isfinite(e) ? std::max(r.error, e) : inf;
        ++r.samples;
        r.comparisons.push_back({std::format("{}: V at its own centre", c.name), vec(v), vec(others) + " (the others)"});
    }
    // A plate edge is singular too: only the uniform flow is left at the vertex.
    const FlowField plate = modelOf({.elements = {UniformFlow{.speed = 1.0, .angle = 15.0 * deg},
                                                  SectorFlow{.strength = 1.0, .wall1 = 0.0, .wall2 = 2.0 * pi}}})
                                .field;
    const Vec2 v = plate.velocity({});
    const double e = length(v - exactUniform(1.0, 15.0 * deg));
    r.error = std::isfinite(e) ? std::max(r.error, e) : inf;
    ++r.samples;
    r.comparisons.push_back({"Plate edge: V at the vertex", vec(v), vec(exactUniform(1.0, 15.0 * deg)) + " (the stream)"});
    r.tolerance = k_exactTolerance;
    r.errorLabel = "max |V − V_others| / |V_others| at each singular point";
    r.note = "An element's own term is left out at its centre, so the velocity stays finite there.";
    return finish(std::move(r));
}

Vec2 turned(Vec2 v, double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    return {c * v.x - s * v.y, s * v.x + c * v.y};
}

Vec2 mirroredInX(Vec2 v) { return {v.x, -v.y}; }

// `element` rotated about the origin.
SceneElement turnedElement(SceneElement element, double a) {
    std::visit(
        [a](auto& e) {
            using T = std::decay_t<decltype(e)>;
            if constexpr (std::is_same_v<T, UniformFlow>) {
                e.angle += a;
            } else if constexpr (std::is_same_v<T, SourceSink> || std::is_same_v<T, Vortex>) {
                e.position = turned(e.position, a);
            } else if constexpr (std::is_same_v<T, Doublet>) {
                e.position = turned(e.position, a);
                e.angle += a;
            } else if constexpr (std::is_same_v<T, SectorFlow>) {
                e.vertex = turned(e.vertex, a);
                e.wall1 += a;
                e.wall2 += a;
            } else {
                e.point = turned(e.point, a);
                e.angle += a;
            }
        },
        element);
    return element;
}

// `element` reflected in the x axis. A vortex turns the other way, and a sector's flow runs in along
// what was its first wall, so both change sign.
SceneElement mirroredElement(SceneElement element) {
    std::visit(
        [](auto& e) {
            using T = std::decay_t<decltype(e)>;
            if constexpr (std::is_same_v<T, UniformFlow>) {
                e.angle = -e.angle;
            } else if constexpr (std::is_same_v<T, SourceSink>) {
                e.position = mirroredInX(e.position);
            } else if constexpr (std::is_same_v<T, Vortex>) {
                e.position = mirroredInX(e.position);
                e.circulation = -e.circulation;
            } else if constexpr (std::is_same_v<T, Doublet>) {
                e.position = mirroredInX(e.position);
                e.angle = -e.angle;
            } else if constexpr (std::is_same_v<T, SectorFlow>) {
                e.vertex = mirroredInX(e.vertex);
                e.strength = -e.strength;
                e.wall1 = -std::exchange(e.wall2, -e.wall1);
            } else {
                e.point = mirroredInX(e.point);
                e.angle = -e.angle;
            }
        },
        element);
    return element;
}

// The flow of `scene` with every element mapped, against the flow of `scene` moved the same way.
// Both sides round the positions differently, which a doublet's 1/r² turns into an error growing as
// 1/r, so nodes within 0.01 of a singular point are left out.
GridStats compareMoved(const Scene& scene, const std::function<SceneElement(SceneElement)>& move,
                       const std::function<Vec2(Vec2)>& forward, const std::function<Vec2(Vec2)>& back) {
    const FlowField original = modelOf(scene).field;
    Scene moved = scene;
    for (SceneElement& element : moved.elements) element = move(std::move(element));
    return compareGrid(modelOf(moved).field, [&](Vec2 p) {
        const Vec2 q = back(p);
        if (original.nearSingularity(q, 0.01)) return k_singular;
        return original.inDomain(q) ? value(forward(original.velocity(q))) : k_outside;
    });
}

VerificationResult symmetries() {
    constexpr double a = 33.0 * deg;
    const Scene scenes[] = {inclinedWallScene(), inclinedChannelScene(), obliqueSectorScene()};
    GridStats g;
    for (const Scene& scene : scenes) {
        g += compareMoved(
            scene, [](SceneElement e) { return turnedElement(std::move(e), a); }, [](Vec2 v) { return turned(v, a); },
            [](Vec2 v) { return turned(v, -a); });
        g += compareMoved(scene, [](SceneElement e) { return mirroredElement(std::move(e)); }, mirroredInX, mirroredInX);
    }
    VerificationResult r = gridResult(g, {{"Flows compared", "6", "3 scenes, turned 33° and mirrored in the x axis"}}, 6);
    r.note += " Nodes within 0.01 of a singular point count as on it: the two sides round their positions differently, "
              "and next to a doublet that is all the difference there is.";
    return r;
}

// F(b) − F(a) = ∫ W dz along the straight line, from the velocity the plot uses (Simpson's rule).
Complex integrateVelocity(const FlowField& field, Vec2 a, Vec2 b) {
    constexpr int n = 4000;
    const auto w = [&](int k) {
        const Vec2 v = field.velocity(a + (b - a) * (static_cast<double>(k) / n));
        return Complex{v.x, -v.y};
    };
    Complex sum = w(0) + w(n);
    for (int k = 1; k < n; ++k) sum += w(k) * (k % 2 == 1 ? 4.0 : 2.0);
    return sum * toComplex(b - a) / (3.0 * n);
}

VerificationResult lineIntegrals() {
    struct Path {
        const FlowField* field;
        const char* name;
        Vec2 from;
        Vec2 to;
    };
    const FlowField mixed = modelOf(mixedScene()).field;
    const FlowField channel = modelOf(channelSourceScene()).field;
    // None crosses a branch cut of φ or ψ, which run along −x from a source or a vortex, and along
    // the channel from a channel's source.
    const Path paths[] = {{&mixed, "Above every element", {-3.0, 1.5}, {3.0, 1.6}},
                          {&mixed, "Below every element", {-3.0, -2.0}, {3.0, -1.5}},
                          {&mixed, "Right of every element", {2.0, -2.0}, {2.5, 2.0}},
                          {&channel, "Along a channel, below the source", {-3.0, -0.8}, {3.0, 0.2}},
                          {&channel, "Along a channel, above the source", {-3.0, 0.7}, {3.0, 0.9}}};
    VerificationResult r;
    for (const Path& path : paths) {
        const Complex program = path.field->complexPotential(path.to) - path.field->complexPotential(path.from);
        const Complex integral = integrateVelocity(*path.field, path.from, path.to);
        r.error = std::max(r.error, std::abs(program - integral) / 3.0);
        r.samples += 4001;
        r.comparisons.push_back({std::format("{}: Δφ, Δψ", path.name), std::format("{:.9f}, {:.9f}", program.real(), program.imag()),
                                 std::format("{:.9f}, {:.9f} (∫ V)", integral.real(), integral.imag())});
    }
    r.tolerance = 1e-9;
    r.errorLabel = "max |ΔF − ∫ W dz| along each path, relative to the largest strength";
    r.note = "Simpson's rule with 4,000 intervals, along paths that keep clear of the branch cuts of log z.";
    return finish(std::move(r));
}

VerificationResult bernoulliEuler() {
    const FlowModel model = modelOf(mixedScene());
    const FlowField& field = model.field;
    const double stream = model.freeStreamSpeed;
    constexpr double h = 1e-5;
    const auto cp = [&](Vec2 p) { return pressureCoefficient(length(field.velocity(p)), stream); };
    VerificationResult r;
    for (const Vec2 p : scatteredPoints(2000)) {
        const Vec2 ex{h, 0.0};
        const Vec2 ey{0.0, h};
        const Vec2 dx = (field.velocity(p + ex) - field.velocity(p - ex)) / (2.0 * h);
        const Vec2 dy = (field.velocity(p + ey) - field.velocity(p - ey)) / (2.0 * h);
        const Vec2 v = field.velocity(p);
        const Vec2 acceleration = dx * v.x + dy * v.y;
        // −∇p / ρ with p − p∞ = ½ ρ V∞² Cp.
        const Vec2 pressure = Vec2{cp(p + ex) - cp(p - ex), cp(p + ey) - cp(p - ey)} * (-0.25 * stream * stream / h);
        const double scale = length(v) * std::sqrt(lengthSquared(dx) + lengthSquared(dy));
        r.error = std::max(r.error, length(acceleration - pressure) / scale);
        ++r.samples;
    }
    r.tolerance = 1e-6;
    r.errorLabel = "max |(V·∇)V + ∇p/ρ| relative to |V| |∇V|, by central differences";
    r.comparisons = {{"max |(V·∇)V + ∇p/ρ|, relative", num(r.error), "0"}};
    r.note = std::format("{} points spread over the view, with p from the Cp the program shows. The tolerance allows for "
                         "the finite differences (step 10⁻⁵).",
                         count(r.samples));
    return finish(std::move(r));
}

// --- Elementary flows -------------------------------------------------------------------

VerificationResult uniform() {
    const FlowField field = modelOf(uniformScene()).field;
    const FlowField level = modelOf({.elements = {UniformFlow{.speed = 1.0, .angle = 0.0}}}).field;
    const Vec2 v = exactUniform(1.0, 45.0 * deg);
    const Vec2 along = exactUniform(1.0, 0.0);
    GridStats g = compareGrid(field, [v](Vec2) { return value(v); });
    g += compareGrid(level, [along](Vec2) { return value(along); });
    return gridResult(g, {velocityAt(field, {3.0, -2.0}, v), velocityAt(level, {3.0, -2.0}, along)}, 2);
}

VerificationResult sourceSink() {
    const FlowField source = modelOf({.elements = {SourceSink{.strength = 2.0 * pi}}}).field;
    const FlowField sink = modelOf({.elements = {SourceSink{.strength = -2.0 * pi, .position = {1.0, 0.5}}}}).field;
    const auto exactSourceAt = [](Vec2 p) { return at(p, {}) ? k_singular : value(exactSource(2.0 * pi, {}, p)); };
    const auto exactSinkAt = [](Vec2 p) {
        return at(p, {1.0, 0.5}) ? k_singular : value(exactSource(-2.0 * pi, {1.0, 0.5}, p));
    };
    GridStats g = compareGrid(source, exactSourceAt);
    g += compareGrid(sink, exactSinkAt);
    return gridResult(g,
                      {{"Source Λ = 2π: V at (1, 0)", vec(source.velocity({1.0, 0.0})), vec(exactSource(2.0 * pi, {}, {1.0, 0.0}))},
                       {"Sink Λ = −2π: V at (2, 0.5)", vec(sink.velocity({2.0, 0.5})),
                        vec(exactSource(-2.0 * pi, {1.0, 0.5}, {2.0, 0.5}))}},
                      2);
}

VerificationResult vortex() {
    const FlowField clockwise = modelOf(vortexScene()).field;
    const FlowField anticlockwise = modelOf({.elements = {Vortex{.circulation = -3.0, .position = {1.0, -0.5}}}}).field;
    GridStats g = compareGrid(clockwise, [](Vec2 p) { return at(p, {}) ? k_singular : value(exactVortex(2.0 * pi, {}, p)); });
    g += compareGrid(anticlockwise, [](Vec2 p) {
        return at(p, {1.0, -0.5}) ? k_singular : value(exactVortex(-3.0, {1.0, -0.5}, p));
    });
    return gridResult(g,
                      {{"Γ = 2π: V at (1, 0)", vec(clockwise.velocity({1.0, 0.0})), vec(exactVortex(2.0 * pi, {}, {1.0, 0.0}))},
                       {"Γ = 2π: V at (0, 1)", vec(clockwise.velocity({0.0, 1.0})), vec(exactVortex(2.0 * pi, {}, {0.0, 1.0}))}},
                      2);
}

VerificationResult doublet() {
    const FlowField axial = modelOf(doubletScene()).field;
    const Vec2 c{0.3, -0.2};
    const FlowField rotated = modelOf({.elements = {Doublet{.strength = 0.8, .position = c, .angle = 40.0 * deg}}}).field;
    GridStats g = compareGrid(axial, [](Vec2 p) { return at(p, {}) ? k_singular : value(exactDoublet(1.0, 0.0, {}, p)); });
    g += compareGrid(rotated, [c](Vec2 p) { return at(p, c) ? k_singular : value(exactDoublet(0.8, 40.0 * deg, c, p)); });
    return gridResult(g,
                      {{"μ = 1: V at (1, 0)", vec(axial.velocity({1.0, 0.0})), vec(exactDoublet(1.0, 0.0, {}, {1.0, 0.0}))},
                       {"μ = 1: V at (0, 1)", vec(axial.velocity({0.0, 1.0})), vec(exactDoublet(1.0, 0.0, {}, {0.0, 1.0}))}},
                      2);
}

ExactField sectorField(double strength, double wall1, double wall2, Vec2 vertex) {
    return [=](Vec2 p) {
        if (pi / (wall2 - wall1) < 1.0 && at(p, vertex)) return k_singular;
        const std::optional<Vec2> v = exactSector(strength, wall1, wall2, vertex, p);
        return v ? value(*v) : k_outside;
    };
}

VerificationResult corner() {
    const FlowField square = modelOf(cornerScene()).field;
    const Vec2 vertex{0.5, -1.0};
    const FlowField oblique =
        modelOf({.elements = {SectorFlow{.strength = 1.3, .vertex = vertex, .wall1 = 30.0 * deg, .wall2 = 120.0 * deg}}}).field;
    GridStats g = compareGrid(square, sectorField(1.0, 0.0, 90.0 * deg, {}));
    g += compareGrid(oblique, sectorField(1.3, 30.0 * deg, 120.0 * deg, vertex));
    return gridResult(g,
                      {velocityAt(square, {1.0, 1.0}, *exactSector(1.0, 0.0, 90.0 * deg, {}, {1.0, 1.0})),
                       velocityAt(square, {2.0, 0.5}, *exactSector(1.0, 0.0, 90.0 * deg, {}, {2.0, 0.5}))},
                      2);
}

VerificationResult plateEdge() {
    const FlowField field = modelOf(plateScene()).field;
    const GridStats g = compareGrid(field, sectorField(1.0, 0.0, 2.0 * pi, {}));
    return gridResult(g, {velocityAt(field, {-1.0, 0.0}, *exactSector(1.0, 0.0, 2.0 * pi, {}, {-1.0, 0.0})),
                          velocityAt(field, {1.0, 1.0}, *exactSector(1.0, 0.0, 2.0 * pi, {}, {1.0, 1.0})),
                          velocityAt(field, {1.0, -1.0}, *exactSector(1.0, 0.0, 2.0 * pi, {}, {1.0, -1.0}))});
}

// --- Superposition ----------------------------------------------------------------------

VerificationResult cylinder() {
    const FlowField straight = modelOf(cylinderScene()).field;
    const Vec2 c{1.0, -0.5};
    const double alpha = 30.0 * deg;
    const FlowField tilted = modelOf({.elements = {UniformFlow{.speed = 2.0, .angle = alpha},
                                                   Doublet{.strength = 2.0 * 0.49, .position = c, .angle = alpha}}})
                                 .field;
    GridStats g = compareGrid(straight, [](Vec2 p) {
        return at(p, {}) ? k_singular : value(exactCylinder(1.0, 0.0, 1.0, 0.0, {}, p));
    });
    g += compareGrid(tilted, [=](Vec2 p) { return at(p, c) ? k_singular : value(exactCylinder(2.0, alpha, 0.7, 0.0, c, p)); });

    // The circle r = a is a streamline: no flow through it.
    double radial = 0.0;
    for (int k = 0; k < 360; ++k) {
        const Vec2 n{std::cos(k * deg), std::sin(k * deg)};
        radial = std::max({radial, std::abs(dot(straight.velocity(n), n)), std::abs(dot(tilted.velocity(c + n * 0.7), n)) / 2.0});
    }
    VerificationResult r =
        gridResult(g,
                   {velocityAt(straight, {0.0, 1.0}, exactCylinder(1.0, 0.0, 1.0, 0.0, {}, {0.0, 1.0})),
                    velocityAt(straight, {-1.0, 0.0}, exactCylinder(1.0, 0.0, 1.0, 0.0, {}, {-1.0, 0.0})),
                    {"max |V_r| / V∞ on the surface r = a", num(radial), "0"}},
                   2);
    r.error = std::max(r.error, radial);
    return finish(std::move(r));
}

VerificationResult liftingCylinder() {
    const FlowField field = modelOf(liftingCylinderScene()).field;
    const auto exact = [](Vec2 p) { return exactCylinder(1.0, 0.0, 1.0, 2.0 * pi, {}, p); };
    const GridStats g = compareGrid(field, [&](Vec2 p) { return at(p, {}) ? k_singular : value(exact(p)); });
    return gridResult(g, {velocityAt(field, {0.0, 1.0}, exact({0.0, 1.0})), velocityAt(field, {0.0, -1.0}, exact({0.0, -1.0})),
                          velocityAt(field, {2.0, 1.0}, exact({2.0, 1.0}))});
}

VerificationResult cylinderStagnation() {
    struct Case {
        double circulation;
        Vec2 expected[2]; // sorted by x
    };
    // sin θ = -Γ / 4π a V∞: θ = 180° and 0° without circulation, 210° and 330° for Γ = 2π.
    const Case cases[] = {{0.0, {{-1.0, 0.0}, {1.0, 0.0}}},
                          {2.0 * pi, {{std::cos(210.0 * deg), std::sin(210.0 * deg)}, {std::cos(330.0 * deg), std::sin(330.0 * deg)}}}};
    VerificationResult r;
    for (const Case& c : cases) {
        const FlowField field =
            modelOf({.elements = {UniformFlow{}, Doublet{.strength = 1.0}, Vortex{.circulation = c.circulation}}}).field;
        std::vector<Vec2> found = findStagnationPoints(field, {{-3.0, -3.0}, {3.0, 3.0}}, 1.0);
        if (found.size() != 2) return unable(std::format("Found {} stagnation points for Γ = {}, expected 2.", found.size(), num(c.circulation)));
        std::ranges::sort(found, {}, &Vec2::x);
        for (int k = 0; k < 2; ++k) {
            r.error = std::max(r.error, length(found[static_cast<std::size_t>(k)] - c.expected[k]));
            ++r.samples;
            r.comparisons.push_back({std::format("Γ = {}: stagnation point {}", c.circulation == 0.0 ? "0" : "2π", k + 1),
                                     std::format("({:.9f}, {:.9f})", found[static_cast<std::size_t>(k)].x + 0.0,
                                                 found[static_cast<std::size_t>(k)].y + 0.0),
                                     std::format("({:.9f}, {:.9f})", c.expected[k].x + 0.0, c.expected[k].y + 0.0)});
        }
    }
    r.tolerance = 1e-9;
    r.errorLabel = "max distance from the exact stagnation point, in cylinder radii";
    r.note = "Found by Newton's method on W(z) = 0 from a grid of starting points; the exact points lie on "
             "the cylinder where sin θ = −Γ / 4πaV∞.";
    return finish(std::move(r));
}

VerificationResult stagnationOffBody() {
    // With s = Γ / 4πaV∞ > 1 the two roots leave the cylinder along the line θ = −90°, at
    // r = a (s ± √(s² − 1)); at s = 1 they meet at the bottom of it as a double root.
    const auto found = [](double circulation) {
        const FlowField field =
            modelOf({.elements = {UniformFlow{}, Doublet{.strength = 1.0}, Vortex{.circulation = circulation}}}).field;
        std::vector<Vec2> points = findStagnationPoints(field, {{-3.0, -3.0}, {3.0, 3.0}}, 1.0);
        std::ranges::sort(points, {}, &Vec2::y);
        return points;
    };
    const std::vector<Vec2> apart = found(6.0 * pi);
    const std::vector<Vec2> merged = found(4.0 * pi);
    if (apart.size() != 2) return unable(std::format("Found {} stagnation points for Γ = 6π, expected 2.", apart.size()));
    if (merged.size() != 1) return unable(std::format("Found {} stagnation points for Γ = 4π, expected 1.", merged.size()));

    const double s = 1.5;
    const Vec2 outer{0.0, -(s + std::sqrt(s * s - 1.0))};
    const Vec2 inner{0.0, -(s - std::sqrt(s * s - 1.0))};
    const Vec2 bottom{0.0, -1.0};
    const auto shown = [](Vec2 p) { return std::format("({:.9f}, {:.9f})", p.x + 0.0, p.y + 0.0); };
    VerificationResult r;
    r.error = std::max({length(apart[0] - outer), length(apart[1] - inner), length(merged[0] - bottom)});
    r.samples = 3;
    r.tolerance = 1e-6;
    r.errorLabel = "max distance from the exact stagnation point, in cylinder radii";
    r.comparisons = {{"Γ = 6π: in the flow", shown(apart[0]), shown(outer)},
                     {"Γ = 6π: its mirror inside the cylinder", shown(apart[1]), shown(inner)},
                     {"Γ = 4π: double root", shown(merged[0]), shown(bottom)}};
    r.note = "Newton's method converges only linearly at a double root, which sets the tolerance.";
    return finish(std::move(r));
}

VerificationResult dividingLines() {
    struct Case {
        const char* name;
        Scene scene;
        Vec2 point;    // the stagnation point
        std::size_t lines;
        double circulation; // of the cylinder; the lines are ψ = 0
        bool cylinder;      // otherwise two sources on the x axis: the lines are the axes
    };
    const Case cases[] = {{"Cylinder", cylinderScene(), {1.0, 0.0}, 4, 0.0, true},
                          {"Lifting cylinder, Γ = 2π", liftingCylinderScene(), {std::cos(330.0 * deg), std::sin(330.0 * deg)}, 4, 2.0 * pi, true},
                          {"Double point, Γ = 4π", {.elements = {UniformFlow{}, Doublet{.strength = 1.0}, Vortex{.circulation = 4.0 * pi}}}, {0.0, -1.0}, 6, 4.0 * pi, true},
                          {"Between two equal sources", twoSourcesScene(), {0.0, 0.0}, 4, 0.0, false}};
    VerificationResult r;
    double spread = 0.0; // largest departure from evenly spaced directions at the point
    for (const Case& c : cases) {
        const FlowField field = modelOf(c.scene).field;
        const std::vector<Vec2> stagnation = findStagnationPoints(field, k_view, 1.0);
        const auto nearest = std::ranges::min_element(stagnation, {}, [&](Vec2 q) { return length(q - c.point); });
        if (nearest == stagnation.end() || length(*nearest - c.point) > 1e-5) return unable(std::format("No stagnation point found: {}.", c.name));
        const StreamlineSet set = dividingStreamlines(field, *nearest, stagnation, k_view);
        if (set.lines.size() != c.lines) {
            return unable(std::format("{}: {} dividing streamlines, expected {}.", c.name, set.lines.size(), c.lines));
        }
        double offset = 0.0;
        std::vector<double> directions;
        for (const auto& line : set.lines) {
            for (const Vec2 p : line) {
                // Near the centre ψ grows like 1/r, so its rounding says nothing about the line.
                if (c.cylinder && length(p) < 0.3) continue;
                offset = std::max(offset, c.cylinder ? std::abs(cylinderPsi(1.0, 1.0, c.circulation, p))
                                                     : std::min(std::abs(p.x), std::abs(p.y)));
                ++r.samples;
            }
            // The direction a line takes at the point, a short way along it, at each end that is at
            // the point: both, for a line that goes around a body and comes back.
            const auto direction = [&](auto first, auto last) {
                if (length(*first - c.point) > 1e-2) return;
                auto q = first;
                while (q + 1 != last && length(*q - c.point) < 0.03) ++q;
                directions.push_back(std::atan2(q->y - c.point.y, q->x - c.point.x));
            };
            direction(line.begin(), line.end());
            direction(line.rbegin(), line.rend());
        }
        // A line around a body is traced from both of its ends, so each direction may come twice.
        std::ranges::sort(directions);
        const auto duplicates = std::ranges::unique(directions, [](double a, double b) { return b - a < 5.0 * deg; });
        directions.erase(duplicates.begin(), duplicates.end());
        if (directions.size() != c.lines) {
            return unable(std::format("{}: {} directions at the point, expected {}.", c.name, directions.size(), c.lines));
        }
        directions.push_back(directions.front() + 2.0 * pi);
        double worst = 0.0;
        for (std::size_t k = 0; k + 1 < directions.size(); ++k) {
            worst = std::max(worst, std::abs(directions[k + 1] - directions[k] - 2.0 * pi / static_cast<double>(c.lines)));
        }
        spread = std::max(spread, worst);
        r.error = std::max(r.error, offset);
        r.comparisons.push_back({std::format("{}: lines through the point", c.name), std::to_string(c.lines), std::to_string(c.lines)});
        r.comparisons.push_back({std::format("{}: largest step between directions", c.name),
                                 std::format("{:.2f}°", (2.0 * pi / static_cast<double>(c.lines) + worst) / deg),
                                 std::format("{:.0f}°", 360.0 / static_cast<double>(c.lines))});
    }
    r.tolerance = 1e-5;
    r.errorLabel = "max |ψ| along the lines of the cylinders (distance from the axes for the two sources)";
    r.note = std::format("{} points of the lines drawn when a stagnation point is clicked. Directions are read 0.03 from "
                         "the point, where a line on a circle of radius 1 has turned by about 1°, and must be within 2° of even spacing.",
                         count(r.samples));
    if (spread > 2.0 * deg) r.error = inf;
    return finish(std::move(r));
}

// The one closed (or open) body found in `scene`, or nullopt.
std::optional<BodySurface> onlyBody(const FlowModel& model, bool closed) {
    const std::vector<BodySurface> bodies = findBodySurfaces(model, k_view);
    if (bodies.size() != 1 || bodies.front().closed != closed) return std::nullopt;
    return bodies.front();
}

VerificationResult rankineOval() {
    const FlowModel model = modelOf(rankineOvalScene());
    const std::optional<BodySurface> body = onlyBody(model, true);
    if (!body) return unable("Expected one closed body.");
    if (body->stagnationPoints.size() != 2) return unable("Expected two stagnation points on the oval.");

    // ψ = V∞ y + (Λ / 2π)(θ₁ − θ₂), with θ₁ and θ₂ the angles seen from the source and the sink.
    const auto psi = [](Vec2 p) { return p.y + (std::atan2(p.y, p.x + 1.0) - std::atan2(p.y, p.x - 1.0)); };
    double psiError = 0.0;
    double halfHeight = 0.0;
    for (const Vec2 p : body->points) {
        psiError = std::max(psiError, std::abs(psi(p)));
        halfHeight = std::max(halfHeight, std::abs(p.y));
    }
    // Stagnation points at x = ±√(b² + Λb / πV∞) = ±√3 (Anderson eq. 3.82); the half-height h solves
    // V∞ h = (Λ / 2π)(π − 2 atan(h / b)).
    double stagnation = 0.0;
    for (const Vec2 s : body->stagnationPoints) stagnation = std::max({stagnation, std::abs(std::abs(s.x) - std::sqrt(3.0)), std::abs(s.y)});
    const BodyForces forces = bodyForces(*body, model);

    VerificationResult r;
    r.error = std::max(psiError, stagnation);
    r.tolerance = 1e-6;
    r.errorLabel = "max |ψ| over the traced surface (should be 0), and the stagnation-point error, per unit b";
    r.samples = static_cast<int>(body->points.size());
    r.comparisons = {{"Stagnation points x", std::format("±{:.9f}", std::abs(body->stagnationPoints[0].x)), std::format("±√3 = ±{:.9f}", std::sqrt(3.0))},
                     {"Chord", std::format("{:.6f}", forces.chord), std::format("2√3 = {:.6f}", 2.0 * std::sqrt(3.0))},
                     {"Half-height h", std::format("{:.6f}", halfHeight), std::format("{:.6f}", ovalHalfHeight())},
                     {"max |ψ| on the surface", num(psiError), "0 (dividing streamline)"}};
    r.note = std::format("The surface is traced from the stagnation points with RK4 ({} points). The half-height "
                         "is the highest traced point, so it may fall short of h by a little.",
                         count(r.samples));
    return finish(std::move(r));
}

VerificationResult halfBody() {
    const FlowModel model = modelOf(halfBodyScene());
    const std::optional<BodySurface> body = onlyBody(model, false);
    if (!body) return unable("Expected one open body.");
    if (body->stagnationPoints.size() != 1) return unable("Expected one stagnation point.");

    // ψ = V∞ y + Λ θ / 2π about the source: ±Λ/2 on the upper and lower halves of the body.
    double psiError = 0.0;
    for (const Vec2 p : body->points) {
        const double psi = p.y + std::atan2(p.y, p.x + 1.0);
        psiError = std::max(psiError, std::abs(std::abs(psi) - pi));
    }
    const Vec2 s = body->stagnationPoints.front();
    const Vec2 end = body->points.back();
    VerificationResult r;
    r.error = std::max(psiError, length(s - Vec2{-2.0, 0.0}));
    r.tolerance = 1e-6;
    r.errorLabel = "max ||ψ| − Λ/2| over the traced surface, and the stagnation-point error, per unit Λ/2πV∞";
    r.samples = static_cast<int>(body->points.size());
    r.comparisons = {{"Stagnation point", std::format("({:.9f}, {:.9f})", s.x, s.y + 0.0), "(−2, 0): Λ/2πV∞ upstream of the source"},
                     {std::format("Half-width at x = {:.3g}", end.x), std::format("{:.6f}", std::abs(end.y)),
                      std::format("{:.6f}; → π far downstream", pi - std::atan2(end.y, end.x + 1.0))},
                     {"max ||ψ| − Λ/2| on the surface", num(psiError), "0 (dividing streamline)"}};
    r.note = std::format("{} traced surface points checked against ψ = ±Λ/2.", count(r.samples));
    return finish(std::move(r));
}

VerificationResult sourceNearWall() {
    const FlowModel model = modelOf(sourceNearWallScene());
    const auto exact = [](Vec2 p) { return exactSource(2.0 * pi, {0.0, 1.0}, p) + exactSource(2.0 * pi, {0.0, -1.0}, p); };
    const GridStats g = compareGrid(model.field, [&](Vec2 p) {
        return at(p, {0.0, 1.0}) || at(p, {0.0, -1.0}) ? k_singular : value(exact(p));
    });
    VerificationResult r = gridResult(g, {velocityAt(model.field, {1.0, 0.0}, exact({1.0, 0.0})),
                                          velocityAt(model.field, {0.0, 2.0}, exact({0.0, 2.0})),
                                          {"Image elements", std::to_string(model.images.size()), "1: a source Λ = 2π at (0, −1)"}});
    if (model.images.size() != 1) r.error = inf;
    return finish(std::move(r));
}

// --- Symmetry walls ---------------------------------------------------------------------

// Largest |v · n| / max |v| at `samples` points along the wall inside a 10 x 10 box.
double wallLeak(const FlowField& field, const SymmetryWall& wall, int samples, int& counted) {
    double leak = 0.0;
    double speed = 0.0;
    for (int k = 0; k < samples; ++k) {
        const Vec2 p = wall.point + wall.direction() * (-5.0 + 10.0 * (k + 0.5) / samples);
        const Vec2 v = field.velocity(p);
        leak = std::max(leak, std::abs(dot(v, wall.normal())));
        speed = std::max(speed, length(v));
        ++counted;
    }
    return speed > 0.0 ? leak / speed : inf;
}

VerificationResult inclinedWall() {
    const Scene scene = inclinedWallScene();
    const FlowModel model = modelOf(scene);
    const SymmetryWall wall = std::get<SymmetryWall>(scene.elements.back());
    VerificationResult r;
    const double leak = wallLeak(model.field, wall, 400, r.samples);
    const double along = std::cos(15.0 * deg); // the 10° stream keeps its component along the 25° wall
    r.error = std::max(leak, std::abs(model.freeStreamSpeed - along));
    r.tolerance = k_exactTolerance;
    r.errorLabel = "max |v·n| / max |v| along the wall";
    r.comparisons = {{"max |v·n| / max |v| on the wall", num(leak), "0"},
                     {"Image elements", std::to_string(model.images.size()), "3 (source, vortex, doublet)"},
                     {"Free stream left along the wall", std::format("{:.9f}", model.freeStreamSpeed), std::format("cos 15° = {:.9f}", along)}};
    r.note = "A wall at 25° through (0.3, −0.8), with a uniform flow, a source, a vortex and a doublet on one side. "
             "400 points along the wall.";
    if (model.images.size() != 3) r.error = inf;
    return finish(std::move(r));
}

VerificationResult cornerWalls() {
    const FlowModel model = modelOf(cornerWallsScene());
    // The images of a vortex in the corner of two perpendicular walls: -Γ in the neighbouring
    // quadrants, +Γ in the opposite one.
    const Vec2 p0{1.0, 2.0};
    const auto exact = [p0](Vec2 p) {
        return exactVortex(1.0, p0, p) + exactVortex(-1.0, {-p0.x, p0.y}, p) + exactVortex(-1.0, {p0.x, -p0.y}, p) +
               exactVortex(1.0, -p0, p);
    };
    const GridStats g = compareGrid(model.field, [&](Vec2 p) {
        return at(p, p0) || at(p, {-p0.x, p0.y}) || at(p, {p0.x, -p0.y}) || at(p, -p0) ? k_singular : value(exact(p));
    });
    int counted = 0;
    const double floor = wallLeak(model.field, SymmetryWall{}, 200, counted);
    const double side = wallLeak(model.field, SymmetryWall{.angle = 90.0 * deg}, 200, counted);
    VerificationResult r = gridResult(g, {{"Image elements", std::to_string(model.images.size()), "3 vortices: −Γ, −Γ, +Γ"},
                                          {"max |v·n| / max |v| on y = 0", num(floor), "0"},
                                          {"max |v·n| / max |v| on x = 0", num(side), "0"},
                                          velocityAt(model.field, {2.0, 0.0}, exact({2.0, 0.0}))});
    r.error = std::max({r.error, floor, side});
    if (model.images.size() != 3 || !model.imagesExact) r.error = inf;
    return finish(std::move(r));
}

// Velocity of one row of sources Λ at (x0, y0 + 2hn), n over all integers, summed by hand:
// u = Λ/4h sinh 2X / D, v = Λ/4h sin 2Y / D, X = π(x − x0)/2h, Y = π(y − y0)/2h, with
// D = cosh 2X − cos 2Y = 2 (sinh²X + sin²Y), which keeps its digits next to the source.
Vec2 exactSourceRow(double strength, double h, Vec2 centre, Vec2 p) {
    const double x = pi * (p.x - centre.x) / (2.0 * h);
    const double y = pi * (p.y - centre.y) / (2.0 * h);
    const double sx = std::sinh(x);
    const double sy = std::sin(y);
    const double scale = strength / (8.0 * h) / (sx * sx + sy * sy);
    return {scale * std::sinh(2.0 * x), scale * std::sin(2.0 * y)};
}

VerificationResult channelSource() {
    const Scene scene = channelSourceScene();
    const FlowModel model = modelOf(scene);
    const double h = 2.0;
    const double strength = 2.0 * pi;
    const Vec2 own{0.0, 0.4};
    const Vec2 image{0.0, -2.4}; // across y = −1
    const auto exact = [&](Vec2 p) {
        return exactSourceRow(strength, h, own, p) + exactSourceRow(strength, h, image, p);
    };
    const GridStats g = compareGrid(model.field, [&](Vec2 p) {
        return at(p, own) || at(p, image) || at(p, image + Vec2{0.0, 2.0 * h}) ? k_singular : value(exact(p));
    });
    int counted = 0;
    const double lower = wallLeak(model.field, std::get<SymmetryWall>(scene.elements[1]), 200, counted);
    const double upper = wallLeak(model.field, std::get<SymmetryWall>(scene.elements[2]), 200, counted);
    const double far = length(model.field.velocity({30.0, 0.0}));
    const double endSpeed = strength / (2.0 * h);
    VerificationResult r = gridResult(g, {{"Speed at (30, 0)", std::format("{:.9f}", far),
                                           std::format("Λ/2h = {:.9f}", endSpeed)},
                                          {"max |v·n| / max |v| on y = −1", num(lower), "0"},
                                          {"max |v·n| / max |v| on y = 1", num(upper), "0"},
                                          velocityAt(model.field, {1.0, 0.0}, exact({1.0, 0.0}))});
    r.error = std::max({r.error, lower, upper, std::abs(far - endSpeed) / endSpeed});
    if (!model.imagesExact || model.field.rows.size() != 2) r.error = inf;
    return finish(std::move(r));
}

VerificationResult inclinedChannel() {
    const Scene scene = inclinedChannelScene();
    const FlowModel model = modelOf(scene);
    const SymmetryWall lower = std::get<SymmetryWall>(scene.elements[4]);
    const SymmetryWall upper = std::get<SymmetryWall>(scene.elements[5]);
    const Vec2 t = lower.direction();
    const Vec2 n = lower.normal();
    const double h = 2.0;
    const double phi = 25.0 * deg;
    const double stream = std::cos(15.0 * deg);

    // Every copy in cell m of the row, written out with the polar forms: the element at z + 2hm n,
    // its mirror image across the lower wall at the same distance on the other side.
    const auto cell = [&](Vec2 p, int m) {
        const Vec2 shift = n * (2.0 * h * m);
        const auto across = [&](Vec2 q) { return q - n * (2.0 * dot(q - lower.point, n)); };
        const Vec2 source{1.0, 0.0};
        const Vec2 vortex{-1.5, -0.5};
        const Vec2 doublet{0.5, 0.9};
        const double beta = 70.0 * deg;
        return exactSource(4.0, source + shift, p) + exactSource(4.0, across(source) + shift, p) +
               exactVortex(-3.0, vortex + shift, p) + exactVortex(3.0, across(vortex) + shift, p) +
               exactDoublet(0.7, beta, doublet + shift, p) +
               exactDoublet(0.7, 2.0 * phi - beta, across(doublet) + shift, p);
    };
    // The partial sums over |m| <= N approach the limit as c / (N + 1/2), with the next term of order
    // 1 / N^3, so two of them give the limit far more closely than either.
    constexpr int shortSum = 1023;
    constexpr int longSum = 2047;
    const auto seriesSum = [&](Vec2 p) {
        Vec2 total = cell(p, 0);
        Vec2 first{};
        for (int m = 1; m <= longSum; ++m) {
            total += cell(p, m) + cell(p, -m);
            if (m == shortSum) first = total;
        }
        return (total * (longSum + 0.5) - first * (shortSum + 0.5)) / static_cast<double>(longSum - shortSum);
    };

    VerificationResult r;
    double worst = 0.0;
    double scale = 0.0;
    Vec2 sampleAt;
    Vec2 sample;
    for (int j = 0; j < 8; ++j) {
        for (int i = 0; i < 16; ++i) {
            const Vec2 p = lower.point + t * (-4.0 + 8.0 * (i + 0.5) / 16.0) + n * (h * (j + 0.5) / 8.0);
            const Vec2 sum = seriesSum(p) + t * stream;
            worst = std::max(worst, length(model.field.velocity(p) - sum));
            scale = std::max(scale, length(sum));
            sampleAt = p;
            sample = sum;
            ++r.samples;
        }
    }
    const double series = worst / scale;
    const double lowerLeak = wallLeak(model.field, lower, 200, r.samples);
    const double upperLeak = wallLeak(model.field, upper, 200, r.samples);
    r.error = std::max({series, lowerLeak, upperLeak, std::abs(model.freeStreamSpeed - stream)});
    r.tolerance = 1e-7; // what is left of the series after extrapolating, not the closed form
    r.errorLabel = "max |V − V_series| / max |V| in the channel, and |v·n| / max |v| on the walls";
    r.comparisons = {{"max |v·n| / max |v| on the lower wall", num(lowerLeak), "0"},
                     {"max |v·n| / max |v| on the upper wall", num(upperLeak), "0"},
                     {"max |V − V_series| / max |V|", num(series), "0 (to the series' own accuracy)"},
                     {"Rows of images", std::to_string(model.field.rows.size()), "6 (each element and its mirror)"},
                     velocityAt(model.field, sampleAt, sample)};
    r.note = std::format("{} copies of each element summed at {} points in a channel at 25°, 2 wide, with a uniform "
                         "flow, a source, a vortex and a doublet, and 400 points along the walls.",
                         count(2 * (2 * longSum + 1)), r.samples - 400);
    if (!model.imagesExact || model.field.rows.size() != 6) r.error = inf;
    return finish(std::move(r));
}

VerificationResult wedgeCorner() {
    // ζ = z³ opens the wedge 0 < θ < 60° onto the upper half plane, where one image at the conjugate
    // point closes the flow, so F(z) = G(z³) with no sum over images at all:
    // W = (Λ / 2π) 3z² [1 / (ζ − ζ₀) + 1 / (ζ − ζ̄₀)], and for a vortex iΓ/2π and a minus sign.
    const Complex z0 = toComplex(k_wedgeElement);
    const Complex zeta0 = z0 * z0 * z0;
    const auto exactW = [zeta0](Vec2 p, Complex strength, double image) {
        const Complex z = toComplex(p);
        const Complex zeta = z * z * z;
        return strength / (2.0 * pi) * 3.0 * z * z * (1.0 / (zeta - zeta0) + image / (zeta - std::conj(zeta0)));
    };
    const auto toVelocity = [](Complex w) { return Vec2{w.real(), -w.imag()}; };
    const auto exactSourcePair = [&](Vec2 p) { return toVelocity(exactW(p, 2.0 * pi, 1.0)); };
    const auto exactVortexPair = [&](Vec2 p) { return toVelocity(exactW(p, Complex{0.0, 2.0 * pi}, -1.0)); };
    // The element and its five images: its point turned through 120° and 240°, and their mirrors.
    std::vector<Vec2> singular;
    for (const double a : {0.0, 120.0 * deg, 240.0 * deg}) {
        singular.push_back(turned(k_wedgeElement, a));
        singular.push_back(mirroredInX(turned(k_wedgeElement, a)));
    }
    const auto onCopy = [&](Vec2 p) { return std::ranges::any_of(singular, [p](Vec2 q) { return at(p, q); }); };

    const Scene withSource = wedgeSourceScene();
    Scene withVortex = withSource;
    withVortex.elements[0] = Vortex{.circulation = 2.0 * pi, .position = k_wedgeElement};
    const FlowModel sourceModel = modelOf(withSource);
    const FlowModel vortexModel = modelOf(withVortex);
    GridStats g = compareGrid(sourceModel.field, [&](Vec2 p) { return onCopy(p) ? k_singular : value(exactSourcePair(p)); });
    g += compareGrid(vortexModel.field, [&](Vec2 p) { return onCopy(p) ? k_singular : value(exactVortexPair(p)); });
    const SymmetryWall floor = std::get<SymmetryWall>(withSource.elements[1]);
    const SymmetryWall side = std::get<SymmetryWall>(withSource.elements[2]);
    int counted = 0;
    const double leak = std::max({wallLeak(sourceModel.field, floor, 200, counted), wallLeak(sourceModel.field, side, 200, counted),
                                  wallLeak(vortexModel.field, floor, 200, counted), wallLeak(vortexModel.field, side, 200, counted)});
    VerificationResult r = gridResult(g,
                                      {{"Image elements", std::to_string(sourceModel.images.size()), "5"},
                                       {"max |v·n| / max |v| on both walls", num(leak), "0"},
                                       velocityAt(sourceModel.field, {1.0, 0.0}, exactSourcePair({1.0, 0.0})),
                                       velocityAt(vortexModel.field, {1.0, 0.0}, exactVortexPair({1.0, 0.0}))},
                                      2);
    r.error = std::max(r.error, leak);
    if (sourceModel.images.size() != 5 || vortexModel.images.size() != 5 || !sourceModel.imagesExact || !vortexModel.imagesExact) {
        r.error = inf;
    }
    return finish(std::move(r));
}

VerificationResult channelVortex() {
    const Scene scene = channelVortexScene();
    const FlowModel model = modelOf(scene);
    const double h = 2.0;
    const double gamma = 2.0 * pi;
    // Images of alternating sign every h along y sum to W = iΓ / (2h sinh(πz/h)). In real form, with
    // a = πx/h, b = πy/h: u = Γ/2h cosh a sin b / D, v = −Γ/2h sinh a cos b / D, D = sinh²a + sin²b.
    const auto exact = [=](Vec2 p) {
        const double a = pi * p.x / h;
        const double b = pi * p.y / h;
        const double sa = std::sinh(a);
        const double sb = std::sin(b);
        const double d = sa * sa + sb * sb;
        return Vec2{gamma / (2.0 * h) * std::cosh(a) * sb / d, -gamma / (2.0 * h) * sa * std::cos(b) / d};
    };
    const GridStats g = compareGrid(model.field, [&](Vec2 p) {
        const bool copy = at(p, {}) || at(p, {0.0, 2.0}) || at(p, {0.0, -2.0});
        return copy ? k_singular : value(exact(p));
    });
    int counted = 0;
    const double lower = wallLeak(model.field, std::get<SymmetryWall>(scene.elements[1]), 200, counted);
    const double upper = wallLeak(model.field, std::get<SymmetryWall>(scene.elements[2]), 200, counted);
    // On the wall the speed falls off as a hyperbolic secant, far faster than a lone vortex's 1/r.
    const auto wallSpeed = [=](double x) { return gamma / (2.0 * h) / std::cosh(pi * x / h); };
    double slip = 0.0;
    for (int k = 0; k <= 40; ++k) {
        const double x = -4.0 + 0.2 * k;
        slip = std::max(slip, std::abs(model.field.velocity({x, 1.0}).x - wallSpeed(x)) / wallSpeed(0.0));
        ++counted;
    }
    VerificationResult r = gridResult(g, {{"Speed on the wall at x = 0", std::format("{:.9f}", length(model.field.velocity({0.0, 1.0}))), std::format("Γ/2h = {:.9f}", wallSpeed(0.0))},
                                          {"Speed on the wall at x = 2", std::format("{:.9f}", length(model.field.velocity({2.0, 1.0}))), std::format("Γ/2h sech π = {:.9f}", wallSpeed(2.0))},
                                          {"max |v·n| / max |v| on y = ±1", num(std::max(lower, upper)), "0"}});
    r.error = std::max({r.error, lower, upper, slip});
    if (!model.imagesExact || model.field.rows.size() != 2) r.error = inf;
    return finish(std::move(r));
}

// --- Pressure and forces ----------------------------------------------------------------

VerificationResult surfacePressure() {
    VerificationResult r;
    for (const double gamma : {0.0, 2.0 * pi}) {
        const FlowModel model = modelOf({.elements = {UniformFlow{}, Doublet{.strength = 1.0}, Vortex{.circulation = gamma}}});
        const std::optional<BodySurface> body = onlyBody(model, true);
        if (!body) return unable("Expected one closed body for each cylinder.");
        // Cp = 1 − (2 sin θ + Γ / 2πaV∞)² on the cylinder r = a = 1.
        double cpError = 0.0;
        double radius = 0.0;
        for (std::size_t k = 0; k < body->points.size(); ++k) {
            const Vec2 p = body->points[k];
            const double q = 2.0 * std::sin(std::atan2(p.y, p.x)) + gamma / (2.0 * pi);
            cpError = std::max(cpError, std::abs(body->cp[k] - (1.0 - q * q)));
            radius = std::max(radius, std::abs(length(p) - 1.0));
        }
        r.error = std::max({r.error, cpError, radius});
        r.samples += static_cast<int>(body->points.size());
        const char* name = gamma == 0.0 ? "Γ = 0" : "Γ = 2π";
        r.comparisons.push_back({std::format("{}: lowest Cp", name), std::format("{:.6f}", *std::ranges::min_element(body->cp)),
                                 gamma == 0.0 ? "−3 (|V| = 2V∞ at the top)" : "−8 (|V| = 3V∞ at the top)"});
        r.comparisons.push_back({std::format("{}: highest Cp", name), std::format("{:.6f}", *std::ranges::max_element(body->cp)),
                                 "1 (stagnation)"});
        r.comparisons.push_back({std::format("{}: max |Cp − Cp_exact|", name), num(cpError), "0"});
    }
    r.tolerance = 1e-6;
    r.errorLabel = "max |Cp − Cp_exact| and max |r − a| / a over the traced surfaces";
    r.note = std::format("{} surface points of the traced dividing streamlines, with and without circulation.", count(r.samples));
    return finish(std::move(r));
}

VerificationResult rankinePressure() {
    VerificationResult r;
    // The oval: Bernoulli with the velocity written out; the top of it is the fastest point, where
    // the source and the sink each add Λ/2π · b / (b² + h²) to V∞.
    const FlowModel ovalModel = modelOf(rankineOvalScene());
    const std::optional<BodySurface> oval = onlyBody(ovalModel, true);
    if (!oval) return unable("Expected one closed body.");
    for (std::size_t k = 0; k < oval->points.size(); ++k) {
        const Vec2 p = oval->points[k];
        const Vec2 v = Vec2{1.0, 0.0} + exactSource(2.0 * pi, {-1.0, 0.0}, p) + exactSource(-2.0 * pi, {1.0, 0.0}, p);
        r.error = std::max(r.error, std::abs(oval->cp[k] - (1.0 - lengthSquared(v))));
        ++r.samples;
    }
    const double h = ovalHalfHeight();
    const double fastest = 1.0 + 2.0 / (1.0 + h * h);
    const double ovalLowest = *std::ranges::min_element(oval->cp);
    const double ovalExact = 1.0 - fastest * fastest;
    r.error = std::max(r.error, std::abs(ovalLowest - ovalExact));

    // The half-body: on its surface r = b (π − θ) / sin θ about the source (b = Λ/2πV∞ = 1), so Cp
    // depends on θ alone: Cp = −[sin²θ / (π − θ)² + 2 sin θ cos θ / (π − θ)].
    const auto cpAt = [](double theta) {
        const double q = std::sin(theta) / (pi - theta);
        return -(q * q + 2.0 * q * std::cos(theta));
    };
    const FlowModel halfModel = modelOf(halfBodyScene());
    const std::optional<BodySurface> half = onlyBody(halfModel, false);
    if (!half) return unable("Expected one open body.");
    for (std::size_t k = 0; k < half->points.size(); ++k) {
        const double theta = std::abs(std::atan2(half->points[k].y, half->points[k].x + 1.0));
        if (pi - theta < 1e-3) continue; // 0 / 0 at the stagnation point itself
        r.error = std::max(r.error, std::abs(half->cp[k] - cpAt(theta)));
        ++r.samples;
    }
    double halfExact = 1.0;
    for (int k = 1; k < 20000; ++k) halfExact = std::min(halfExact, cpAt(pi * k / 20000.0));
    const double halfLowest = *std::ranges::min_element(half->cp);
    r.error = std::max(r.error, std::abs(halfLowest - halfExact));

    r.tolerance = 1e-5;
    r.errorLabel = "max |Cp − Cp_exact| over the traced surfaces, and the error of the lowest Cp";
    r.comparisons = {{"Oval: lowest Cp, at the top", std::format("{:.6f}", ovalLowest), std::format("1 − (1 + 2b²/(b² + h²))² = {:.6f}", ovalExact)},
                     {"Half-body: lowest Cp", std::format("{:.6f}", halfLowest), std::format("{:.6f}", halfExact)},
                     {"Half-body: fastest speed", std::format("{:.4f} V∞", std::sqrt(1.0 - halfLowest)), "1.26 V∞ (textbook)"}};
    r.note = std::format("{} surface points. The half-body's formula depends only on the angle seen from the source, so "
                         "it tests where the surface was traced as well as the pressure.",
                         count(r.samples));
    return finish(std::move(r));
}

// F_x − i F_y = (iρ/2) ∮ W² dz counter-clockwise around a circle, with ρ = 1 (trapezoidal rule,
// which converges geometrically on a circle clear of every singular point).
Complex blasiusForce(const FlowField& field, Vec2 centre, double radius) {
    constexpr int n = 4096;
    Complex sum{};
    for (int k = 0; k < n; ++k) {
        const Complex offset = std::polar(radius, 2.0 * pi * k / n);
        const Complex w = field.complexVelocity(centre + Vec2{offset.real(), offset.imag()});
        sum += w * w * Complex{0.0, 1.0} * offset;
    }
    return Complex{0.0, 0.5} * sum * (2.0 * pi / n);
}

VerificationResult blasius() {
    struct Case {
        const char* name;
        Scene scene;
        Vec2 centre;
        double radius;
        double speed;
        double alpha;
        double strength;    // total source strength inside the circle
        double circulation; // total clockwise circulation inside the circle
    };
    const Case cases[] = {{"Lifting cylinder at 20°", liftScene(), {0.5, -0.3}, 2.0, 2.0, 20.0 * deg, 0.0, 2.0 * pi},
                          {"Rankine oval", rankineOvalScene(), {}, 3.0, 1.0, 0.0, 0.0, 0.0},
                          {"Half-body", halfBodyScene(), {-1.0, 0.0}, 1.5, 1.0, 0.0, 2.0 * pi, 0.0},
                          {"Every element", mixedScene(), {}, 3.0, 1.0, 15.0 * deg, 0.7, 3.0}};
    VerificationResult r;
    for (const Case& c : cases) {
        const FlowField field = modelOf(c.scene).field;
        const Complex program = blasiusForce(field, c.centre, c.radius);
        // Only the 1/z term of W² has a residue: F_x − i F_y = −ρ W∞ (ΣΛ + iΣΓ).
        const Complex exact = -std::polar(c.speed, -c.alpha) * Complex{c.strength, c.circulation};
        r.error = std::max(r.error, std::abs(program - exact) / std::max(std::abs(exact), 1.0));
        r.samples += 4096;
        const auto force = [](Complex f) {
            const auto zero = [](double x) { return std::abs(x) < 1e-11 ? 0.0 : x; };
            return std::format("({}, {})", num(zero(f.real())), num(zero(-f.imag())));
        };
        r.comparisons.push_back({std::format("{}: (F_x, F_y) / ρ", c.name), force(program), force(exact)});
    }
    r.tolerance = 1e-10;
    r.errorLabel = "max |F − F_exact| / max(|F_exact|, 1), F from ∮ W² dz on a circle around the elements";
    r.note = "No surface is traced: the force comes from the flow on a circle well away from the body. A net source feels "
             "a force −ρV∞Λ along the stream (Lagally), and circulation the Kutta–Joukowski lift.";
    return finish(std::move(r));
}

VerificationResult sweepForces() {
    const Scene scene = liftingCylinderScene();
    const ParameterRef circulation{.element = 2, .parameter = 0};
    const std::vector<double> values = sweepValues(-3.0 * pi, 3.0 * pi, 7);
    VerificationResult r;
    double liftError = 0.0;
    double dragError = 0.0;
    double stagnationError = 0.0;
    for (const double gamma : values) {
        const std::optional<SweepSample> sample = sampleSweep(scene, circulation, gamma, k_view);
        if (!sample) return unable("The sweep refused a value of Γ.");
        const auto closed = std::ranges::find_if(sample->bodies, &SweepBody::closed);
        if (closed == sample->bodies.end()) return unable(std::format("No closed body at Γ = {}.", num(gamma)));
        if (sample->stagnationPoints.size() != 2) return unable(std::format("{} stagnation points at Γ = {}, expected 2.", sample->stagnationPoints.size(), num(gamma)));
        const double exactLift = 2.0 * gamma; // 2Γ/V∞
        const double scale = std::max(std::abs(exactLift), 1.0);
        liftError = std::max({liftError, std::abs(closed->lift - exactLift) / scale, std::abs(closed->exactLift - exactLift) / scale});
        dragError = std::max(dragError, std::abs(closed->drag) / closed->chord);
        for (const Vec2 p : sample->stagnationPoints) {
            // On the circle, at sin θ = −Γ / 4πaV∞.
            stagnationError = std::max({stagnationError, std::abs(length(p) - 1.0), std::abs(p.y + gamma / (4.0 * pi))});
        }
        r.samples += 3;
        r.comparisons.push_back({std::format("Γ = {:.4f}: L′/q∞", gamma), std::format("{:.6f}", closed->lift), std::format("{:.6f}", exactLift)});
    }
    r.error = std::max({liftError, dragError, stagnationError});
    r.tolerance = 1e-4;
    r.errorLabel = "max lift and drag errors, and the stagnation-point error, over the sweep";
    r.comparisons.push_back({"Largest |D′/q∞| / c", num(dragError), "0"});
    r.comparisons.push_back({"Largest stagnation-point error", num(stagnationError), "0"});
    r.note = "The numbers the sweep exports, for Γ from −3π to 3π on a cylinder a = 1 in V∞ = 1. Lift is linear in Γ, "
             "and the stagnation points slide around the circle.";
    return finish(std::move(r));
}

VerificationResult dAlembert() {
    struct Case {
        const char* name;
        Scene scene;
    };
    const Case cases[] = {{"Cylinder", cylinderScene()}, {"Lifting cylinder, Γ = 2π", liftingCylinderScene()},
                          {"Rankine oval", rankineOvalScene()}};
    VerificationResult r;
    for (const Case& c : cases) {
        const FlowModel model = modelOf(c.scene);
        const std::optional<BodySurface> body = onlyBody(model, true);
        if (!body) return unable(std::format("Expected one closed body: {}.", c.name));
        const BodyForces forces = bodyForces(*body, model);
        r.error = std::max(r.error, std::abs(forces.drag) / forces.chord);
        r.samples += static_cast<int>(body->points.size());
        r.comparisons.push_back({std::format("{}: D′/q∞", c.name), num(forces.drag), "0"});
    }
    r.tolerance = 1e-5;
    r.errorLabel = "max |D′/q∞| / c, the drag from integrating −Cp n around the surface";
    r.note = "The trapezoidal rule over the traced surface; its error, not the flow, sets how close to 0 the drag comes.";
    return finish(std::move(r));
}

VerificationResult kuttaJoukowski() {
    struct Case {
        const char* name;
        Scene scene;
        double gamma;
        double speed;
    };
    const Case cases[] = {
        {"V∞ = 2 at 20°, Γ = 2π", liftScene(), 2.0 * pi, 2.0},
        {"V∞ = 1.5, Γ = −π (downforce)",
         {.elements = {UniformFlow{.speed = 1.5}, Doublet{.strength = 1.5}, Vortex{.circulation = -pi}}},
         -pi,
         1.5},
        {"V∞ = 1, Γ = 6π (stagnation point off the body)", {.elements = {UniformFlow{}, Doublet{.strength = 1.0}, Vortex{.circulation = 6.0 * pi}}}, 6.0 * pi, 1.0},
    };
    VerificationResult r;
    for (const Case& c : cases) {
        const FlowModel model = modelOf(c.scene);
        const std::vector<BodySurface> bodies = findBodySurfaces(model, k_view);
        const auto closed = std::ranges::find_if(bodies, &BodySurface::closed);
        if (closed == bodies.end()) return unable(std::format("No closed body found: {}.", c.name));
        const BodyForces forces = bodyForces(*closed, model);
        const double exact = 2.0 * c.gamma / c.speed; // L′ = ρ V∞ Γ, so L′ / q∞ = 2Γ / V∞
        r.error = std::max(r.error, std::abs(forces.lift - exact) / std::abs(exact));
        r.samples += static_cast<int>(closed->points.size());
        r.comparisons.push_back({std::format("{}: L′/q∞", c.name), std::format("{:.6f}", forces.lift), std::format("2Γ/V∞ = {:.6f}", exact)});
    }
    r.tolerance = 1e-4;
    r.errorLabel = "max |L′ − L′_exact| / |L′_exact|, the lift from integrating −Cp n around the surface";
    r.note = "With Γ > 4πaV∞ the body found is the closed dividing streamline around the cylinder, which carries "
             "the same lift. The trapezoidal rule over the traced surface sets the tolerance.";
    return finish(std::move(r));
}

// --- Streamlines and equipotentials -----------------------------------------------------

// Lines over the default view spaced as in the plot, with arrows.
StreamlineSettings plotSettings(double arrowSpacing = 0.7) {
    StreamlineSettings settings;
    settings.bounds = k_view.expanded(0.05);
    settings.separation = k_view.size().y / PlotSettings{}.streamlineDensity;
    settings.singularRadius = 0.02;
    settings.arrowSpacing = arrowSpacing;
    return settings;
}

VerificationResult uniformLines() {
    VerificationResult r;
    double crossing = 0.0; // |sin| of the angle between a step and the stream
    for (const double alpha : {0.0, 45.0 * deg}) {
        const FlowField field = modelOf({.elements = {UniformFlow{.speed = 1.0, .angle = alpha}}}).field;
        const Vec2 stream = exactUniform(1.0, alpha);
        const StreamlineSet set = computeStreamlines(field, plotSettings());
        bool downstream = !set.lines.empty();
        for (const auto& line : set.lines) {
            for (std::size_t k = 0; k + 1 < line.size(); ++k) {
                const Vec2 d = line[k + 1] - line[k];
                crossing = std::max(crossing, std::abs(cross(d, stream)) / length(d));
                downstream = downstream && dot(d, stream) > 0.0;
                ++r.samples;
            }
        }
        double arrows = 0.0;
        for (const Arrow& a : set.arrows) arrows = std::max(arrows, length(a.direction - stream));
        r.error = std::max(r.error, arrows);
        if (!downstream || set.arrows.empty()) r.error = inf;
        r.comparisons.push_back({std::format("α = {:.0f}°: {} lines, steps run", alpha / deg, set.lines.size()),
                                 downstream ? "downstream" : "upstream", "downstream"});
        r.comparisons.push_back({std::format("α = {:.0f}°: largest arrow error", alpha / deg), num(arrows), "0"});
    }
    r.error = std::max(r.error, crossing);
    r.comparisons.push_back({"Largest angle between a line and the stream", std::format("{}°", num(std::asin(std::min(crossing, 1.0)) / deg)), "0°"});
    r.tolerance = k_exactTolerance;
    r.errorLabel = "max |sin| of the angle between any step and the stream, and max arrow error";
    r.note = std::format("{} steps of the lines the plot draws, at α = 0° and 45°.", count(r.samples));
    return finish(std::move(r));
}

VerificationResult radialLines() {
    const Vec2 centre{0.3, -0.2};
    VerificationResult r;
    double offRay = 0.0; // largest distance from the ray through the line's far end
    for (const double strength : {2.0 * pi, -2.0 * pi}) {
        const FlowField field = modelOf({.elements = {SourceSink{.strength = strength, .position = centre}}}).field;
        const double outward = strength > 0.0 ? 1.0 : -1.0;
        const StreamlineSet set = computeStreamlines(field, plotSettings());
        bool direction = set.lines.size() >= 8;
        for (const auto& line : set.lines) {
            const Vec2 far = length(line.front() - centre) > length(line.back() - centre) ? line.front() : line.back();
            const Vec2 ray = (far - centre) / length(far - centre);
            for (std::size_t k = 0; k < line.size(); ++k) {
                offRay = std::max(offRay, std::abs(cross(line[k] - centre, ray)));
                if (k > 0) direction = direction && outward * (length(line[k] - centre) - length(line[k - 1] - centre)) > 0.0;
                ++r.samples;
            }
        }
        for (const Arrow& a : set.arrows) {
            const Vec2 radial = (a.position - centre) / length(a.position - centre);
            r.error = std::max(r.error, length(a.direction - radial * outward));
        }
        if (!direction || set.arrows.empty()) r.error = inf;
        r.comparisons.push_back({std::format("{}: {} lines run", strength > 0.0 ? "Source" : "Sink", set.lines.size()),
                                 direction ? (strength > 0.0 ? "outwards" : "inwards") : "the wrong way",
                                 strength > 0.0 ? "outwards" : "inwards"});
    }
    r.error = std::max(r.error, offRay);
    r.comparisons.push_back({"Largest distance from a straight ray", num(offRay), "0"});
    r.tolerance = 1e-9;
    r.errorLabel = "max distance of any point from its ray, and max arrow error";
    r.note = std::format("{} points of the lines the plot draws, for Λ = ±2π. RK4 keeps a straight line straight, so only "
                         "rounding and the direction at the seed are left.",
                         count(r.samples));
    return finish(std::move(r));
}

VerificationResult vortexCircles() {
    const FlowField field = modelOf(vortexScene()).field;
    const FlowField reversed = modelOf({.elements = {Vortex{.circulation = -2.0 * pi}}}).field;
    StreamlineSettings settings;
    settings.bounds = {{-2.0, -2.0}, {2.0, 2.0}};
    settings.separation = 0.3;
    const StreamlineSet set = computeStreamlines(field, settings);
    const StreamlineSet reversedSet = computeStreamlines(reversed, settings);
    const auto turning = [](const StreamlineSet& lines, double sign) {
        return !lines.arrows.empty() && std::ranges::all_of(lines.arrows, [sign](const Arrow& a) {
            return sign * cross(a.position, a.direction) > 0.0;
        });
    };
    const bool clockwise = turning(set, -1.0);
    const bool anticlockwise = turning(reversedSet, 1.0);
    VerificationResult r;
    int closed = 0;
    for (const auto& line : set.lines) {
        const double radius = length(line.front());
        for (const Vec2 p : line) r.error = std::max(r.error, std::abs(length(p) - radius) / radius);
        r.samples += static_cast<int>(line.size());
        if (line.front() == line.back() && line.size() > 10) ++closed;
    }
    r.tolerance = 1e-5;
    r.errorLabel = "max |r − r₀| / r₀ along each traced streamline";
    r.comparisons = {{"Streamlines traced", std::to_string(set.lines.size()), "—"},
                     {"Closed on themselves", std::to_string(closed), "every complete circle"},
                     {"Largest radius drift", num(r.error), "0 (circles r = const)"},
                     {"Γ = 2π: every arrow turns", clockwise ? "clockwise" : "the wrong way", "clockwise"},
                     {"Γ = −2π: every arrow turns", anticlockwise ? "counter-clockwise" : "the wrong way", "counter-clockwise"}};
    r.note = "RK4 in arc-length steps; its truncation error over one loop sets the tolerance.";
    if (closed < 3 || !clockwise || !anticlockwise) r.error = inf;
    return finish(std::move(r));
}

// The largest change of `psi` along any line that stays `clearance` away from the origin.
double psiDrift(const StreamlineSet& set, const std::function<double(Vec2)>& psi, double clearance, int& samples, int& lines) {
    double drift = 0.0;
    for (const auto& line : set.lines) {
        if (!std::ranges::all_of(line, [clearance](Vec2 p) { return length(p) > clearance; })) continue;
        ++lines;
        const double start = psi(line.front());
        for (const Vec2 p : line) drift = std::max(drift, std::abs(psi(p) - start));
        samples += static_cast<int>(line.size());
    }
    return drift;
}

VerificationResult streamFunction() {
    const double gamma = 2.0 * pi;
    const StreamlineSettings settings = plotSettings(0.0);
    const auto cylinderStream = [gamma](Vec2 p) { return cylinderPsi(1.0, 1.0, gamma, p); };
    const auto field = [](double circulation) {
        return modelOf({.elements = {UniformFlow{}, Doublet{.strength = 1.0}, Vortex{.circulation = circulation}}}).field;
    };
    VerificationResult r;
    int lines = 0;
    const double cylinderDrift = psiDrift(computeStreamlines(field(gamma), settings), cylinderStream, 0.5, r.samples, lines);

    // ψ = V r² sin 2θ = 2V xy: the streamlines of a corner are hyperbolas xy = const.
    StreamlineSettings cornerSettings = settings;
    cornerSettings.bounds = {{0.0, 0.0}, {4.0, 2.5}};
    const double cornerDrift = psiDrift(computeStreamlines(modelOf(cornerScene()).field, cornerSettings),
                                        [](Vec2 p) { return 2.0 * p.x * p.y; }, 0.0, r.samples, lines);

    // Lines kept from frame to frame during a sweep are still streamlines of the final flow.
    StreamlineSet kept = computeStreamlines(field(0.9 * gamma), settings);
    for (int k = 1; k <= 20; ++k) kept = computeStreamlines(field((0.9 + 0.005 * k) * gamma), settings, kept.seeds);
    const double keptDrift = psiDrift(kept, cylinderStream, 0.5, r.samples, lines);

    r.error = std::max({cylinderDrift, cornerDrift, keptDrift}) / settings.separation;
    r.tolerance = 1e-3;
    r.errorLabel = "max change of ψ along a line, in units of the spacing Δψ ≈ V∞ · separation";
    r.comparisons = {{"Lifting cylinder: largest change of ψ", num(cylinderDrift), "0"},
                     {"Corner: largest change of xy along a line", num(0.5 * cornerDrift), "0 (hyperbolas)"},
                     {"Lines kept through a sweep: largest change of ψ", num(keptDrift), "0"},
                     {"Lines checked", std::to_string(lines), "—"}};
    r.note = "ψ is written out in polar form, independently of the velocity that the lines are traced along. Lines "
             "passing within 0.5 of the cylinder's centre are left out.";
    if (lines < 20) r.error = inf;
    return finish(std::move(r));
}

VerificationResult evenSpacing() {
    VerificationResult r;
    // A uniform flow: neighbouring lines exactly one separation apart.
    const Vec2 normal{-std::sin(30.0 * deg), std::cos(30.0 * deg)};
    const StreamlineSettings settings = plotSettings();
    const double separation = settings.separation;
    const StreamlineSet parallel =
        computeStreamlines(modelOf({.elements = {UniformFlow{.speed = 1.0, .angle = 30.0 * deg}}}).field, settings);
    std::vector<double> offsets;
    for (const auto& line : parallel.lines) offsets.push_back(dot(line.front(), normal));
    std::ranges::sort(offsets);
    double spacing = 0.0;
    for (std::size_t k = 0; k + 1 < offsets.size(); ++k) {
        spacing = std::max(spacing, std::abs(offsets[k + 1] - offsets[k] - separation) / separation);
        ++r.samples;
    }

    // Around a cylinder: no two lines closer than half a separation (where a line stops), and no
    // point of the plot further from a line than a seed would need to fit one.
    const StreamlineSet curved = computeStreamlines(modelOf(liftingCylinderScene()).field, settings);
    double closest = inf;
    for (std::size_t a = 0; a < curved.lines.size(); ++a) {
        for (std::size_t b = a + 1; b < curved.lines.size(); ++b) {
            for (const Vec2 p : curved.lines[a]) {
                for (const Vec2 q : curved.lines[b]) closest = std::min(closest, lengthSquared(p - q));
            }
        }
    }
    closest = std::sqrt(closest);
    double widest = 0.0;
    for (int j = 0; j < 24; ++j) {
        for (int i = 0; i < 40; ++i) {
            const Vec2 p{k_view.min.x + k_view.size().x * (i + 0.5) / 40.0, k_view.min.y + k_view.size().y * (j + 0.5) / 24.0};
            double nearest = inf;
            for (const auto& line : curved.lines) {
                for (const Vec2 q : line) nearest = std::min(nearest, lengthSquared(p - q));
            }
            widest = std::max(widest, std::sqrt(nearest));
            ++r.samples;
        }
    }
    r.error = spacing;
    r.tolerance = 1e-9;
    r.errorLabel = "max |gap − d| / d in a uniform flow (lines closer than d/2 or further than d apart around a cylinder fail outright)";
    r.comparisons = {{"Uniform flow: spacing", std::format("{:.6f}", separation * (1.0 + spacing)), std::format("d = {:.6f}", separation)},
                     {"Cylinder: closest two lines", std::format("{:.3f} d", closest / separation), "at least 0.5 d"},
                     {"Cylinder: furthest point from a line", std::format("{:.3f} d", widest / separation), "at most 1 d"}};
    r.note = std::format("{} lines in the uniform flow and {} around the lifting cylinder, seeded as the plot does "
                         "(Jobard and Lefer), with d the streamline separation.",
                         parallel.lines.size(), curved.lines.size());
    if (closest < 0.5 * separation || widest > separation) r.error = inf;
    return finish(std::move(r));
}

VerificationResult arrows() {
    // On a curved line: along the flow where the arrow stands.
    const FlowField field = modelOf(liftingCylinderScene()).field;
    const StreamlineSettings settings = plotSettings();
    const StreamlineSet curved = computeStreamlines(field, settings);
    double angle = 0.0;
    VerificationResult r;
    for (const Arrow& a : curved.arrows) {
        if (length(a.position) < 0.5) continue;
        const Vec2 v = field.velocity(a.position);
        const double along = dot(a.direction, v) / length(v);
        if (along <= 0.0) r.error = inf;
        angle = std::max(angle, std::abs(cross(a.direction, v)) / length(v));
        ++r.samples;
    }
    // On straight lines: one every arrowSpacing along each.
    const StreamlineSet straight = computeStreamlines(modelOf({.elements = {UniformFlow{}}}).field, settings);
    double spacing = 0.0;
    std::vector<Arrow> sorted = straight.arrows;
    std::ranges::sort(sorted, [](const Arrow& a, const Arrow& b) {
        return a.position.y != b.position.y ? a.position.y < b.position.y : a.position.x < b.position.x;
    });
    for (std::size_t k = 0; k + 1 < sorted.size(); ++k) {
        if (std::abs(sorted[k + 1].position.y - sorted[k].position.y) > 1e-9) continue;
        spacing = std::max(spacing, std::abs(sorted[k + 1].position.x - sorted[k].position.x - settings.arrowSpacing) / settings.arrowSpacing);
        ++r.samples;
    }
    r.error = std::max({r.error, angle, spacing});
    r.tolerance = 0.05;
    r.errorLabel = "max |sin| of the angle between an arrow and the velocity at its foot, and the arrow spacing error";
    r.comparisons = {{"Largest angle from the flow", std::format("{:.3f}°", std::asin(std::min(angle, 1.0)) / deg), "0°, to the chord of a step"},
                     {"Arrows pointing upstream", std::isinf(r.error) ? "some" : "none", "none"},
                     {"Largest spacing error on straight lines", num(spacing), "0"}};
    r.note = "An arrow lies along a straight step of the traced line, so on a curved line it can differ from the local "
             "flow by half the angle the line turns in a step. The tolerance is that angle for a step of d/6 on a line "
             "of radius 0.35, about 3°. Arrows within 0.5 of the centre are left out.";
    return finish(std::move(r));
}

VerificationResult equipotentials() {
    const Scene scene = equipotentialScene();
    const FlowField field = modelOf(scene).field;
    StreamlineSettings settings;
    settings.bounds = {{-3.0, -3.0}, {3.0, 3.0}};
    settings.separation = 0.25;
    settings.arrowSpacing = 0.0;
    settings.orthogonal = true;
    const StreamlineSet set = computeStreamlines(field, settings);

    // Away from the stagnation points on the cylinder, where the direction field turns fastest.
    double angle = 0.0; // |cos| of the angle between a step and the velocity
    double drift = 0.0; // of φ along a line, in units of the spacing Δφ ≈ V∞ · separation
    int lines = 0;
    VerificationResult r;
    for (const auto& line : set.lines) {
        if (!std::ranges::all_of(line, [](Vec2 p) { return length(p) > 1.2; })) continue;
        ++lines;
        const double phi0 = field.complexPotential(line.front()).real();
        for (std::size_t i = 0; i + 1 < line.size(); ++i) {
            const Vec2 d = line[i + 1] - line[i];
            const Vec2 v = field.velocity((line[i] + line[i + 1]) * 0.5);
            if (length(d) <= 0.0) continue;
            angle = std::max(angle, std::abs(dot(d / length(d), v / length(v))));
            drift = std::max(drift, std::abs(field.complexPotential(line[i + 1]).real() - phi0) / settings.separation);
            ++r.samples;
        }
    }
    r.error = std::max(angle, drift);
    r.tolerance = 1e-3;
    r.errorLabel = "max |cos| of the crossing angle, and max φ drift along a line in units of Δφ";
    r.comparisons = {{"Equipotential lines checked", std::to_string(lines), "—"},
                     {"Largest departure from 90°", std::format("{:.4f}°", std::asin(std::min(angle, 1.0)) / deg), "0°"},
                     {"Largest change of φ along a line", std::format("{} Δφ", num(drift)), "0"}};
    r.note = "Lines of constant φ around a cylinder at 20°, traced as curves normal to the velocity.";
    if (lines < 10) r.error = inf;
    return finish(std::move(r));
}

// --- Particles and dye ------------------------------------------------------------------

constexpr double k_plotStep = 8.0 / 300.0;   // the plot's largest particle step in the default view
constexpr double k_frame = 1.0 / 60.0;

VerificationResult particles() {
    VerificationResult r;
    // A vortex Γ = 2π: circles at the angular speed Γ / 2πr² = 1 / r², clockwise.
    const FlowField vortexField = modelOf(vortexScene()).field;
    ParticleSettings settings;
    settings.bounds = {{-3.0, -3.0}, {3.0, 3.0}};
    settings.maxAge = 1e9;
    settings.maxStepLength = k_plotStep;
    ParticleSystem orbiting(1000, settings);
    const auto inRing = [](Vec2 p) { return length(p) > 0.5 && length(p) < 2.9; };
    const auto vortexVelocity = [&](Vec2 p) { return vortexField.velocity(p); };
    orbiting.respawnAll(inRing);
    const std::vector<Particle> start(orbiting.particles().begin(), orbiting.particles().end());
    constexpr int steps = 360;
    const double elapsed = steps * k_frame;
    for (int k = 0; k < steps; ++k) orbiting.step(k_frame, vortexVelocity, inRing);
    double orbit = 0.0;
    const std::span<const Particle> moved = orbiting.particles();
    for (std::size_t k = 0; k < moved.size(); ++k) {
        const auto [radius, theta] = polarAbout(start[k].pos, {});
        const double later = theta - elapsed / (radius * radius);
        const Vec2 exact{radius * std::cos(later), radius * std::sin(later)};
        orbit = std::max(orbit, length(moved[k].pos - exact) / radius);
        ++r.samples;
    }

    // A cylinder: each particle stays on its streamline, whatever the speed along it. Those that
    // left the box are back at age 0 and do not count.
    const FlowField cylinderField = modelOf(cylinderScene()).field;
    settings.bounds = {{-4.0, -2.5}, {4.0, 2.5}};
    ParticleSystem streaming(1000, settings);
    const auto outside = [](Vec2 p) { return length(p) > 1.05; };
    const auto cylinderVelocity = [&](Vec2 p) { return cylinderField.velocity(p); };
    streaming.respawnAll(outside);
    const std::vector<Particle> first(streaming.particles().begin(), streaming.particles().end());
    for (int k = 0; k < steps; ++k) streaming.step(k_frame, cylinderVelocity, outside);
    double psi = 0.0;
    int followed = 0;
    const std::span<const Particle> last = streaming.particles();
    for (std::size_t k = 0; k < last.size(); ++k) {
        if (last[k].age < elapsed) continue;
        psi = std::max(psi, std::abs(cylinderPsi(1.0, 1.0, 0.0, last[k].pos) - cylinderPsi(1.0, 1.0, 0.0, first[k].pos)));
        ++followed;
        ++r.samples;
    }
    r.error = std::max(orbit, psi);
    r.tolerance = 5e-3;
    r.errorLabel = "max position error per unit radius after 6 s around a vortex, and max change of ψ past a cylinder";
    r.comparisons = {{"Vortex: largest error, relative to r", num(orbit), "0"},
                     {"Cylinder: largest change of ψ", num(psi), "0"},
                     {"Particles followed past the cylinder", std::to_string(followed), "—"}};
    r.note = "Midpoint RK2 with the plot's limits on the step: the phase error is about (ωh)²/6 per radian turned, which "
             "comes to 0.4% for a particle at r = 0.5 over these 6 s. That sets the tolerance.";
    if (followed < 100) r.error = inf;
    return finish(std::move(r));
}

VerificationResult dyeStreaks() {
    VerificationResult r;
    DyeSettings settings;
    settings.bounds = k_view.expanded(0.25);
    settings.maxAge = 15.0;
    settings.maxStepLength = k_plotStep;

    // A uniform flow: every point at V∞ × its age from the emitter, a release interval apart, however
    // unevenly the frames fall: a frame is rarely a whole number of release intervals long.
    const double speed = 1.5;
    const Vec2 stream = exactUniform(speed, 30.0 * deg);
    const FlowField uniformField = modelOf({.elements = {UniformFlow{.speed = speed, .angle = 30.0 * deg}}}).field;
    const DyeEmitter emitter{.position = {-3.0, -1.5}};
    DyeSystem straight;
    straight.settings() = settings;
    const auto always = [](Vec2) { return true; };
    const auto uniformVelocity = [&](Vec2 p) { return uniformField.velocity(p); };
    const double frames[] = {1.0 / 60.0, 1.0 / 37.0, 1.0 / 20.0};
    for (int k = 0; k < 120; ++k) straight.step(frames[k % 3], std::span(&emitter, 1), uniformVelocity, always);
    if (straight.streaks().size() != 1) return unable("Expected one streak.");
    double position = 0.0;
    double gap = 0.0;
    const std::deque<DyePoint>& points = straight.streaks().front().points;
    for (std::size_t k = 0; k < points.size(); ++k) {
        position = std::max(position, length(points[k].pos - (emitter.position + stream * points[k].age)));
        if (k > 0) gap = std::max(gap, std::abs(length(points[k].pos - points[k - 1].pos) - speed * settings.releaseInterval));
        ++r.samples;
    }

    // A lifting cylinder: in a steady flow a streakline is the streamline through its emitter.
    const Scene scene = dyeScene();
    const FlowField cylinderField = modelOf(scene).field;
    std::vector<DyeEmitter> emitters;
    for (const Vec2 p : emitterPositions(scene.rakes.front())) emitters.push_back({.position = p});
    DyeSystem curved;
    curved.settings() = settings;
    const auto velocity = [&](Vec2 p) { return cylinderField.velocity(p); };
    const auto valid = [&](Vec2 p) { return cylinderField.inDomain(p) && !cylinderField.nearSingularity(p, 0.03); };
    for (int k = 0; k < 480; ++k) curved.step(k_frame, emitters, velocity, valid);
    double psi = 0.0;
    bool past = false; // some dye got past the cylinder
    for (std::size_t i = 0; i < curved.streaks().size() && i < emitters.size(); ++i) {
        const double start = cylinderPsi(1.0, 1.0, 2.0 * pi, emitters[i].position);
        for (const DyePoint& p : curved.streaks()[i].points) {
            if (!p.alive) continue;
            psi = std::max(psi, std::abs(cylinderPsi(1.0, 1.0, 2.0 * pi, p.pos) - start));
            past = past || p.pos.x > 0.0;
            ++r.samples;
        }
    }
    r.error = std::max({position, gap, psi});
    r.tolerance = 1e-3;
    r.errorLabel = "max position and gap errors in a uniform flow, and max change of ψ along the streaks past a cylinder";
    r.comparisons = {{"Uniform flow: position error", num(position), "0"},
                     {"Uniform flow: spacing of the dye", std::format("{:.6f}", length(points[1].pos - points[0].pos)),
                      std::format("V∞ Δt = {:.6f}", speed * settings.releaseInterval)},
                     {"Cylinder: largest change of ψ along a streak", num(psi), "0"}};
    r.note = "Dye moves by midpoint RK2 with the plot's step; its second-order error sets the tolerance. Because dye is "
             "released at a steady rate, its spacing shows the speed.";
    if (!past) r.error = inf;
    return finish(std::move(r));
}


// --- The list -----------------------------------------------------------------------------

constexpr const char* k_foundations = "Foundations";
constexpr const char* k_elementary = "Elementary flows";
constexpr const char* k_superposition = "Superposition and bodies";
constexpr const char* k_walls = "Symmetry walls";
constexpr const char* k_forces = "Pressure and forces";
constexpr const char* k_lines = "Streamlines and equipotentials";
constexpr const char* k_tracers = "Particles and dye";

const VerificationCheck k_checks[] = {
    {k_foundations, "Laplace's equation",
     "Potential flow is incompressible and irrotational: div V = 0 and curl V = 0 everywhere except at the singular points.",
     "Differentiate the velocity of a flow with every kind of element numerically, at 2,000 points.", mixedScene, laplace},
    {k_foundations, "Complex velocity is dF/dz",
     "W = u − iv is the derivative of the complex potential F = φ + iψ, so u = ∂φ/∂x = ∂ψ/∂y and v = ∂φ/∂y = −∂ψ/∂x "
     "(Cauchy–Riemann).",
     "Differentiate each element's F numerically and compare with the W the program uses for the velocity.", mixedScene,
     cauchyRiemann},
    {k_foundations, "Mass and circulation balance",
     "The volume flow out of any closed curve is the total source strength inside it, and the circulation around it is "
     "the total vortex strength inside it.",
     "Integrate v·n and v·t around six circles: around each element, around all of them and around none.", mixedScene,
     fluxAndCirculation},
    {k_foundations, "Superposition of every element",
     "The velocity of a sum of potential flows is the sum of their velocities.",
     "Compare every node of a flow with a uniform stream, a source, a sink, a vortex and a doublet with the five exact "
     "velocities added by hand.",
     mixedScene, superposition},
    {k_foundations, "Singular points stay finite",
     "A source, a vortex or a doublet has infinite velocity at its centre, so it must never be evaluated there: the "
     "velocity at the centre is that of the other elements.",
     "Read the velocity exactly at the centre of each element, and at the edge of a plate.", mixedScene, singularPoints},
    {k_foundations, "Turning or mirroring a scene",
     "The equations do not depend on the axes: turn every element about the origin and the flow turns with it; mirror "
     "them and it mirrors, with vortices and sector flows running the other way.",
     "Compare every node for three scenes (a wall, a channel and a sector, with every element), turned 33° and mirrored "
     "in the x axis.",
     inclinedWallScene, symmetries},
    {k_foundations, "Velocity integrates to φ and ψ",
     "dF = W dz, so φ(b) − φ(a) = ∫ V·dl along any path between the two points, and ψ(b) − ψ(a) = ∫ V·n ds is the volume "
     "flow across it.",
     "Integrate the velocity along five paths, two of them in a channel, and compare with the program's complex potential.",
     mixedScene, lineIntegrals},
    {k_foundations, "Bernoulli satisfies Euler's equation",
     "Steady inviscid flow obeys (V·∇)V = −∇p/ρ. Bernoulli's p + ½ρ|V|² = const is its integral when the flow is "
     "irrotational, so the pressure the program shows must drive the flow it shows.",
     "Differentiate the velocity and the Cp numerically at 2,000 points and compare both sides.", mixedScene, bernoulliEuler},

    {k_elementary, "Uniform flow",
     "V = V∞ (cos α, sin α) everywhere.", "Compare every node of the view for V∞ = 1, at α = 45° and α = 0°.", uniformScene,
     uniform},
    {k_elementary, "Source and sink",
     "V_r = Λ / 2πr, V_θ = 0: outwards for Λ > 0, inwards for Λ < 0.",
     "Compare every node for a source Λ = 2π at the origin and a sink Λ = −2π at (1, 0.5).", sourceSinkScene, sourceSink},
    {k_elementary, "Point vortex",
     "V_r = 0, V_θ = −Γ / 2πr: clockwise for Γ > 0 (Anderson's convention).",
     "Compare every node for Γ = 2π at the origin and Γ = −3 at (1, −0.5).", vortexScene, vortex},
    {k_elementary, "Doublet",
     "V_r = −μ cos(θ − β) / r², V_θ = −μ sin(θ − β) / r², with μ = κ / 2π.",
     "Compare every node for μ = 1 along x, and μ = 0.8 with its axis at β = 40°.", doubletScene, doublet},
    {k_elementary, "Flow in a corner",
     "F = V zⁿ with n = π / (θ₂ − θ₁): V_r = nV rⁿ⁻¹ cos nθ, V_θ = −nV rⁿ⁻¹ sin nθ, defined only inside the wedge.",
     "Compare every node for a 90° corner (n = 2) and a 90° wedge at 30° to 120°; the nodes outside must be masked.",
     cornerScene, corner},
    {k_elementary, "Flow around a plate edge",
     "The sector flow with walls at 0° and 360° (n = ½) is the flow around the edge of a semi-infinite flat plate; "
     "the speed grows as r^(−½) towards the edge.",
     "Compare every node, on both sides of the plate.", plateScene, plateEdge},

    {k_superposition, "Cylinder",
     "Uniform flow plus a doublet μ = V∞a² is the flow past a cylinder of radius a: "
     "V_r = V∞ (1 − a²/r²) cos θ, V_θ = −V∞ (1 + a²/r²) sin θ.",
     "Compare every node for a = 1, and for a = 0.7 in a stream V∞ = 2 at 30°; check that no flow crosses r = a.",
     cylinderScene, cylinder},
    {k_superposition, "Lifting cylinder",
     "Adding a vortex at the centre adds V_θ = −Γ / 2πr: faster over the top (|V| = 3V∞ for Γ = 2π) and slower "
     "underneath.",
     "Compare every node for a = 1, V∞ = 1, Γ = 2π.", liftingCylinderScene, liftingCylinder},
    {k_superposition, "Stagnation points on the cylinder",
     "The flow stops on the cylinder where sin θ = −Γ / 4πaV∞: at 0° and 180° without circulation, at 210° and 330° "
     "for Γ = 2π.",
     "Find the stagnation points the plot marks, for Γ = 0 and Γ = 2π.", liftingCylinderScene, cylinderStagnation},
    {k_superposition, "Stagnation points off the cylinder",
     "For Γ > 4πaV∞ the stagnation point leaves the cylinder and stands in the flow below it, at r = a (s + √(s² − 1)) with "
     "s = Γ / 4πaV∞, with its mirror inside. At Γ = 4πaV∞ the two meet at the bottom as a double root.",
     "Find the stagnation points for Γ = 6π and Γ = 4π.", offBodyStagnationScene, stagnationOffBody},
    {k_superposition, "Dividing streamlines",
     "Through a stagnation point run the streamlines that divide the flow: four at right angles at a simple point, six "
     "at 60° at a double one. Around a cylinder they all lie on ψ = 0.",
     "Trace the lines drawn when a stagnation point is clicked, on a cylinder, on a lifting cylinder, at a double point and "
     "between two equal sources (where they are the x and y axes).",
     liftingCylinderScene, dividingLines},
    {k_superposition, "Rankine oval",
     "A source and an equal sink in a uniform flow make a closed oval: the streamline ψ = 0, with stagnation points at "
     "x = ±√(b² + Λb/πV∞).",
     "Trace the body the way the surface-pressure view does, and check ψ at every point of it.", rankineOvalScene,
     rankineOval},
    {k_superposition, "Rankine half-body",
     "A source in a uniform flow makes an open body: the streamline ψ = Λ/2, starting Λ/2πV∞ upstream of the source "
     "and widening to Λ/V∞ far downstream.",
     "Trace the body and check ψ at every point of it.", halfBodyScene, halfBody},
    {k_superposition, "Source near a wall",
     "A wall is replaced by a mirror image of the source: the pair is the textbook solution, with no flow through the wall.",
     "Compare every node of the program's flow (the source, the wall and the image the program adds) with the source pair "
     "written out by hand.",
     sourceNearWallScene, sourceNearWall},

    {k_walls, "Inclined wall, every element",
     "A symmetry wall is a streamline: the velocity has no component normal to it anywhere along it.",
     "A wall at 25° with a uniform flow, a source, a vortex and a doublet; sample v·n along the wall.", inclinedWallScene,
     inclinedWall},
    {k_walls, "Corner of two walls",
     "Two walls at 90° need three images (−Γ, −Γ, +Γ for a vortex), and then close exactly.",
     "Compare the field with the four-vortex solution and sample v·n along both walls.", cornerWallsScene, cornerWalls},
    {k_walls, "Source in a channel",
     "Two parallel walls a distance h apart need two infinite rows of images, which sum to Λ/4h [coth π(z − z₀)/2h + "
     "coth π(z − z̄₀)/2h]. Half of the source's flow leaves each way, at speed Λ/2h far down the channel.",
     "Compare every node with the same sum in real form (sinh, cosh and sin only), the speed far down the channel, and "
     "v·n along both walls.",
     channelSourceScene, channelSource},
    {k_walls, "Inclined channel, every element",
     "The closed form of a channel is the limit of the series of images, for every kind of element.",
     "A channel at 25° with a uniform flow, a source, a vortex and a doublet: compare the flow with the direct sum over "
     "thousands of copies, and sample v·n along both walls.",
     inclinedChannelScene, inclinedChannel},
    {k_walls, "Corner of 60°, by conformal map",
     "ζ = z³ opens a 60° corner onto a half plane, where one image closes the flow, so the solution is known without any "
     "sum of images: W = (Λ/2π) 3z² [1/(ζ − ζ₀) + 1/(ζ − ζ̄₀)].",
     "Compare every node of a source and of a vortex with the program's five images, and sample v·n along both walls.",
     wedgeSourceScene, wedgeCorner},
    {k_walls, "Vortex in a channel",
     "A vortex midway between walls h apart has images of alternating sign, which sum to W = iΓ / (2h sinh(πz/h)). The speed "
     "on the wall is Γ/2h sech(πx/h), dying away within a few channel widths.",
     "Compare every node with the same sum in real form, and the speed along the wall with the secant.", channelVortexScene,
     channelVortex},

    {k_forces, "Pressure on a cylinder",
     "Bernoulli: Cp = 1 − |V|²/V∞² = 1 − (2 sin θ + Γ / 2πaV∞)² on the cylinder: 1 at the stagnation points, −3 at the "
     "top without circulation.",
     "Trace the surface and compare Cp at every point, with and without circulation.", surfacePressureScene,
     surfacePressure},
    {k_forces, "Pressure on the Rankine bodies",
     "The speed peaks at the top of an oval, at V∞ + (Λ/π) b / (b² + h²). On a half-body Cp depends only on the angle θ seen "
     "from the source, and the fastest speed is 1.26 V∞.",
     "Trace both surfaces and compare Cp at every point, and the lowest Cp.", halfBodyPressureScene, rankinePressure},
    {k_forces, "d'Alembert's paradox",
     "A closed body in steady, inviscid potential flow feels no drag.",
     "Integrate the surface pressure around a cylinder, a lifting cylinder and a Rankine oval.", dragScene, dAlembert},
    {k_forces, "Kutta–Joukowski lift",
     "The lift per unit span is L′ = ρ∞ V∞ Γ, so L′/q∞ = 2Γ/V∞, whatever the body's shape.",
     "Integrate the surface pressure for three lifting cylinders: at incidence, with negative Γ, and with Γ > 4πaV∞.",
     liftScene, kuttaJoukowski},
    {k_forces, "Forces from the far field",
     "Blasius' theorem gives the force from the flow on any curve around the body: F_x − iF_y = (iρ/2) ∮ W² dz = −ρ W∞ "
     "(ΣΛ + iΣΓ). So a net source feels −ρV∞Λ along the stream (Lagally) and circulation feels the Kutta–Joukowski lift.",
     "Integrate W² on a circle around four sets of elements. No surface is traced, so this is independent of the "
     "surface-pressure checks.",
     liftScene, blasius},
    {k_forces, "Forces along a parameter sweep",
     "Lift is linear in the circulation, L′/q∞ = 2Γ/V∞, the drag stays 0, and the stagnation points slide around the "
     "cylinder as sin θ = −Γ / 4πaV∞.",
     "Sweep Γ from −3π to 3π with the code behind the sweep export, and compare what it would write.", surfacePressureScene,
     sweepForces},

    {k_lines, "Uniform flow streamlines",
     "The streamlines of a uniform flow are straight lines along α, running downstream.",
     "Trace the lines the plot draws at α = 0° and 45°, and check every step's direction and every arrow.", uniformScene,
     uniformLines},
    {k_lines, "Source and sink streamlines are rays",
     "The streamlines of a source or a sink are straight rays from its centre: outwards for Λ > 0, inwards for Λ < 0.",
     "Trace the lines the plot draws for both signs, and check that each stays on its ray and runs the right way.",
     radialScene, radialLines},
    {k_lines, "Vortex streamlines are circles",
     "Around a point vortex the streamlines are the circles r = const, clockwise for Γ > 0 and counter-clockwise for Γ < 0.",
     "Trace the streamlines the plot draws, measure how far each strays from its starting radius, and check which way the "
     "arrows turn for both signs.",
     vortexScene, vortexCircles},
    {k_lines, "ψ is constant along streamlines",
     "A streamline is a line of constant ψ. Around a lifting cylinder ψ = V∞ (r − a²/r) sin θ + (Γ/2π) ln(r/a); in a 90° "
     "corner ψ = 2Vxy, so the streamlines are hyperbolas.",
     "Evaluate ψ along every traced line, in both flows and for lines kept from frame to frame during a sweep.",
     liftingCylinderScene, streamFunction},
    {k_lines, "Streamlines are evenly spaced",
     "Seeded by the method of Jobard and Lefer, lines stay at least half a separation d apart and leave no gap wider than d. "
     "In a uniform flow they are exactly d apart.",
     "Measure the spacing in a uniform flow at 30°, and the closest lines and the widest gap around a lifting cylinder.",
     liftingCylinderScene, evenSpacing},
    {k_lines, "Arrows show the flow direction",
     "The arrows on the streamlines must point along the local velocity, never upstream, and be evenly spaced along a line.",
     "Compare every arrow around a lifting cylinder with the velocity at its foot, and measure the spacing on straight lines.",
     liftingCylinderScene, arrows},
    {k_lines, "Equipotentials cross streamlines at 90°",
     "V = grad φ, so lines of constant φ are everywhere perpendicular to the flow.",
     "Trace the equipotential lines the plot draws around a cylinder; check the angle of every step and φ along each line.",
     equipotentialScene, equipotentials},

    {k_tracers, "Tracer particles move with the flow",
     "A particle at radius r around a vortex Γ turns at the angular speed Γ / 2πr², so its position at any time is known "
     "exactly. In any steady flow a particle stays on its streamline.",
     "Run the particle system for 6 s around a vortex, and past a cylinder, and compare the positions with the exact orbit "
     "and with ψ.",
     vortexScene, particles},
    {k_tracers, "Dye streaks are streamlines",
     "In a steady flow a streakline is the streamline through its emitter. In a uniform flow the dye is a distance V∞Δt "
     "apart for a release interval Δt, so the spacing shows the speed.",
     "Run the dye system in a uniform flow and past a lifting cylinder, and check positions, spacing and ψ along the streaks.",
     dyeScene, dyeStreaks},
};

} // namespace

std::span<const VerificationCheck> verificationChecks() { return k_checks; }

} // namespace pfl
