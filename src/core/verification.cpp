#include "core/verification.hpp"

#include "core/flow_element.hpp"
#include "core/flow_field.hpp"
#include "core/image_system.hpp"
#include "core/streamline.hpp"
#include "core/surface_pressure.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <format>
#include <functional>
#include <limits>
#include <numbers>
#include <optional>
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

// Exact sum of the elements of mixedScene().
Vec2 exactMixed(Vec2 p) {
    return exactUniform(1.0, 15.0 * deg) + exactSource(2.5, {-1.5, 0.6}, p) + exactSource(-1.8, {1.6, -0.4}, p) +
           exactVortex(3.0, {0.2, 1.1}, p) + exactDoublet(0.3, 35.0 * deg, {0.4, -1.2}, p);
}

// --- Foundations ------------------------------------------------------------------------

VerificationResult laplace() {
    const FlowField field = modelOf(mixedScene()).field;
    const std::vector<Vec2> singular = {{-1.5, 0.6}, {1.6, -0.4}, {0.2, 1.1}, {0.4, -1.2}};
    constexpr double h = 1e-5;
    double divergence = 0.0;
    double curl = 0.0;
    int samples = 0;
    std::uint32_t state = 20'240'601u; // a fixed sequence, so every run checks the same points
    const auto uniform = [&state] {
        state = state * 1'664'525u + 1'013'904'223u;
        return static_cast<double>(state >> 8) / static_cast<double>(1u << 24);
    };
    while (samples < 2000) {
        const Vec2 p{k_view.min.x + k_view.size().x * uniform(), k_view.min.y + k_view.size().y * uniform()};
        if (std::ranges::any_of(singular, [p](Vec2 s) { return length(p - s) < 0.25; })) continue;
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

// --- Elementary flows -------------------------------------------------------------------

VerificationResult uniform() {
    const FlowField field = modelOf(uniformScene()).field;
    const Vec2 v = exactUniform(1.0, 45.0 * deg);
    const GridStats g = compareGrid(field, [v](Vec2) { return value(v); });
    return gridResult(g, {velocityAt(field, {3.0, -2.0}, v), velocityAt(field, {-1.0, 1.5}, v)});
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
    double lo = 0.0;
    double hi = 5.0;
    for (int k = 0; k < 200; ++k) {
        const double mid = 0.5 * (lo + hi);
        (mid - (pi - 2.0 * std::atan(mid)) > 0.0 ? hi : lo) = mid;
    }
    const BodyForces forces = bodyForces(*body, model);

    VerificationResult r;
    r.error = std::max(psiError, stagnation);
    r.tolerance = 1e-6;
    r.errorLabel = "max |ψ| over the traced surface (should be 0), and the stagnation-point error, per unit b";
    r.samples = static_cast<int>(body->points.size());
    r.comparisons = {{"Stagnation points x", std::format("±{:.9f}", std::abs(body->stagnationPoints[0].x)), std::format("±√3 = ±{:.9f}", std::sqrt(3.0))},
                     {"Chord", std::format("{:.6f}", forces.chord), std::format("2√3 = {:.6f}", 2.0 * std::sqrt(3.0))},
                     {"Half-height h", std::format("{:.6f}", halfHeight), std::format("{:.6f}", 0.5 * (lo + hi))},
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

VerificationResult vortexCircles() {
    const FlowField field = modelOf(vortexScene()).field;
    StreamlineSettings settings;
    settings.bounds = {{-2.0, -2.0}, {2.0, 2.0}};
    settings.separation = 0.3;
    const StreamlineSet set = computeStreamlines(field, settings);
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
                     {"Largest radius drift", num(r.error), "0 (circles r = const)"}};
    r.note = "RK4 in arc-length steps; its truncation error over one loop sets the tolerance.";
    if (closed < 3) r.error = inf;
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

// --- The list -----------------------------------------------------------------------------

constexpr const char* k_foundations = "Foundations";
constexpr const char* k_elementary = "Elementary flows";
constexpr const char* k_superposition = "Superposition and bodies";
constexpr const char* k_walls = "Symmetry walls";
constexpr const char* k_forces = "Pressure and forces";
constexpr const char* k_lines = "Streamlines and equipotentials";

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

    {k_elementary, "Uniform flow",
     "V = V∞ (cos α, sin α) everywhere.", "Compare every node of the view for V∞ = 1, α = 45°.", uniformScene, uniform},
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

    {k_forces, "Pressure on a cylinder",
     "Bernoulli: Cp = 1 − |V|²/V∞² = 1 − (2 sin θ + Γ / 2πaV∞)² on the cylinder: 1 at the stagnation points, −3 at the "
     "top without circulation.",
     "Trace the surface and compare Cp at every point, with and without circulation.", surfacePressureScene,
     surfacePressure},
    {k_forces, "d'Alembert's paradox",
     "A closed body in steady, inviscid potential flow feels no drag.",
     "Integrate the surface pressure around a cylinder, a lifting cylinder and a Rankine oval.", dragScene, dAlembert},
    {k_forces, "Kutta–Joukowski lift",
     "The lift per unit span is L′ = ρ∞ V∞ Γ, so L′/q∞ = 2Γ/V∞, whatever the body's shape.",
     "Integrate the surface pressure for three lifting cylinders: at incidence, with negative Γ, and with Γ > 4πaV∞.",
     liftScene, kuttaJoukowski},

    {k_lines, "Vortex streamlines are circles",
     "Around a point vortex the streamlines are the circles r = const.",
     "Trace the streamlines the plot draws and measure how far each strays from its starting radius.", vortexScene,
     vortexCircles},
    {k_lines, "Equipotentials cross streamlines at 90°",
     "V = grad φ, so lines of constant φ are everywhere perpendicular to the flow.",
     "Trace the equipotential lines the plot draws around a cylinder; check the angle of every step and φ along each line.",
     equipotentialScene, equipotentials},
};

} // namespace

std::span<const VerificationCheck> verificationChecks() { return k_checks; }

} // namespace pfl
