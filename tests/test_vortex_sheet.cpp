#include "core/parameters.hpp"
#include "core/scene.hpp"
#include "core/streamline.hpp"
#include "core/surface_pressure.hpp"
#include "core/sweep.hpp"
#include "core/vortex_sheet.hpp"
#include "io/data_export.hpp"
#include "io/scene_io.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <vector>

using namespace pfl;

namespace {

constexpr double pi = std::numbers::pi;
constexpr double deg = pi / 180.0;
const Rect k_view{{-4.0, -2.5}, {4.0, 2.5}};

bool contains(const std::vector<std::string>& warnings, const std::string& text) {
    return std::ranges::any_of(warnings, [&](const std::string& w) { return w.find(text) != std::string::npos; });
}

// An airfoil of chord 2 at the origin at angle of attack alpha in a unit stream along +x.
Scene airfoilScene(CamberLine camber, double alpha, int panels = 60) {
    return {.elements = {UniformFlow{}, ThinAirfoil{.camber = camber, .chord = 2.0, .angle = -alpha, .panels = panels}}};
}

bool segmentsCross(Vec2 a, Vec2 b, Vec2 c, Vec2 d) {
    const double o1 = cross(b - a, c - a);
    const double o2 = cross(b - a, d - a);
    const double o3 = cross(d - c, a - c);
    const double o4 = cross(d - c, b - c);
    return ((o1 > 0.0) != (o2 > 0.0)) && ((o3 > 0.0) != (o4 > 0.0));
}

} // namespace

TEST_CASE("thin airfoil parameters are validated") {
    CHECK_FALSE(thinAirfoilError(ThinAirfoil{}));
    CHECK(thinAirfoilError(ThinAirfoil{.chord = 0.0}));
    CHECK(thinAirfoilError(ThinAirfoil{.panels = k_minSheetPanels - 1}));
    CHECK(thinAirfoilError(ThinAirfoil{.panels = k_maxSheetPanels + 1}));
    CHECK(thinAirfoilError(ThinAirfoil{.angle = std::nan("")}));
    CHECK(thinAirfoilError(ThinAirfoil{.nacaCamber = 10}));
    CHECK(thinAirfoilError(ThinAirfoil{.nacaPosition = 0}));
    CHECK_FALSE(thinAirfoilError(ThinAirfoil{.nacaCamber = 0, .nacaPosition = 0})); // a flat line needs no position
    CHECK(thinAirfoilError(ThinAirfoil{.camber = CamberLine::Arc, .arcCamber = 0.3}));
    CHECK_FALSE(thinAirfoilError(ThinAirfoil{.camber = CamberLine::Arc, .arcCamber = -0.1}));
    CHECK(airfoilName(ThinAirfoil{}) == "NACA 24xx mean line");
    CHECK(airfoilName(ThinAirfoil{.camber = CamberLine::Flat}) == "flat plate");
}

TEST_CASE("the camber line: vertices from the leading to the trailing edge, cosine-spaced") {
    const ThinAirfoil airfoil{.position = {0.5, -0.2}, .chord = 3.0, .angle = 20.0 * deg, .panels = 40};
    const std::vector<Vec2> v = sheetVertices(airfoil);
    REQUIRE(v.size() == 41);
    const Vec2 e = chordDirection(airfoil);
    CHECK(length(v.front() - (airfoil.position - e * 1.5)) < 1e-14);
    CHECK(length(v.back() - (airfoil.position + e * 1.5)) < 1e-14);
    // Short panels at both edges, the longest at mid-chord.
    CHECK(length(v[1] - v[0]) < 0.05 * length(v[21] - v[20]));
    CHECK(length(v[40] - v[39]) < 0.05 * length(v[21] - v[20]));

    // NACA 2412: largest camber 2% of the chord at 40%; a circular arc: h at mid-chord.
    CHECK(camberHeight(ThinAirfoil{}, 0.4) == doctest::Approx(0.02));
    CHECK(camberSlope(ThinAirfoil{}, 0.4) == doctest::Approx(0.0));
    const ThinAirfoil arc{.camber = CamberLine::Arc, .arcCamber = 0.05};
    CHECK(camberHeight(arc, 0.5) == doctest::Approx(0.05));
    CHECK(std::abs(camberHeight(arc, 0.0)) < 1e-15);
    CHECK(std::abs(camberHeight(arc, 1.0)) < 1e-15);
    for (const double x : {0.1, 0.3, 0.7}) {
        const double h = 1e-6;
        CHECK(camberSlope(arc, x) == doctest::Approx((camberHeight(arc, x + h) - camberHeight(arc, x - h)) / (2.0 * h)).epsilon(1e-7));
    }
    CHECK(angleOfAttack(ThinAirfoil{.angle = -5.0 * deg}, {1.0, 0.0}) == doctest::Approx(5.0 * deg));
}

TEST_CASE("thin-airfoil theory: flat plate, NACA mean line and circular arc") {
    const ThinAirfoilTheory plate = thinAirfoilTheory(ThinAirfoil{.camber = CamberLine::Flat}, 5.0 * deg);
    CHECK(plate.A[0] == doctest::Approx(5.0 * deg));
    CHECK(plate.A[1] == 0.0);
    CHECK(plate.cl == doctest::Approx(2.0 * pi * 5.0 * deg));
    CHECK(plate.cmQuarter == 0.0);
    CHECK(plate.alphaZeroLift == 0.0);
    for (const double x : {0.1, 0.5, 0.9}) CHECK(plate.deltaCp(x) == doctest::Approx(4.0 * 5.0 * deg * std::sqrt((1.0 - x) / x)));
    CHECK(plate.deltaCp(1.0) == 0.0);

    // NACA 2412 (Anderson, example 4.6): alpha_L0 = -2.077 deg, c_m,c/4 = -0.0532.
    const ThinAirfoilTheory naca = thinAirfoilTheory(ThinAirfoil{}, 4.0 * deg);
    CHECK(naca.alphaZeroLift / deg == doctest::Approx(-2.077).epsilon(1e-3));
    CHECK(naca.cmQuarter == doctest::Approx(-0.0532).epsilon(3e-3));
    CHECK(naca.cl == doctest::Approx(2.0 * pi * (4.0 * deg - naca.alphaZeroLift)));
    // The closed form against the trapezoidal rule in theta0 (second order across the kink at x = p).
    constexpr int samples = 200000;
    double a1 = 0.0;
    for (int k = 0; k <= samples; ++k) {
        const double theta = pi * k / samples;
        const double w = (k == 0 || k == samples) ? 0.5 : 1.0;
        a1 += w * camberSlope(ThinAirfoil{}, 0.5 * (1.0 - std::cos(theta))) * std::cos(theta);
    }
    a1 *= 2.0 / samples; // (2/pi) * pi/samples
    CHECK(naca.A[1] == doctest::Approx(a1).epsilon(1e-8));

    // A shallow arc is nearly the parabola z = 4h x(1 - x): alpha_L0 = -2h/c, c_m,c/4 = -pi h/c.
    const ThinAirfoilTheory arc = thinAirfoilTheory(ThinAirfoil{.camber = CamberLine::Arc, .arcCamber = 0.01}, 0.0);
    CHECK(arc.alphaZeroLift == doctest::Approx(-0.02).epsilon(1e-3));
    CHECK(arc.cmQuarter == doctest::Approx(-pi * 0.01).epsilon(1e-3));
}

TEST_CASE("vortex sheet: Kutta condition, no flow through the line, and the flat plate's lift") {
    const FlowModel model = buildFlowModel(airfoilScene(CamberLine::Flat, 5.0 * deg, 80));
    REQUIRE(model.field.vortexSheets.size() == 1);
    const VortexSheet& sheet = model.field.vortexSheets.front();
    CHECK(sheet.gamma.size() == 81);
    CHECK(sheet.gamma.back() == 0.0);
    CHECK(sheetNormalFlowResidual(sheet, model.field) < 1e-13);
    double gamma = 0.0;
    for (std::size_t j = 0; j < sheet.size(); ++j) gamma += sheet.lengths[j] * sheet.strengthAt(j);
    CHECK(sheet.circulation == doctest::Approx(gamma));
    REQUIRE(sheet.unitCirculation);

    const std::optional<AirfoilLoading> loading = airfoilLoading(model.field, 0, model.freeStream);
    REQUIRE(loading);
    CHECK(loading->alpha == doctest::Approx(5.0 * deg));
    CHECK(loading->clKuttaJoukowski == doctest::Approx(2.0 * pi * std::sin(5.0 * deg)).epsilon(1e-4));
    CHECK(loading->clPressure == doctest::Approx(loading->clKuttaJoukowski).epsilon(1e-12));
    CHECK(std::abs(loading->cmQuarter) < 1e-4);
    CHECK(std::abs(*loading->alphaZeroLift) < 1e-12);
    CHECK(*loading->liftSlope == doctest::Approx(2.0 * pi * std::cos(5.0 * deg)).epsilon(1e-4));
    // The upper side is faster: suction above.
    CHECK(loading->cpUpper[40] < loading->cpLower[40]);
    CHECK(loading->deltaCp[40] > 0.0);
}

TEST_CASE("vortex sheet: the velocity jumps by gamma across it, and the far-field series joins the panel sums") {
    const FlowModel model = buildFlowModel(airfoilScene(CamberLine::Naca, 4.0 * deg, 50));
    const VortexSheet& sheet = model.field.vortexSheets.front();
    for (const std::size_t i : {std::size_t{5}, std::size_t{25}, std::size_t{45}}) {
        const Vec2 off = sheet.normals[i] * (1e-7 * sheet.lengths[i]);
        const Vec2 above = model.field.velocity(sheet.controlPoints[i] + off);
        const Vec2 below = model.field.velocity(sheet.controlPoints[i] - off);
        CHECK(dot(above - below, sheet.tangents[i]) == doctest::Approx(sheet.strengthAt(i)).epsilon(1e-6));
        CHECK(dot(above - below, sheet.normals[i]) == doctest::Approx(0.0).epsilon(1e-6));
        // On the sheet: the mean of the two sides.
        const Vec2 mean = model.field.velocity(sheet.controlPoints[i]);
        CHECK(length(mean - (above + below) * 0.5) < 1e-6);
    }
    // Either side of the switch to the Laurent series, against the panels summed one by one.
    for (const double r : {1.999, 2.001}) {
        for (const double t : {0.3, 2.0, 4.0}) {
            const Vec2 p = sheet.centre + Vec2{std::cos(t), std::sin(t)} * (r * sheet.radius);
            Complex direct{};
            for (std::size_t j = 0; j < sheet.size(); ++j) {
                const PanelInfluence w = unitVortexPanelVelocity(sheet.vertices[j], sheet.vertices[j + 1], p);
                direct += w.first * sheet.gamma[j] + w.second * sheet.gamma[j + 1];
            }
            CHECK(std::abs(sheet.complexVelocity(p) - direct) < 1e-12);
        }
    }
}

TEST_CASE("vortex sheet: dF/dz = W, and phi jumps by Gamma only across the sheet and the wake ray") {
    const FlowModel model = buildFlowModel(airfoilScene(CamberLine::Naca, 4.0 * deg, 40));
    const VortexSheet& sheet = model.field.vortexSheets.front();
    for (const Vec2 p : {Vec2{0.3, 0.5}, Vec2{-1.2, -0.3}, Vec2{0.5, -0.1}, Vec2{-0.9, 0.1}, Vec2{2.5, 1.0}, Vec2{-3.0, -2.0}}) {
        const double h = 1e-5;
        const Complex derivative = (sheet.complexPotential(p + Vec2{h, 0.0}) - sheet.complexPotential(p - Vec2{h, 0.0})) / (2.0 * h);
        CHECK(std::abs(derivative - sheet.complexVelocity(p)) < 1e-7);
    }
    const auto phi = [&](Vec2 p) { return sheet.complexPotential(p).real(); };
    const Vec2 te = sheet.trailingEdge();
    const Vec2 wake = chordDirection(sheet.airfoil);
    const Vec2 across = perp(wake) * 1e-8;
    // Behind the trailing edge, on the ray: phi above minus below is Gamma, the circulation of a loop
    // from below the ray round the leading edge to above it.
    const Vec2 behind = te + wake * 0.7;
    CHECK(phi(behind + across) - phi(behind - across) == doctest::Approx(sheet.circulation).epsilon(1e-6));
    // Ahead of the leading edge and beside the line, continuous.
    for (const Vec2 p : {sheet.leadingEdge() - wake * 0.5, Vec2{0.0, 1.0}, Vec2{0.0, -1.0}}) {
        CHECK(std::abs(phi(p + across) - phi(p - across)) < 1e-6);
    }
}

TEST_CASE("vortex sheet: solved with source-panel bodies; alpha_L0 only in a uniform stream") {
    Scene scene = airfoilScene(CamberLine::Naca, 4.0 * deg);
    const double alone = buildFlowModel(scene).field.vortexSheets.front().circulation;
    scene.elements.push_back(PanelBody{.position = {0.0, -1.2}, .semiAxisA = 1.5, .semiAxisB = 0.3, .panels = 60});
    const FlowModel model = buildFlowModel(scene);
    REQUIRE(model.field.vortexSheets.size() == 1);
    REQUIRE(model.field.panelBodies.size() == 1);
    CHECK(model.field.vortexSheets.front().circulation != doctest::Approx(alone).epsilon(1e-3));
    CHECK(model.field.vortexSheets.front().unitCirculation); // still only uniform flow outside the panels
    CHECK(sheetNormalFlowResidual(model.field.vortexSheets.front(), model.field) < 1e-12);
    CHECK(normalFlowResidual(model.field.panelBodies.front(), model.field) < 1e-12);

    scene.elements.push_back(SourceSink{.strength = 1.0, .position = {-3.0, 1.5}});
    const FlowModel withSource = buildFlowModel(scene);
    CHECK_FALSE(withSource.field.vortexSheets.front().unitCirculation);
    CHECK_FALSE(airfoilLoading(withSource.field, 0, withSource.freeStream)->alphaZeroLift);
}

TEST_CASE("vortex sheet: a mirrored airfoil has the opposite circulation") {
    const Scene upright{.elements = {UniformFlow{.speed = 1.0, .angle = 10.0 * deg},
                                     ThinAirfoil{.position = {0.3, 0.4}, .angle = 5.0 * deg}}};
    const Scene mirrored{.elements = {UniformFlow{.speed = 1.0, .angle = -10.0 * deg},
                                      ThinAirfoil{.position = {0.3, -0.4}, .angle = -5.0 * deg, .nacaCamber = -2}}};
    CHECK(buildFlowModel(mirrored).field.vortexSheets.front().circulation ==
          doctest::Approx(-buildFlowModel(upright).field.vortexSheets.front().circulation));
}

TEST_CASE("vortex sheet: streamlines do not cut through it, and its stagnation point is under the leading edge") {
    const FlowModel model = buildFlowModel(airfoilScene(CamberLine::Flat, 8.0 * deg, 40));
    const VortexSheet& sheet = model.field.vortexSheets.front();
    StreamlineSettings settings;
    settings.bounds = k_view;
    settings.separation = k_view.size().y / 30.0;
    settings.singularRadius = 0.01;
    const StreamlineSet lines = computeStreamlines(model.field, settings);
    REQUIRE(lines.lines.size() > 10);
    for (const std::vector<Vec2>& line : lines.lines) {
        for (std::size_t k = 0; k + 1 < line.size(); ++k) CHECK_FALSE(sheet.crossedBy(line[k], line[k + 1]));
    }
    CHECK(model.field.nearSingularity(sheet.leadingEdge(), 1e-6));
    CHECK_FALSE(model.field.nearSingularity(sheet.trailingEdge(), 1e-6));

    const std::vector<Vec2> stagnation = findStagnationPoints(model.field, k_view, 1.0);
    REQUIRE(stagnation.size() == 1);
    const Vec2 local = stagnation.front() - sheet.leadingEdge();
    const double along = dot(local, chordDirection(sheet.airfoil)) / sheet.airfoil.chord;
    CHECK(along > 0.0);
    CHECK(along < 0.05); // the exact plate: x/c = sin² α ≈ 0.019
    CHECK(dot(local, perp(chordDirection(sheet.airfoil))) < 0.0); // below
}

TEST_CASE("vortex sheet: scene files, parameters and warnings") {
    Scene scene{.elements = {UniformFlow{},
                             ThinAirfoil{.camber = CamberLine::Arc, .position = {1.0, 0.5}, .chord = 1.5, .angle = -3.0 * deg,
                                         .panels = 33, .arcCamber = -0.04},
                             ThinAirfoil{.camber = CamberLine::Naca, .panels = 20, .nacaCamber = 4, .nacaPosition = 3}}};
    const SceneLoadResult loaded = sceneFromJson(sceneToJson(scene));
    REQUIRE(loaded.scene);
    for (std::size_t i = 1; i < 3; ++i) {
        const auto& a = std::get<ThinAirfoil>(scene.elements[i]);
        const auto& b = std::get<ThinAirfoil>(loaded.scene->elements[i]);
        CHECK(a.camber == b.camber);
        CHECK(a.position == b.position);
        CHECK(a.chord == b.chord);
        CHECK(a.angle == doctest::Approx(b.angle));
        CHECK(a.panels == b.panels);
        CHECK(a.arcCamber == b.arcCamber);
        CHECK(a.nacaCamber == b.nacaCamber);
        CHECK(a.nacaPosition == b.nacaPosition);
    }
    CHECK_FALSE(sceneFromJson(R"({"format": "pfl-scene", "version": 1, "elements": [{"type": "vortex_sheet", "camber": "wavy",
        "x": 0, "y": 0, "chord": 1, "panels": 10}]})").scene);

    // Sweeping alpha = -theta.
    const std::span<const ElementParameter> parameters = elementParameters(ElementKind::ThinAirfoil);
    const auto alpha = std::ranges::find_if(parameters, [](const ElementParameter& p) { return std::string(p.label) == "α"; });
    REQUIRE(alpha != parameters.end());
    SceneElement element = ThinAirfoil{.angle = 0.1};
    CHECK(alpha->get(element) == doctest::Approx(-0.1));
    alpha->set(element, 0.2);
    CHECK(std::get<ThinAirfoil>(element).angle == doctest::Approx(-0.2));

    // Warnings: no stream; an airfoil into a body; two airfoils crossing.
    const Scene still{.elements = {ThinAirfoil{}}};
    CHECK(contains(sceneWarnings(still, buildFlowModel(still, false)), "no uniform flow"));
    const Scene overlap{.elements = {UniformFlow{}, ThinAirfoil{}, PanelBody{.position = {0.9, 0.0}, .semiAxisA = 0.3, .semiAxisB = 0.3}}};
    CHECK(contains(sceneWarnings(overlap, buildFlowModel(overlap, false)), "runs into source-panel body"));
    const Scene crossing{.elements = {UniformFlow{}, ThinAirfoil{}, ThinAirfoil{.angle = pi / 2.0}}};
    CHECK(contains(sceneWarnings(crossing, buildFlowModel(crossing, false)), "cross each other"));
    const Scene steep = airfoilScene(CamberLine::Flat, 20.0 * deg);
    CHECK(contains(sceneWarnings(steep, buildFlowModel(steep, false)), "stalls"));
}

TEST_CASE("sweeping an airfoil's angle of attack: a c_l-alpha line in the sweep CSV") {
    const Scene scene = airfoilScene(CamberLine::Naca, 0.0, 40);
    const auto alpha = std::ranges::find_if(elementParameters(ElementKind::ThinAirfoil),
                                            [](const ElementParameter& p) { return std::string(p.label) == "α"; });
    const ParameterRef ref{1, static_cast<std::size_t>(alpha - elementParameters(ElementKind::ThinAirfoil).begin())};
    std::vector<SweepSample> samples;
    for (const double value : sweepValues(-2.0 * deg, 6.0 * deg, 5)) {
        const std::optional<SweepSample> sample = sampleSweep(scene, ref, value, k_view);
        REQUIRE(sample);
        REQUIRE(sample->airfoils.size() == 1);
        const SweepAirfoil& a = sample->airfoils.front();
        CHECK(a.alpha == doctest::Approx(value));
        CHECK(a.exactLift / a.chord == doctest::Approx(a.clTheory).epsilon(0.01)); // the panels near thin-airfoil theory
        CHECK(a.cmQuarter == doctest::Approx(a.cmTheory).epsilon(0.02));
        samples.push_back(*sample);
    }
    // Linear in alpha: equal steps of c_l.
    const double step = samples[1].airfoils[0].exactLift - samples[0].airfoils[0].exactLift;
    CHECK(samples[4].airfoils[0].exactLift - samples[3].airfoils[0].exactLift == doctest::Approx(step).epsilon(0.01));
    const std::string csv = sweepCsv(samples);
    std::size_t rows = 0;
    for (std::size_t at = csv.find(",airfoil,"); at != std::string::npos; at = csv.find(",airfoil,", at + 1)) ++rows;
    CHECK(rows == 5);
}

TEST_CASE("an airfoil makes no closed body of its own, and makes every other body non-isolated") {
    const Scene scene{.elements = {UniformFlow{}, Doublet{.strength = 0.25, .position = {0.0, -1.6}}, ThinAirfoil{.position = {0.0, 0.8}}}};
    const FlowModel model = buildFlowModel(scene);
    const std::vector<BodySurface> bodies = findBodySurfaces(model, k_view);
    REQUIRE(bodies.size() == 1); // the cylinder
    CHECK_FALSE(bodyForces(bodies.front(), model).isolated);
    CHECK_FALSE(cylinderReference(bodies.front(), model));
}
