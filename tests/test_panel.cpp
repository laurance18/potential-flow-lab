#include "core/panel_method.hpp"
#include "core/parameters.hpp"
#include "core/scene.hpp"
#include "core/surface_pressure.hpp"
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

double signedArea(const std::vector<Vec2>& v) {
    double twice = 0.0;
    for (std::size_t k = 0; k < v.size(); ++k) twice += cross(v[k], v[(k + 1) % v.size()]);
    return 0.5 * twice;
}

bool contains(const std::vector<std::string>& warnings, const std::string& text) {
    return std::ranges::any_of(warnings, [&](const std::string& w) { return w.find(text) != std::string::npos; });
}

const PanelBody k_ellipse{.position = {0.4, -0.3}, .semiAxisA = 1.5, .semiAxisB = 0.6, .angle = 25.0 * deg, .panels = 60};

Scene ellipseScene() { return {.elements = {UniformFlow{.speed = 1.2, .angle = 10.0 * deg}, k_ellipse}}; }

} // namespace

TEST_CASE("panel vertices go clockwise around the shape") {
    SUBCASE("ellipse") {
        const std::vector<Vec2> v = panelVertices(k_ellipse);
        REQUIRE(v.size() == 60);
        CHECK(signedArea(v) < 0.0);
        for (const Vec2 p : v) {
            const Vec2 d = p - k_ellipse.position;
            const double c = std::cos(-k_ellipse.angle);
            const double s = std::sin(-k_ellipse.angle);
            const Vec2 q{c * d.x - s * d.y, s * d.x + c * d.y};
            CHECK(std::pow(q.x / 1.5, 2) + std::pow(q.y / 0.6, 2) == doctest::Approx(1.0));
        }
        // The first control point faces -x in the body's frame (the front, for a stream along its axis).
        const std::vector<Vec2> circle = panelVertices(PanelBody{.panels = 8});
        const Vec2 front = (circle[0] + circle[1]) * 0.5;
        CHECK(front.x == doctest::Approx(-std::cos(pi / 8.0)));
        CHECK(front.y == doctest::Approx(0.0).epsilon(1e-12));
    }
    SUBCASE("rectangle") {
        const PanelBody body{.shape = PanelShape::Rectangle, .semiAxisA = 2.0, .semiAxisB = 0.5, .panels = 20};
        const std::vector<Vec2> v = panelVertices(body);
        REQUIRE(v.size() == 20);
        CHECK(signedArea(v) == doctest::Approx(-4.0 * 2.0 * 0.5));
        for (const Vec2 corner : {Vec2{-2.0, -0.5}, Vec2{-2.0, 0.5}, Vec2{2.0, 0.5}, Vec2{2.0, -0.5}}) {
            CHECK(std::ranges::any_of(v, [&](Vec2 p) { return length(p - corner) < 1e-12; }));
        }
        for (const Vec2 p : v) CHECK((std::abs(std::abs(p.x) - 2.0) < 1e-12 || std::abs(std::abs(p.y) - 0.5) < 1e-12));
    }
}

TEST_CASE("panel body parameters are validated") {
    CHECK_FALSE(panelBodyError(PanelBody{}));
    CHECK(panelBodyError(PanelBody{.semiAxisA = 0.0}));
    CHECK(panelBodyError(PanelBody{.semiAxisB = -1.0}));
    CHECK(panelBodyError(PanelBody{.panels = k_minPanels - 1}));
    CHECK(panelBodyError(PanelBody{.panels = k_maxPanels + 1}));
    CHECK(panelBodyError(PanelBody{.shape = PanelShape::Rectangle, .panels = 21}));
    CHECK_FALSE(panelBodyError(PanelBody{.shape = PanelShape::Rectangle, .panels = 22}));
    CHECK(panelBodyError(PanelBody{.angle = std::nan("")}));
    CHECK(sceneElementError(SceneElement{PanelBody{.panels = 2}}));
}

TEST_CASE("the panel count is a whole number, and even on a rectangle") {
    const ElementParameter& n = elementParameters(ElementKind::PanelBody).back();
    REQUIRE(n.kind == ParameterKind::Count);
    SceneElement ellipse = PanelBody{};
    n.set(ellipse, 33.4);
    CHECK(std::get<PanelBody>(ellipse).panels == 33);
    SceneElement rectangle = PanelBody{.shape = PanelShape::Rectangle};
    n.set(rectangle, 33.4);
    CHECK(std::get<PanelBody>(rectangle).panels == 34);
    n.set(rectangle, std::nan(""));
    CHECK(sceneElementError(rectangle));
}

TEST_CASE("the solved panels let no fluid through, and the inside is outside the flow domain") {
    const FlowModel model = buildFlowModel(ellipseScene());
    REQUIRE(model.field.panelBodies.size() == 1);
    const SourcePanels& body = model.field.panelBodies.front();
    CHECK(body.element == 1);
    CHECK(normalFlowResidual(body, model.field) < 1e-13);
    CHECK_FALSE(model.field.inDomain(k_ellipse.position));
    CHECK(model.field.panelBodyAt(k_ellipse.position) == &body);
    CHECK(model.field.inDomain(k_ellipse.position + Vec2{0.0, 2.0}));
    CHECK(model.field.panelBodyAt({3.0, 2.0}) == nullptr);
    // Not solved on request (the parameter dialog only needs the warnings).
    CHECK(buildFlowModel(ellipseScene(), false).field.panelBodies.empty());
}

TEST_CASE("panel bodies feel each other and the other elements") {
    Scene scene{.elements = {UniformFlow{}, PanelBody{.position = {-1.0, 0.0}, .semiAxisA = 0.5, .semiAxisB = 0.5},
                             PanelBody{.shape = PanelShape::Rectangle, .position = {1.0, 0.3}, .semiAxisA = 0.4,
                                       .semiAxisB = 0.3, .panels = 30},
                             SourceSink{.strength = 2.0, .position = {0.0, -1.5}},
                             Vortex{.circulation = 1.5, .position = {0.0, 1.5}}}};
    const FlowModel model = buildFlowModel(scene);
    REQUIRE(model.field.panelBodies.size() == 2);
    for (const SourcePanels& body : model.field.panelBodies) CHECK(normalFlowResidual(body, model.field) < 1e-13);

    const Scene alone{.elements = {UniformFlow{}, scene.elements[1]}};
    const FlowModel single = buildFlowModel(alone);
    double difference = 0.0;
    for (std::size_t j = 0; j < single.field.panelBodies.front().size(); ++j) {
        difference = std::max(difference, std::abs(single.field.panelBodies.front().strengths[j] -
                                                   model.field.panelBodies.front().strengths[j]));
    }
    CHECK(difference > 0.05);
}

TEST_CASE("the complex potential of the panels is continuous around the body and its derivative is W") {
    const FlowModel model = buildFlowModel(ellipseScene());
    const SourcePanels& body = model.field.panelBodies.front();
    // On a circle around the body: no jumps (the logarithms are cut inside the body).
    const double radius = 1.1 * body.radius;
    constexpr int samples = 4000;
    Complex previous = body.complexPotential(body.centre + Vec2{radius, 0.0});
    double largestStep = 0.0;
    for (int k = 1; k <= samples; ++k) {
        const double t = 2.0 * pi * k / samples;
        const Complex f = body.complexPotential(body.centre + Vec2{std::cos(t), std::sin(t)} * radius);
        largestStep = std::max(largestStep, std::abs(f - previous));
        previous = f;
    }
    CHECK(largestStep < 0.01);
    // dF/dz = W, near and far.
    for (const Vec2 p : {Vec2{2.3, 0.4}, Vec2{-1.2, 1.1}, Vec2{0.4, -1.4}, Vec2{6.0, 5.0}}) {
        constexpr double h = 1e-5;
        const Complex derivative = (body.complexPotential(p + Vec2{h, 0.0}) - body.complexPotential(p - Vec2{h, 0.0})) / (2.0 * h);
        CHECK(std::abs(derivative - body.complexVelocity(p)) < 1e-7);
    }
}

TEST_CASE("the Laurent series continues the panel sum beyond two radii") {
    const FlowModel model = buildFlowModel(ellipseScene());
    const SourcePanels& body = model.field.panelBodies.front();
    for (int k = 0; k < 24; ++k) {
        const double t = 2.0 * pi * k / 24.0;
        const Vec2 p = body.centre + Vec2{std::cos(t), std::sin(t)} * (2.0 * body.radius * 1.001);
        Complex direct{};
        for (std::size_t j = 0; j < body.size(); ++j) {
            direct += body.vertexWeights[j] * std::log(toComplex(p - body.vertices[j]) / toComplex(p - body.centre));
        }
        CHECK(std::abs(body.complexVelocity(p) - direct / (2.0 * pi)) < 1e-14);
    }
}

TEST_CASE("the exact ellipse flow is the cylinder flow for a = b, and runs along the surface") {
    const Vec2 stream{1.1, 0.4};
    for (const Vec2 p : {Vec2{2.0, 0.5}, Vec2{-1.5, -1.0}, Vec2{0.2, 1.8}}) {
        // Cylinder of radius 1 at the origin: W = W∞ − conj(W∞)/z² + iΓ/2πz.
        const Complex z = toComplex(p);
        const Complex w = Complex{stream.x, -stream.y} - Complex{stream.x, stream.y} / (z * z) + Complex{0.0, 2.0} / (2.0 * pi * z);
        CHECK(std::abs(exactEllipseVelocity(1.0, 1.0, {}, 0.3, stream, 2.0, p) - w) < 1e-14);
    }
    // On the ellipse the velocity is tangential and its size is exactEllipseSurfaceSpeed.
    const double a = 1.7;
    const double b = 0.5;
    const double angle = 0.4;
    for (int k = 0; k < 12; ++k) {
        const double eta = 2.0 * pi * (k + 0.3) / 12.0;
        const Vec2 local{a * std::cos(eta), b * std::sin(eta)};
        const Vec2 tangent{-a * std::sin(eta), b * std::cos(eta)};
        const auto turn = [angle](Vec2 v) {
            return Vec2{std::cos(angle) * v.x - std::sin(angle) * v.y, std::sin(angle) * v.x + std::cos(angle) * v.y};
        };
        const Vec2 p = Vec2{0.5, -0.2} + turn(local * 1.0000001);
        const Vec2 v = velocityFromW(exactEllipseVelocity(a, b, {0.5, -0.2}, angle, stream, 1.3, p));
        CHECK(std::abs(cross(v, turn(tangent))) / (length(v) * length(tangent) + 1e-300) < 1e-5);
        CHECK(length(v) == doctest::Approx(exactEllipseSurfaceSpeed(a, b, angle, stream, 1.3, eta)).epsilon(1e-5));
    }
}

TEST_CASE("a panel body's surface comes from its control points") {
    const Scene scene{.elements = {UniformFlow{}, PanelBody{.panels = 40}}};
    const FlowModel model = buildFlowModel(scene);
    const std::vector<BodySurface> bodies = findBodySurfaces(model, k_view);
    REQUIRE(bodies.size() == 1);
    const BodySurface& body = bodies.front();
    CHECK(body.closed);
    REQUIRE(body.panelBody == std::optional<std::size_t>{0});
    CHECK(body.points.size() == 41);
    CHECK(body.points.front() == body.points.back());
    // It starts at the front stagnation point, and the rear one is on it too.
    CHECK(body.points.front().x < -0.9);
    REQUIRE(body.stagnationPoints.size() == 2);
    std::vector<Vec2> stagnation = body.stagnationPoints;
    std::ranges::sort(stagnation, {}, &Vec2::x);
    CHECK(stagnation[0].x == doctest::Approx(-std::cos(pi / 40.0)).epsilon(1e-6));
    CHECK(stagnation[1].x == doctest::Approx(std::cos(pi / 40.0)).epsilon(1e-6));
    CHECK(model.field.inDomain(stagnation[0]));
    CHECK(body.cp.front() == doctest::Approx(1.0));
    // Exact for a circle at the control points; no force.
    REQUIRE(body.cpExact.size() == body.cp.size());
    for (std::size_t k = 0; k < body.cp.size(); ++k) CHECK(body.cp[k] == doctest::Approx(body.cpExact[k]).epsilon(1e-12));
    const BodyForces forces = bodyForces(body, model);
    CHECK(forces.isolated);
    CHECK(std::abs(forces.drag) < 1e-12);
    CHECK(std::abs(forces.lift) < 1e-12);
    CHECK_FALSE(cylinderReference(body, model));
}

TEST_CASE("a panel body's surface carries the panel strengths, and the exact sheet where it is known") {
    SUBCASE("circle: lambda tends to -2 V cos(theta - alpha)") {
        const FlowModel model = buildFlowModel({.elements = {UniformFlow{.speed = 1.5}, PanelBody{.panels = 400}}});
        const std::vector<BodySurface> bodies = findBodySurfaces(model, k_view);
        REQUIRE(bodies.size() == 1);
        const BodySurface& body = bodies.front();
        REQUIRE(body.strength.size() == body.points.size());
        REQUIRE(body.strengthExact.size() == body.points.size());
        for (std::size_t k = 0; k < body.points.size(); ++k) {
            const SourcePanels& panels = model.field.panelBodies.front();
            CHECK(body.strength[k] == panels.strengths[(body.firstPanel + k) % panels.size()]);
            const double theta = std::atan2(body.points[k].y, body.points[k].x);
            CHECK(body.strengthExact[k] == doctest::Approx(-2.0 * 1.5 * std::cos(theta)).epsilon(1e-9));
            CHECK(std::abs(body.strength[k] - body.strengthExact[k]) < 0.02);
        }
    }
    SUBCASE("a vortex at a circle's centre leaves the sheet as it is; one off the centre of an ellipse does not") {
        const FlowModel centred =
            buildFlowModel({.elements = {UniformFlow{}, PanelBody{.panels = 64}, Vortex{.circulation = 3.0}}});
        CHECK_FALSE(findBodySurfaces(centred, k_view).front().strengthExact.empty());
        const FlowModel offCentre = buildFlowModel(
            {.elements = {UniformFlow{}, PanelBody{.semiAxisA = 1.5, .semiAxisB = 0.8, .panels = 64}, Vortex{.position = {0.2, 0.1}}}});
        const BodySurface body = findBodySurfaces(offCentre, k_view).front();
        CHECK(body.strengthExact.empty());
        CHECK_FALSE(body.cpExact.empty()); // the outside flow is still known
    }
    SUBCASE("traced bodies have none") {
        const FlowModel model = buildFlowModel({.elements = {UniformFlow{}, Doublet{}}});
        const BodySurface body = findBodySurfaces(model, k_view).front();
        CHECK(body.strength.empty());
        CHECK(body.strengthExact.empty());
    }
}

TEST_CASE("the panel strength view is saved with the scene") {
    Scene scene{.elements = {UniformFlow{}, PanelBody{}}};
    CHECK_FALSE(scene.plot.showPanelStrengths);
    scene.plot.showPanelStrengths = true;
    const SceneLoadResult loaded = sceneFromJson(sceneToJson(scene));
    REQUIRE(loaded.scene);
    CHECK(loaded.scene->plot.showPanelStrengths);
}

TEST_CASE("a stagnation point exactly at a control point is found") {
    // 8 panels: control points at 0° and 180°, where the tangential speed is 0 to rounding.
    const FlowModel model = buildFlowModel({.elements = {UniformFlow{}, PanelBody{.panels = 8}}});
    const std::vector<BodySurface> bodies = findBodySurfaces(model, k_view);
    REQUIRE(bodies.size() == 1);
    std::vector<Vec2> points = bodies.front().stagnationPoints;
    REQUIRE(points.size() == 2);
    std::ranges::sort(points, {}, &Vec2::x);
    CHECK(points[0].x == doctest::Approx(-std::cos(pi / 8.0)));
    CHECK(points[1].x == doctest::Approx(std::cos(pi / 8.0)));
    CHECK(findStagnationPoints(model.field, k_view, 1.0).size() == 2);
}

TEST_CASE("stagnation points on a panel body are found and followed") {
    const FlowModel model = buildFlowModel({.elements = {UniformFlow{}, PanelBody{.panels = 40}}});
    std::vector<Vec2> points = findStagnationPoints(model.field, k_view, 1.0);
    REQUIRE(points.size() == 2);
    std::ranges::sort(points, {}, &Vec2::x);
    CHECK(points[0].x == doctest::Approx(-1.0).epsilon(1e-2));
    CHECK(points[1].x == doctest::Approx(1.0).epsilon(1e-2));
    const std::optional<Vec2> followed = followStagnationPoint(model.field, points[0] + Vec2{0.0, 0.01}, 1.0, 1.0);
    REQUIRE(followed);
    CHECK(length(*followed - points[0]) < 1e-12);
}

TEST_CASE("a vortex inside a panel circle makes a lifting cylinder") {
    const FlowModel model = buildFlowModel({.elements = {UniformFlow{}, PanelBody{.panels = 128}, Vortex{.circulation = 2.0 * pi}}});
    const std::vector<BodySurface> bodies = findBodySurfaces(model, k_view);
    REQUIRE(bodies.size() == 1); // the vortex does not make a second, traced body
    const BodySurface& body = bodies.front();
    REQUIRE_FALSE(body.cpExact.empty());
    const BodyForces forces = bodyForces(body, model);
    CHECK(forces.isolated);
    CHECK(forces.circulation == doctest::Approx(2.0 * pi));
    CHECK(forces.lift == doctest::Approx(4.0 * pi).epsilon(1e-3));
    CHECK(std::abs(forces.drag) < 1e-10);
    // Stagnation points at 210° and 330°.
    REQUIRE(body.stagnationPoints.size() == 2);
    for (const Vec2 p : body.stagnationPoints) CHECK(p.y == doctest::Approx(-0.5).epsilon(1e-3));
}

TEST_CASE("only the circulation of what lies inside a panel body shows outside it") {
    const Scene plain{.elements = {UniformFlow{}, PanelBody{.panels = 160}}};
    Scene withSource = plain;
    withSource.elements.push_back(SourceSink{.strength = 1.5, .position = {0.2, 0.1}});
    withSource.elements.push_back(Doublet{.strength = 0.2, .position = {-0.3, 0.0}, .angle = 0.7});
    const FlowModel a = buildFlowModel(plain);
    const FlowModel b = buildFlowModel(withSource);
    double worst = 0.0;
    for (int k = 0; k < 36; ++k) {
        const double t = 2.0 * pi * k / 36.0;
        const Vec2 p{2.0 * std::cos(t), 2.0 * std::sin(t)};
        worst = std::max(worst, length(a.field.velocity(p) - b.field.velocity(p)));
    }
    CHECK(worst < 0.02);
    CHECK(contains(sceneWarnings(withSource, b), "lies inside source-panel body #2"));
}

TEST_CASE("panel body warnings: walls, overlaps") {
    const Scene walled{.elements = {UniformFlow{}, PanelBody{.position = {0.0, 1.5}}, SymmetryWall{}}};
    CHECK(contains(sceneWarnings(walled, buildFlowModel(walled)), "has no images in the walls"));
    const Scene overlapping{.elements = {UniformFlow{}, PanelBody{}, PanelBody{.position = {1.5, 0.0}}}};
    CHECK(contains(sceneWarnings(overlapping, buildFlowModel(overlapping)), "overlap"));
    const Scene apart{.elements = {UniformFlow{}, PanelBody{}, PanelBody{.position = {2.5, 0.0}}}};
    CHECK(sceneWarnings(apart, buildFlowModel(apart)).empty());
}

TEST_CASE("panel bodies are other bodies for the traced ones") {
    const FlowModel model = buildFlowModel(
        {.elements = {UniformFlow{}, Doublet{.strength = 0.25, .position = {-2.0, 0.0}}, PanelBody{.position = {1.5, 0.0}}}});
    const std::vector<BodySurface> bodies = findBodySurfaces(model, k_view);
    REQUIRE(bodies.size() == 2);
    CHECK(bodies[0].panelBody);
    CHECK_FALSE(bodies[1].panelBody);
    CHECK(bodies[0].cpExact.empty()); // not alone in the stream
    CHECK_FALSE(bodyForces(bodies[0], model).isolated);
    CHECK_FALSE(bodyForces(bodies[1], model).isolated);
    CHECK_FALSE(cylinderReference(bodies[1], model));
}

TEST_CASE("panel bodies are described, moved, saved and reopened") {
    SceneElement circle = PanelBody{.position = {1.0, 2.0}, .panels = 32};
    CHECK(describe(circle) == "circle a = 1 at (1, 2), N = 32");
    CHECK(describe(SceneElement{PanelBody{.shape = PanelShape::Rectangle, .semiAxisA = 1.0, .semiAxisB = 0.25,
                                          .angle = 90.0 * deg, .panels = 20}}) == "2 × 0.5 rectangle at (0, 0), θ = 90°, N = 20");
    CHECK(anchorPoint(circle) == Vec2{1.0, 2.0});
    CHECK(std::get<PanelBody>(movedTo(circle, {-1.0, 0.5})).position == Vec2{-1.0, 0.5});
    CHECK(kindOf(circle) == ElementKind::PanelBody);

    const Scene scene{.elements = {UniformFlow{}, k_ellipse,
                                   PanelBody{.shape = PanelShape::Rectangle, .position = {-2.0, 1.0}, .semiAxisA = 0.3,
                                             .semiAxisB = 0.2, .angle = 15.0 * deg, .panels = 24}}};
    const SceneLoadResult loaded = sceneFromJson(sceneToJson(scene));
    REQUIRE(loaded.scene);
    REQUIRE(loaded.scene->elements.size() == 3);
    for (std::size_t i = 1; i < 3; ++i) {
        const PanelBody& before = std::get<PanelBody>(scene.elements[i]);
        const PanelBody& after = std::get<PanelBody>(loaded.scene->elements[i]);
        CHECK(after.shape == before.shape);
        CHECK(after.panels == before.panels);
        CHECK(after.position == before.position);
        CHECK(after.semiAxisA == before.semiAxisA);
        CHECK(after.semiAxisB == before.semiAxisB);
        CHECK(after.angle == doctest::Approx(before.angle));
    }
    const std::string bad = R"({"format": "pfl-scene", "version": 1, "elements": [
        {"type": "source_panels", "shape": "star", "x": 0, "y": 0, "a": 1, "b": 1, "panels": 20}]})";
    CHECK(sceneFromJson(bad).error.find("unknown panel body shape") != std::string::npos);
    const std::string odd = R"({"format": "pfl-scene", "version": 1, "elements": [
        {"type": "source_panels", "shape": "rectangle", "x": 0, "y": 0, "a": 1, "b": 1, "panels": 21}]})";
    CHECK(sceneFromJson(odd).error.find("even number") != std::string::npos);
}
