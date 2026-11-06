#include "core/scene.hpp"
#include "core/surface_pressure.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

using namespace pfl;

namespace {

constexpr double pi = std::numbers::pi;
constexpr double deg = pi / 180.0;
const Rect k_view{{-4.0, -2.5}, {4.0, 2.5}};

std::vector<BodySurface> surfacesOf(const Scene& scene) { return findBodySurfaces(buildFlowModel(scene), k_view); }

double minCp(const BodySurface& body) { return *std::ranges::min_element(body.cp); }
double maxCp(const BodySurface& body) { return *std::ranges::max_element(body.cp); }

} // namespace

TEST_CASE("stagnation points of the cylinder are at (+-a, 0)") {
    const Scene scene{.elements = {UniformFlow{}, Doublet{.strength = 4.0}}}; // a = 2
    const FlowModel model = buildFlowModel(scene);
    std::vector<Vec2> points = findStagnationPoints(model.field, {{-5.0, -5.0}, {5.0, 5.0}}, 1.0);
    REQUIRE(points.size() == 2);
    std::ranges::sort(points, {}, &Vec2::x);
    CHECK(points[0].x == doctest::Approx(-2.0));
    CHECK(points[1].x == doctest::Approx(2.0));
    CHECK(points[0].y == doctest::Approx(0.0).epsilon(1e-9));
}

TEST_CASE("a body in a channel is found, and is not treated as isolated") {
    // A cylinder of radius about 1 between walls at y = +-3: blocked a little, so faster at its sides.
    const Scene scene{.elements = {UniformFlow{}, Doublet{}, Vortex{.circulation = 1.0}, SymmetryWall{.point = {0.0, -3.0}},
                                   SymmetryWall{.point = {0.0, 3.0}}}};
    const FlowModel model = buildFlowModel(scene);
    const std::vector<BodySurface> bodies = findBodySurfaces(model, k_view);
    const auto body = std::ranges::find_if(bodies, [](const BodySurface& b) { return length(b.centroid) < 0.5; });
    REQUIRE(body != bodies.end());
    CHECK(body->closed);
    const BodyForces forces = bodyForces(*body, model);
    CHECK_FALSE(forces.isolated);
    CHECK(forces.circulation == doctest::Approx(1.0));
    CHECK_FALSE(cylinderReference(*body, model).has_value());
}

TEST_CASE("cylinder: surface Cp matches 1 - 4 sin^2(theta), no drag, no lift") {
    const Scene scene{.elements = {UniformFlow{}, Doublet{}}};
    const FlowModel model = buildFlowModel(scene);
    const std::vector<BodySurface> bodies = findBodySurfaces(model, k_view);
    REQUIRE(bodies.size() == 1);
    const BodySurface& body = bodies.front();
    CHECK(body.closed);
    CHECK(body.stagnationPoints.size() == 2);
    CHECK(body.points.front() == body.points.back());
    CHECK(body.points.front().x == doctest::Approx(-1.0)); // starts at the front stagnation point
    CHECK(body.points[body.points.size() / 4].y > 0.0);    // clockwise: upper side first
    CHECK(body.arcLength.back() == doctest::Approx(2.0 * pi).epsilon(1e-4));
    CHECK(maxCp(body) == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(minCp(body) == doctest::Approx(-3.0).epsilon(1e-4));
    CHECK(body.centroid.x == doctest::Approx(0.0).epsilon(1e-5));

    const std::optional<CylinderReference> exact = cylinderReference(body, model);
    REQUIRE(exact.has_value());
    CHECK(exact->radius == doctest::Approx(1.0));
    double error = 0.0;
    for (std::size_t k = 0; k < body.points.size(); ++k) {
        const Vec2 r = body.points[k] - exact->centre;
        error = std::max(error, std::abs(body.cp[k] - exact->cp(std::atan2(r.y, r.x))));
        CHECK(std::abs(length(r) - 1.0) < 1e-6);
    }
    CHECK(error < 1e-5);

    const BodyForces forces = bodyForces(body, model);
    CHECK(forces.isolated);
    CHECK(forces.chord == doctest::Approx(2.0).epsilon(1e-6));
    CHECK(std::abs(forces.drag) < 1e-5);
    CHECK(std::abs(forces.lift) < 1e-5);
}

TEST_CASE("lifting cylinder at an angle: Kutta-Joukowski lift, zero drag") {
    const double alpha = 20.0 * deg;
    const double gamma = 2.0 * pi;
    const Scene scene{.elements = {UniformFlow{.speed = 2.0, .angle = alpha},
                                   Doublet{.strength = 2.0, .position = {0.5, -0.3}, .angle = alpha},
                                   Vortex{.circulation = gamma, .position = {0.5, -0.3}}}};
    const FlowModel model = buildFlowModel(scene);
    const std::vector<BodySurface> bodies = findBodySurfaces(model, k_view);
    REQUIRE(bodies.size() == 1);
    const BodySurface& body = bodies.front();
    REQUIRE(body.closed);

    const std::optional<CylinderReference> exact = cylinderReference(body, model);
    REQUIRE(exact.has_value());
    CHECK(exact->circulation == doctest::Approx(gamma));
    for (std::size_t k = 0; k < body.points.size(); k += 50) {
        const Vec2 r = body.points[k] - exact->centre;
        CHECK(body.cp[k] == doctest::Approx(exact->cp(std::atan2(r.y, r.x))).epsilon(1e-5));
    }

    const BodyForces forces = bodyForces(body, model);
    CHECK(forces.isolated);
    CHECK(forces.circulation == doctest::Approx(gamma));
    CHECK(std::abs(forces.drag) < 1e-4);
    CHECK(forces.lift == doctest::Approx(2.0 * gamma / 2.0).epsilon(1e-4)); // 2 Gamma / V_inf
}

TEST_CASE("strong circulation: the dividing streamline loop carries the Kutta-Joukowski lift") {
    // Gamma > 4 pi a V_inf: one stagnation point, off the cylinder.
    const double gamma = 6.0 * pi;
    const Scene scene{.elements = {UniformFlow{}, Doublet{}, Vortex{.circulation = gamma}}};
    const FlowModel model = buildFlowModel(scene);
    const std::vector<BodySurface> bodies = findBodySurfaces(model, k_view);
    const auto closed = std::ranges::find_if(bodies, &BodySurface::closed);
    REQUIRE(closed != bodies.end());
    CHECK(closed->stagnationPoints.size() == 1);
    CHECK_FALSE(cylinderReference(*closed, model).has_value());
    const BodyForces forces = bodyForces(*closed, model);
    CHECK(forces.isolated);
    CHECK(std::abs(forces.drag) < 1e-3);
    CHECK(forces.lift == doctest::Approx(2.0 * gamma).epsilon(1e-3));
}

TEST_CASE("Rankine oval is closed, the half-body is open") {
    const Scene oval{.elements = {UniformFlow{}, SourceSink{.strength = 2.0 * pi, .position = {-1.0, 0.0}},
                                  SourceSink{.strength = -2.0 * pi, .position = {1.0, 0.0}}}};
    std::vector<BodySurface> bodies = surfacesOf(oval);
    REQUIRE(bodies.size() == 1);
    CHECK(bodies[0].closed);
    // Stagnation points at x = +-sqrt(b^2 + Lambda b / (pi V)) = +-sqrt(3) (Anderson eq. 3.82).
    CHECK(std::abs(bodies[0].stagnationPoints[0].x) == doctest::Approx(std::sqrt(3.0)));
    const BodyForces forces = bodyForces(bodies[0], buildFlowModel(oval));
    CHECK(std::abs(forces.drag) < 1e-4);
    CHECK(std::abs(forces.lift) < 1e-4);
    CHECK(forces.chord == doctest::Approx(2.0 * std::sqrt(3.0)).epsilon(1e-6));

    const Scene half{.elements = {UniformFlow{}, SourceSink{.strength = 2.0 * pi, .position = {-1.0, 0.0}}}};
    bodies = surfacesOf(half);
    REQUIRE(bodies.size() == 1);
    const BodySurface& body = bodies[0];
    CHECK_FALSE(body.closed);
    // Stagnation point Lambda / (2 pi V) = 1 upstream of the source; it is where s = 0.
    CHECK(body.stagnationPoints[0].x == doctest::Approx(-2.0));
    const auto front = std::ranges::find(body.arcLength, 0.0);
    REQUIRE(front != body.arcLength.end());
    CHECK(body.cp[static_cast<std::size_t>(front - body.arcLength.begin())] == doctest::Approx(1.0).epsilon(1e-6));
    CHECK(body.points.back().y > 0.0); // positive s on the upper branch
    // Far downstream the half-body is Lambda / V = 2 pi wide and Cp tends to 0.
    CHECK(std::abs(body.points.back().y) < pi);
    CHECK(std::abs(body.cp.back()) < 0.2);
}

TEST_CASE("open surfaces of a source and a sink both have the body on the right") {
    for (const double strength : {2.0 * pi, -2.0 * pi}) {
        INFO("Lambda = " << strength);
        const Vec2 centre{0.5, 0.2};
        const Scene scene{.elements = {UniformFlow{.speed = 1.0, .angle = 30.0 * deg},
                                       SourceSink{.strength = strength, .position = centre}}};
        const std::vector<BodySurface> bodies = surfacesOf(scene);
        REQUIRE(bodies.size() == 1);
        const BodySurface& body = bodies[0];
        REQUIRE_FALSE(body.closed);
        for (std::size_t k = 0; k + 1 < body.points.size(); k += 200) {
            const Vec2 d = body.points[k + 1] - body.points[k];
            // The singularity is inside the body, i.e. on the right of the direction of travel.
            if (length(body.points[k] - centre) < 2.0) CHECK(cross(d, centre - body.points[k]) < 0.0);
        }
    }
}

TEST_CASE("a Rankine oval at incidence is not closed") {
    // psi is odd about the oval's centre, so the stagnation points lie on psi and -psi: the dividing
    // streamlines miss each other by far more than the closure tolerance.
    const Scene scene{.elements = {UniformFlow{.speed = 1.0, .angle = 10.0 * deg},
                                   SourceSink{.strength = 2.0 * pi, .position = {-1.0, 0.0}},
                                   SourceSink{.strength = -2.0 * pi, .position = {1.0, 0.0}}}};
    CHECK_FALSE(std::ranges::any_of(surfacesOf(scene), &BodySurface::closed));
}

TEST_CASE("no free stream, no surfaces") {
    const Scene scene{.elements = {SourceSink{.strength = 1.0}, Vortex{.position = {1.0, 0.0}}}};
    CHECK(surfacesOf(scene).empty());
}

TEST_CASE("a body near a wall is not isolated") {
    const Scene scene{.elements = {UniformFlow{}, Doublet{.position = {0.0, 2.0}}, SymmetryWall{}}};
    const FlowModel model = buildFlowModel(scene);
    const std::vector<BodySurface> bodies = findBodySurfaces(model, {{-4.0, -1.0}, {4.0, 4.0}});
    const auto upper = std::ranges::find_if(bodies, [](const BodySurface& b) { return b.closed && b.centroid.y > 0.0; });
    REQUIRE(upper != bodies.end());
    CHECK_FALSE(bodyForces(*upper, model).isolated);
    CHECK_FALSE(cylinderReference(*upper, model).has_value());
}


TEST_CASE("two cylinders: each nearly closed, the gap reported and small") {
    const Scene scene{.elements = {UniformFlow{}, Doublet{.strength = 0.25, .position = {-2.0, 0.0}},
                                   Doublet{.strength = 0.25, .position = {2.0, 0.5}}}};
    const FlowModel model = buildFlowModel(scene);
    const std::vector<BodySurface> bodies = findBodySurfaces(model, {{-5.0, -3.0}, {5.0, 3.0}});
    REQUIRE(bodies.size() == 2);
    for (const BodySurface& body : bodies) {
        CHECK(body.closed);
        CHECK(body.closureGap > 0.0);
        CHECK(body.closureGap < 0.05 * length(body.stagnationPoints[1] - body.stagnationPoints[0]));
        // Close to a cylinder of radius 0.5 (each sees the other's weak perturbation).
        const double chord = bodyForces(body, model).chord;
        CHECK(chord == doctest::Approx(1.0).epsilon(0.05));
        CHECK_FALSE(bodyForces(body, model).isolated);
    }
}

