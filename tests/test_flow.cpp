#include "core/flow_element.hpp"
#include "core/flow_field.hpp"
#include "core/vec2.hpp"

#include <doctest/doctest.h>

#include <cmath>
#include <numbers>

using namespace pfl;

namespace {

constexpr double pi = std::numbers::pi;
constexpr double deg = pi / 180.0;

Vec2 onCircle(Vec2 centre, double r, double theta) { return centre + Vec2{r * std::cos(theta), r * std::sin(theta)}; }

void checkVelocity(const FlowField& field, Vec2 p, Vec2 expected) {
    INFO("at (" << p.x << ", " << p.y << ")");
    const Vec2 v = field.velocity(p);
    CHECK(v.x == doctest::Approx(expected.x));
    CHECK(v.y == doctest::Approx(expected.y));
}

} // namespace

// Section 5 of the reference equations (V_inf = 1, a = 1).
TEST_CASE("periodic rows: the closed form is the limit of the series, dF/dz is W, and nothing overflows") {
    const Vec2 period = Vec2{-0.8, 1.7} * 1.5;
    const FlowElement elements[] = {SourceSink{.strength = 1.7, .position = {0.2, -0.3}},
                                    Vortex{.circulation = -2.3, .position = {0.2, -0.3}},
                                    Doublet{.strength = 0.9, .position = {0.2, -0.3}, .angle = 50.0 * deg}};
    for (const FlowElement& element : elements) {
        const PeriodicRow row{.element = element, .period = period};

        // Partial sums approach the limit as c / (N + 1/2), so two of them cancel c.
        const Vec2 p{1.1, 0.4};
        constexpr int n1 = 2000;
        constexpr int n2 = 4000;
        Complex total = complexVelocity(element, p);
        Complex first{};
        for (int n = 1; n <= n2; ++n) {
            total += complexVelocity(translated(element, period * n), p) + complexVelocity(translated(element, period * -n), p);
            if (n == n1) first = total;
        }
        const Complex limit = (total * (n2 + 0.5) - first * (n1 + 0.5)) / (n2 - n1 + 0.0);
        CHECK(std::abs(complexVelocity(row, p) - limit) < 1e-8);

        // Close to any copy the row is that element, plus the smooth field of the others.
        const Vec2 centre = Vec2{0.2, -0.3} + period * 3.0;
        CHECK(length(nearestCopy(row, centre + period * 0.4) - centre) < 1e-12);
        const Vec2 near = centre + Vec2{1e-4, 1e-4};
        CHECK(std::abs(complexVelocity(row, near) - complexVelocity(translated(element, period * 3.0), near)) < 5.0);

        // Along the row direction it repeats; at right angles it stays finite to any distance.
        const Vec2 along = perp(period) * 2.0;
        for (const double distance : {0.3, 1.0, 5.0, 19.0, 21.0, 40.0, 1e4}) {
            for (const double side : {1.0, -1.0}) {
                const Vec2 q = Vec2{0.2, -0.3} + along * (side * distance) + period * 0.3;
                const double h = 1e-5;
                const Complex dFdx = (complexPotential(row, q + Vec2{h, 0.0}) - complexPotential(row, q - Vec2{h, 0.0})) / (2.0 * h);
                const Complex dFdy = (complexPotential(row, q + Vec2{0.0, h}) - complexPotential(row, q - Vec2{0.0, h})) / (2.0 * h);
                const Complex w = complexVelocity(row, q);
                INFO("distance " << distance << ", side " << side);
                REQUIRE(std::isfinite(std::abs(w)));
                // At huge distances the potential is too large for the difference to resolve W.
                const double tolerance = distance > 100.0 ? 1e-3 : 1e-7;
                CHECK(std::abs(dFdx - w) < tolerance);
                CHECK(std::abs(dFdy - Complex{0.0, 1.0} * w) < tolerance);
                CHECK(std::abs(complexVelocity(row, q + period * 3.0) - w) < 1e-9 * (1.0 + std::abs(w)));
            }
        }
    }
}

TEST_CASE("reference unit-test values") {
    checkVelocity({.elements = {UniformFlow{.speed = 1.0, .angle = 45.0 * deg}}}, {3.0, -2.0}, {0.707107, 0.707107});

    checkVelocity({.elements = {SourceSink{.strength = 2.0 * pi}}}, {1.0, 0.0}, {1.0, 0.0});

    const FlowField vortex{.elements = {Vortex{.circulation = 2.0 * pi}}};
    checkVelocity(vortex, {1.0, 0.0}, {0.0, -1.0}); // clockwise
    checkVelocity(vortex, {0.0, 1.0}, {1.0, 0.0});

    const FlowField doublet{.elements = {Doublet{.strength = 1.0, .angle = 0.0}}};
    checkVelocity(doublet, {1.0, 0.0}, {-1.0, 0.0});
    checkVelocity(doublet, {0.0, 1.0}, {1.0, 0.0});

    FlowField cylinder{.elements = {UniformFlow{.speed = 1.0}, Doublet{.strength = 1.0}}};
    checkVelocity(cylinder, {0.0, 1.0}, {2.0, 0.0});

    cylinder.elements.push_back(Vortex{.circulation = 2.0 * pi});
    checkVelocity(cylinder, {0.0, 1.0}, {3.0, 0.0});
    checkVelocity(cylinder, onCircle({}, 1.0, 210.0 * deg), {0.0, 0.0});
    checkVelocity(cylinder, onCircle({}, 1.0, 330.0 * deg), {0.0, 0.0});

    checkVelocity({.elements = {SectorFlow{.strength = 1.0, .wall1 = 0.0, .wall2 = 90.0 * deg}}}, {1.0, 1.0},
                  {2.0, -2.0});
}

TEST_CASE("sector walls are streamlines and the field is masked outside the wedge") {
    const Vec2 vertex{0.5, -1.0};
    const SectorFlow sector{.strength = 1.3, .vertex = vertex, .wall1 = 30.0 * deg, .wall2 = 120.0 * deg};
    const FlowField field{.elements = {sector}};

    for (const double wall : {sector.wall1, sector.wall2}) {
        const Vec2 normal = perp(Vec2{std::cos(wall), std::sin(wall)});
        for (const double r : {0.25, 1.0, 2.5}) {
            const Vec2 p = onCircle(vertex, r, wall);
            INFO("wall " << wall / deg << " deg, r = " << r);
            CHECK(field.inDomain(p));
            CHECK(dot(field.velocity(p), normal) == doctest::Approx(0.0));
            CHECK(length(field.velocity(p)) > 0.0);
        }
    }
    CHECK(field.inDomain(onCircle(vertex, 1.0, 75.0 * deg)));
    CHECK_FALSE(field.inDomain(onCircle(vertex, 1.0, 200.0 * deg)));
    CHECK_FALSE(field.inDomain(onCircle(vertex, 1.0, 10.0 * deg)));

    // A 360-degree sector (flow around the edge of a plate) admits the whole plane.
    const FlowField plate{.elements = {SectorFlow{.wall1 = 0.0, .wall2 = 2.0 * pi}}};
    CHECK(plate.inDomain({-1.0, -1e-6}));
    CHECK(plate.inDomain({1.0, -1e-12}));
}

TEST_CASE("complex velocity is the derivative of the complex potential") {
    const Vec2 c{0.3, -0.2};
    const FlowElement elements[] = {
        UniformFlow{.speed = 1.7, .angle = 25.0 * deg},
        SourceSink{.strength = -2.4, .position = c},
        Vortex{.circulation = 3.1, .position = c},
        Doublet{.strength = 0.8, .position = c, .angle = 40.0 * deg},
        SectorFlow{.strength = 1.2, .vertex = c, .wall1 = 30.0 * deg, .wall2 = 160.0 * deg},  // n ~ 1.38
        SectorFlow{.strength = 0.9, .vertex = c, .wall1 = -45.0 * deg, .wall2 = 225.0 * deg}, // n = 2/3
    };
    // Inside every wedge and away from the log branch cut along -x from c.
    const Vec2 points[] = {{1.1, 0.9}, {0.5, 1.4}, {-0.4, 0.8}};
    constexpr double h = 1e-6;

    for (const FlowElement& e : elements) {
        for (const Vec2 p : points) {
            INFO("element " << e.index() << " at (" << p.x << ", " << p.y << ")");
            const Complex w = complexVelocity(e, p);
            const Complex dFdx = (complexPotential(e, p + Vec2{h, 0.0}) - complexPotential(e, p - Vec2{h, 0.0})) / (2.0 * h);
            const Complex dFdy = (complexPotential(e, p + Vec2{0.0, h}) - complexPotential(e, p - Vec2{0.0, h})) / (2.0 * h);
            // dF/dz = dF/dx = -i dF/dy  (i.e. u = phi_x = psi_y, v = phi_y = -psi_x).
            CHECK(std::abs(dFdx - w) < 1e-6 * (1.0 + std::abs(w)));
            CHECK(std::abs(dFdy - Complex{0.0, 1.0} * w) < 1e-6 * (1.0 + std::abs(w)));
        }
    }
}

TEST_CASE("rotated doublet in uniform flow at the same angle is a cylinder") {
    const double alpha = 30.0 * deg;
    const double speed = 2.0;
    const double radius = 0.7;
    const Vec2 centre{1.0, -0.5};
    const FlowField field{.elements = {UniformFlow{.speed = speed, .angle = alpha},
                                       Doublet{.strength = speed * radius * radius, .position = centre, .angle = alpha}}};

    for (double theta = 0.0; theta < 2.0 * pi; theta += 17.0 * deg) {
        const Vec2 p = onCircle(centre, radius, theta);
        const Vec2 radial{std::cos(theta), std::sin(theta)};
        INFO("theta = " << theta / deg << " deg");
        CHECK(dot(field.velocity(p), radial) == doctest::Approx(0.0)); // V_r = 0 on the surface
        CHECK(dot(field.velocity(p), perp(radial)) == doctest::Approx(-2.0 * speed * std::sin(theta - alpha)));
    }
    checkVelocity(field, onCircle(centre, radius, alpha), {0.0, 0.0});
    checkVelocity(field, onCircle(centre, radius, alpha + pi), {0.0, 0.0});
}

TEST_CASE("singular terms are skipped at their centres") {
    const FlowField field{.elements = {UniformFlow{.speed = 1.0}, SourceSink{.strength = 1.0}, Vortex{.circulation = 1.0},
                                       Doublet{.strength = 1.0}}};
    checkVelocity(field, {0.0, 0.0}, {1.0, 0.0});
    CHECK(field.nearSingularity({0.01, 0.0}, 0.05));
    CHECK_FALSE(field.nearSingularity({0.1, 0.0}, 0.05));
    CHECK_FALSE(FlowField{.elements = {UniformFlow{}}}.nearSingularity({0.0, 0.0}, 1.0));

    // A sector vertex is a stagnation point for n > 1 and singular for n < 1.
    CHECK_FALSE(singularPoint(SectorFlow{.wall1 = 0.0, .wall2 = 90.0 * deg}).has_value());
    CHECK(singularPoint(SectorFlow{.wall1 = 0.0, .wall2 = 270.0 * deg}).has_value());
    checkVelocity({.elements = {SectorFlow{.wall1 = 0.0, .wall2 = 90.0 * deg}}}, {0.0, 0.0}, {0.0, 0.0});
}

TEST_CASE("invalid parameters are reported") {
    CHECK_FALSE(validationError(SectorFlow{.wall1 = 0.0, .wall2 = 2.0 * pi}).has_value());
    CHECK(validationError(SectorFlow{.wall1 = 30.0 * deg, .wall2 = 30.0 * deg}).has_value());
    CHECK(validationError(SectorFlow{.wall1 = 90.0 * deg, .wall2 = 0.0}).has_value());
    CHECK(validationError(SectorFlow{.wall1 = 0.0, .wall2 = 400.0 * deg}).has_value());
    CHECK(validationError(UniformFlow{.speed = -1.0}).has_value());
    CHECK(validationError(Vortex{.circulation = 1.0, .position = {std::nan(""), 0.0}}).has_value());
    CHECK_FALSE(validationError(SourceSink{.strength = -3.0, .position = {1.0, 2.0}}).has_value());
}
