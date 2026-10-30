#include "core/parameters.hpp"
#include "core/scene.hpp"
#include "core/sweep.hpp"
#include "io/data_export.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>

using namespace pfl;

namespace {

constexpr double pi = std::numbers::pi;
const Rect k_view{{-4.0, -2.5}, {4.0, 2.5}};

// Uniform flow V = 1 past a cylinder of radius a = 1, with a vortex of circulation Gamma at its centre.
Scene liftingCylinder(double gamma) {
    return Scene{.elements = {UniformFlow{}, Doublet{.strength = 1.0}, Vortex{.circulation = gamma}}};
}

std::size_t countLines(const std::string& text, const std::string& needle) {
    std::size_t n = 0;
    for (std::size_t at = text.find(needle); at != std::string::npos; at = text.find(needle, at + 1)) ++n;
    return n;
}

} // namespace

TEST_CASE("every element kind lists its parameters, and they read and write the element") {
    for (const ElementKind kind : k_elementKinds) {
        SceneElement element = defaultElement(kind, {0.25, -0.5});
        const std::span<const ElementParameter> parameters = elementParameters(kind);
        REQUIRE_FALSE(parameters.empty());
        for (const ElementParameter& p : parameters) {
            const double before = p.get(element);
            p.set(element, before + 0.125);
            CHECK(p.get(element) == doctest::Approx(before + 0.125));
            p.set(element, before);
        }
        if (const std::optional<Vec2> anchor = anchorPoint(element)) {
            const auto x = std::ranges::find(parameters, ParameterKind::X, &ElementParameter::kind);
            const auto y = std::ranges::find(parameters, ParameterKind::Y, &ElementParameter::kind);
            REQUIRE(x != parameters.end());
            REQUIRE(y != parameters.end());
            CHECK(x->get(element) == anchor->x);
            CHECK(y->get(element) == anchor->y);
        }
    }
}

TEST_CASE("withParameter rejects invalid values and references") {
    const Scene scene{.elements = {UniformFlow{}, SectorFlow{}}};
    CHECK(withParameter(scene, {0, 0}, 2.0).has_value());   // U = 2
    CHECK_FALSE(withParameter(scene, {0, 0}, -1.0));        // negative speed
    CHECK_FALSE(withParameter(scene, {1, 4}, -0.1));        // theta2 < theta1
    CHECK_FALSE(withParameter(scene, {2, 0}, 1.0));         // no element #3
    CHECK_FALSE(withParameter(scene, {0, 7}, 1.0));         // no such parameter
    CHECK(parameterName(scene, {1, 3}) == "#2 Sector (wedge) flow: θ₁");
}

TEST_CASE("sweep values include both ends") {
    const std::vector<double> values = sweepValues(-1.0, 3.0, 5);
    REQUIRE(values.size() == 5);
    CHECK(values.front() == -1.0);
    CHECK(values[2] == doctest::Approx(1.0));
    CHECK(values.back() == 3.0);
    CHECK(sweepValues(2.0, 4.0, 1) == std::vector<double>{2.0});
}

TEST_CASE("sweeping the circulation of a lifting cylinder: Kutta-Joukowski lift, stagnation points move down") {
    const Scene scene = liftingCylinder(0.0);
    const ParameterRef gamma{2, 0};
    for (const double value : sweepValues(0.0, 3.0 * pi, 4)) {
        CAPTURE(value);
        const std::optional<SweepSample> sample = sampleSweep(scene, gamma, value, k_view);
        REQUIRE(sample);
        CHECK(sample->value == value);
        // On the surface at sin(theta) = -Gamma / (4 pi a V).
        REQUIRE(sample->stagnationPoints.size() == 2);
        for (const Vec2 p : sample->stagnationPoints) {
            CHECK(p.y == doctest::Approx(-value / (4.0 * pi)).epsilon(1e-6));
            CHECK(length(p) == doctest::Approx(1.0).epsilon(1e-6));
        }
        REQUIRE(sample->bodies.size() == 1);
        const SweepBody& body = sample->bodies.front();
        CHECK(body.closed);
        CHECK(body.chord == doctest::Approx(2.0).epsilon(1e-4));
        CHECK(body.circulation == doctest::Approx(value));
        CHECK(body.exactLift == doctest::Approx(2.0 * value));
        CHECK(body.lift == doctest::Approx(2.0 * value).epsilon(1e-4));
        CHECK(std::abs(body.drag) < 1e-4);
        CHECK(body.cpMin == doctest::Approx(1.0 - std::pow(2.0 + value / (2.0 * pi), 2.0)).epsilon(1e-3));
    }
}

TEST_CASE("an invalid sweep value gives no sample") {
    const Scene scene{.elements = {UniformFlow{}}};
    CHECK_FALSE(sampleSweep(scene, {0, 0}, -1.0, k_view));
}

TEST_CASE("sweep CSV: stagnation and body rows for every step") {
    const Scene scene = liftingCylinder(0.0);
    std::vector<SweepSample> samples;
    for (const double value : sweepValues(0.0, pi, 3)) samples.push_back(*sampleSweep(scene, {2, 0}, value, k_view));
    const std::string csv = sweepCsv(samples);
    CHECK(csv.starts_with("step,value,kind,index,x,y,closed,chord,circulation,lift_q,drag_q,cl,cd,cp_min,lift_q_exact\n"));
    CHECK(countLines(csv, ",stagnation,") == 6);
    CHECK(countLines(csv, ",body,") == 3);
    CHECK(countLines(csv, "\n") == 1 + 6 + 3);
    CHECK(csv.find("3,3.14159") != std::string::npos); // step 3 at Gamma = pi
}
