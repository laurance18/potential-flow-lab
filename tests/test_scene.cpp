#include "core/image_system.hpp"
#include "core/scene.hpp"
#include "core/streamline.hpp"
#include "io/data_export.hpp"
#include "io/examples.hpp"
#include "io/scene_io.hpp"

#include <doctest/doctest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <numbers>
#include <sstream>
#include <string>
#include <vector>

using namespace pfl;

namespace {

constexpr double pi = std::numbers::pi;
constexpr double deg = pi / 180.0;

// Largest |v . n| relative to the largest |v| along the wall within a 10 x 10 box.
double wallResidual(const FlowModel& model, const SymmetryWall& wall) {
    const WallCheck check = checkWall(model.field, wall, {{-5.0, -5.0}, {5.0, 5.0}}, 200);
    REQUIRE(check.samples > 100);
    REQUIRE(check.maxSpeed > 0.0);
    return check.maxNormalSpeed / check.maxSpeed;
}

} // namespace

TEST_CASE("images across a wall: sources keep, vortices flip, doublets reflect their axis") {
    const SymmetryWall wall{.point = {0.0, -1.0}, .angle = 0.0}; // the line y = -1
    CHECK(reflect({2.0, 1.0}, wall).x == doctest::Approx(2.0));
    CHECK(reflect({2.0, 1.0}, wall).y == doctest::Approx(-3.0));

    const auto source = std::get<SourceSink>(*mirror(SourceSink{.strength = 3.0, .position = {1.0, 0.0}}, wall));
    CHECK(source.strength == 3.0);
    CHECK(source.position.y == doctest::Approx(-2.0));
    const auto vortex = std::get<Vortex>(*mirror(Vortex{.circulation = 2.0, .position = {1.0, 0.0}}, wall));
    CHECK(vortex.circulation == -2.0);
    const auto doublet = std::get<Doublet>(*mirror(Doublet{.strength = 1.0, .angle = 30.0 * deg}, wall));
    CHECK(doublet.angle == doctest::Approx(-30.0 * deg));
    CHECK_FALSE(mirror(UniformFlow{}, wall).has_value());
}

TEST_CASE("a symmetry wall has zero normal velocity") {
    const SymmetryWall wall{.point = {0.3, -0.8}, .angle = 25.0 * deg};
    const Scene scene{.elements = {UniformFlow{.speed = 1.0, .angle = 10.0 * deg}, SourceSink{.strength = 4.0, .position = {1.0, 1.0}},
                                   Vortex{.circulation = -3.0, .position = {-1.5, 0.5}},
                                   Doublet{.strength = 0.7, .position = {0.5, 2.0}, .angle = 70.0 * deg}, wall}};
    const FlowModel model = buildFlowModel(scene);
    CHECK(model.imagesExact);
    CHECK(model.images.size() == 3);
    CHECK(wallResidual(model, wall) < 1e-12);
    // The uniform flow keeps only its component along the wall.
    CHECK(model.freeStreamSpeed == doctest::Approx(std::cos(15.0 * deg)));
}

TEST_CASE("two perpendicular walls close the image system exactly") {
    const SymmetryWall floor{.point = {0.0, 0.0}, .angle = 0.0};
    const SymmetryWall side{.point = {0.0, 0.0}, .angle = 90.0 * deg};
    const Scene scene{.elements = {Vortex{.circulation = 1.0, .position = {1.0, 2.0}}, floor, side}};
    const FlowModel model = buildFlowModel(scene);
    CHECK(model.imagesExact);
    CHECK(model.images.size() == 3);
    CHECK(wallResidual(model, floor) < 1e-12);
    CHECK(wallResidual(model, side) < 1e-12);
}

TEST_CASE("three parallel walls truncate the image series and say so") {
    const Scene scene{.elements = {SourceSink{.strength = 1.0, .position = {0.0, 0.2}}, SymmetryWall{.point = {0.0, -1.0}},
                                   SymmetryWall{.point = {0.0, 1.0}}, SymmetryWall{.point = {0.0, 2.5}}}};
    const FlowModel model = buildFlowModel(scene);
    CHECK_FALSE(model.imagesExact);
    CHECK(model.field.rows.empty());
    CHECK_FALSE(sceneWarnings(scene, model).empty());
}

TEST_CASE("a channel of two parallel walls is exact for every element") {
    for (const double flip : {0.0, 180.0}) { // the second wall may run the other way
        const SymmetryWall lower{.point = {0.3, -0.8}, .angle = 25.0 * deg};
        const SymmetryWall upper{.point = lower.point + lower.normal() * 2.0, .angle = (25.0 + flip) * deg};
        const Scene scene{.elements = {UniformFlow{.speed = 1.0, .angle = 10.0 * deg},
                                       SourceSink{.strength = 4.0, .position = {1.0, 0.0}},
                                       Vortex{.circulation = -3.0, .position = {-1.5, -0.5}},
                                       Doublet{.strength = 0.7, .position = {0.5, 0.9}, .angle = 70.0 * deg}, lower, upper}};
        const FlowModel model = buildFlowModel(scene);
        CHECK(model.imagesExact);
        CHECK(model.images.empty());
        CHECK(model.field.rows.size() == 6);
        CHECK(model.field.elements.size() == 1); // the uniform flow; the others are in their rows
        CHECK(wallResidual(model, lower) < 1e-12);
        CHECK(wallResidual(model, upper) < 1e-12);
        CHECK(model.freeStreamSpeed == doctest::Approx(std::cos(15.0 * deg)));
        CHECK(sceneWarnings(scene, model).size() == 1); // only that the stream is cut to its part along the walls
    }
}

TEST_CASE("channel: half of a source's flow leaves each way") {
    const Scene scene{.elements = {SourceSink{.strength = 3.0, .position = {0.0, 0.2}}, SymmetryWall{.point = {0.0, -1.0}},
                                   SymmetryWall{.point = {0.0, 1.0}}}};
    const FlowField field = buildFlowModel(scene).field;
    REQUIRE_FALSE(field.empty());
    CHECK(field.velocity({30.0, 0.5}).x == doctest::Approx(3.0 / 4.0));
    CHECK(field.velocity({-30.0, -0.5}).x == doctest::Approx(-3.0 / 4.0));
    CHECK(std::abs(field.velocity({30.0, 0.5}).y) < 1e-12);
}

TEST_CASE("channel: elements on a wall behave as they do at a single wall") {
    const SymmetryWall lower{.point = {0.0, -1.0}};
    const SymmetryWall upper{.point = {0.0, 1.0}};
    for (const double y : {-1.0, 1.0}) { // on the second wall the image is one period away
        const FlowModel source = buildFlowModel({.elements = {SourceSink{.strength = 3.0, .position = {0.0, y}}, lower, upper}});
        CHECK(source.field.rows.size() == 1); // its own image, not a second source
        CHECK(source.field.velocity({30.0, 0.0}).x == doctest::Approx(3.0 / 8.0));
        const FlowModel vortex = buildFlowModel({.elements = {Vortex{.circulation = 3.0, .position = {0.0, y}}, lower, upper}});
        CHECK(length(vortex.field.velocity({0.5, 0.3})) == doctest::Approx(0.0)); // cancelled by its image
    }
}

TEST_CASE("channel: an element and its own image in the scene are one pair, not two") {
    const SymmetryWall lower{.point = {0.0, -1.0}};
    const SymmetryWall upper{.point = {0.0, 1.0}};
    const SourceSink a{.strength = 2.0, .position = {0.5, 0.2}};
    const SourceSink b{.strength = 2.0, .position = {0.5, -2.2}}; // a reflected across y = -1
    const FlowModel single = buildFlowModel({.elements = {a, lower, upper}});
    const FlowModel both = buildFlowModel({.elements = {a, b, lower, upper}});
    CHECK(both.field.rows.size() == 2);
    const Vec2 p{1.0, 0.3};
    CHECK(both.field.velocity(p).x == doctest::Approx(single.field.velocity(p).x));
    CHECK(both.field.velocity(p).y == doctest::Approx(single.field.velocity(p).y));
}

TEST_CASE("channel: walls along one line are one wall, not a channel") {
    const FlowModel model = buildFlowModel({.elements = {SourceSink{.position = {0.0, 1.0}}, SymmetryWall{},
                                                         SymmetryWall{.point = {3.0, 0.0}}}});
    CHECK(model.imagesExact);
    CHECK(model.field.rows.empty());
    CHECK(model.images.size() == 1);
}

TEST_CASE("channel: every copy is a singular point, however far away") {
    const FlowField field = buildFlowModel({.elements = {SourceSink{.position = {0.0, 0.2}}, SymmetryWall{.point = {0.0, -1.0}},
                                                         SymmetryWall{.point = {0.0, 1.0}}}})
                                .field;
    CHECK(field.nearSingularity({0.0, 0.2}, 0.01));         // the source
    CHECK(field.nearSingularity({0.0, 0.2 + 40.0}, 0.01));  // ten periods up
    CHECK(field.nearSingularity({0.0, -2.2 - 40.0}, 0.01)); // the mirror row
    CHECK_FALSE(field.nearSingularity({0.5, 0.2}, 0.01));
    // Rows at y = 0.2 + 4n and -2.2 + 4n; the view holds 0.2, 1.8 and -2.2.
    CHECK(field.singularPoints({{-1.0, -2.5}, {1.0, 2.5}}).size() == 3);

    // A very narrow channel seen from far away stays bounded.
    const SymmetryWall lower{};
    const FlowField narrow = buildFlowModel({.elements = {SourceSink{.position = {0.0, 1e-6}}, lower,
                                                          SymmetryWall{.point = {0.0, 2e-6}}}})
                                 .field;
    CHECK(narrow.singularPoints({{-1.0, -1.0}, {1.0, 1.0}}).size() <= 512);
    CHECK(std::isfinite(length(narrow.velocity({1.0, 1e-6})))); // far past where cosh overflows
}

TEST_CASE("channel: streamlines never cross the walls") {
    const FlowField field = buildFlowModel({.elements = {SourceSink{.strength = 6.0, .position = {0.0, 0.2}},
                                                         SymmetryWall{.point = {0.0, -1.0}}, SymmetryWall{.point = {0.0, 1.0}}}})
                                .field;
    StreamlineSettings settings;
    settings.bounds = {{-4.0, -2.5}, {4.0, 2.5}};
    settings.separation = 0.25;
    const StreamlineSet set = computeStreamlines(field, settings);
    REQUIRE(set.lines.size() > 10);
    for (const std::vector<Vec2>& line : set.lines) {
        const auto strip = [](Vec2 p) { return p.y < -1.0 ? 0 : (p.y < 1.0 ? 1 : 2); };
        const bool sameStrip = std::ranges::all_of(line, [&](Vec2 p) {
            return std::abs(p.y + 1.0) < 1e-6 || std::abs(p.y - 1.0) < 1e-6 || strip(p) == strip(line[line.size() / 2]);
        });
        CHECK(sameStrip);
    }
}

TEST_CASE("channel: elements outside it are flagged") {
    const SymmetryWall lower{.point = {0.0, -1.0}};
    const SymmetryWall upper{.point = {0.0, 1.0}};
    const Scene inside{.elements = {SourceSink{.position = {0.0, 0.5}}, lower, upper}};
    CHECK(sceneWarnings(inside, buildFlowModel(inside)).empty());
    const Scene outside{.elements = {SourceSink{.position = {0.0, 1.5}}, lower, upper}};
    const std::vector<std::string> warnings = sceneWarnings(outside, buildFlowModel(outside));
    REQUIRE(warnings.size() == 1);
    CHECK(warnings.front().find("outside the channel") != std::string::npos);
}

TEST_CASE("warnings for coincident singularities and elements on walls") {
    const Scene clean{.elements = {UniformFlow{}, Doublet{.strength = 1.0}}};
    CHECK(sceneWarnings(clean, buildFlowModel(clean)).empty());

    const Scene coincident{.elements = {SourceSink{.position = {1.0, 1.0}}, Vortex{.position = {1.0, 1.0}}}};
    CHECK(sceneWarnings(coincident, buildFlowModel(coincident)).size() == 1);

    const Scene onWall{.elements = {Vortex{.position = {2.0, 0.0}}, SymmetryWall{}}};
    CHECK(sceneWarnings(onWall, buildFlowModel(onWall)).size() == 1);
}

TEST_CASE("streamlines of a uniform flow are straight, evenly spaced and point downstream") {
    const FlowField field{.elements = {UniformFlow{.speed = 2.0, .angle = 45.0 * deg}}};
    StreamlineSettings settings;
    settings.bounds = {{-2.0, -2.0}, {2.0, 2.0}};
    settings.separation = 0.25;
    const StreamlineSet set = computeStreamlines(field, settings);

    REQUIRE(set.lines.size() > 15);
    CHECK_FALSE(set.arrows.empty());
    for (const auto& line : set.lines) {
        REQUIRE(line.size() >= 3);
        const Vec2 d = line.back() - line.front();
        CHECK(cross(d, Vec2{1.0, 1.0}) == doctest::Approx(0.0).epsilon(1e-6).scale(length(d)));
        CHECK(dot(d, Vec2{1.0, 1.0}) > 0.0);
    }
    for (const Arrow& a : set.arrows) {
        CHECK(a.direction.x == doctest::Approx(std::sqrt(0.5)));
        CHECK(a.direction.y == doctest::Approx(std::sqrt(0.5)));
    }
}

TEST_CASE("streamlines of a vortex are closed circles") {
    const FlowField field{.elements = {Vortex{.circulation = 2.0 * pi}}};
    StreamlineSettings settings;
    settings.bounds = {{-2.0, -2.0}, {2.0, 2.0}};
    settings.separation = 0.3;
    const StreamlineSet set = computeStreamlines(field, settings);

    int closedCircles = 0;
    for (const auto& line : set.lines) {
        const double r = length(line.front());
        double drift = 0.0;
        for (const Vec2 p : line) drift = std::max(drift, std::abs(length(p) - r) / r);
        INFO("r = " << r << ", points = " << line.size());
        CHECK(drift < 1e-5); // RK4 truncation over one loop
        if (line.front() == line.back() && line.size() > 10) ++closedCircles;
    }
    CHECK(closedCircles >= 3);
}

namespace {

double nearestSquared(const StreamlineSet& set, Vec2 p) {
    double best = std::numeric_limits<double>::infinity();
    for (const auto& line : set.lines) {
        for (const Vec2 q : line) best = std::min(best, lengthSquared(q - p));
    }
    return best;
}

// Fraction of a's arrows with no arrow of b within `tolerance`.
double arrowsMoved(const StreamlineSet& a, const StreamlineSet& b, double tolerance) {
    std::size_t moved = 0;
    for (const Arrow& x : a.arrows) {
        double best = std::numeric_limits<double>::infinity();
        for (const Arrow& y : b.arrows) best = std::min(best, lengthSquared(x.position - y.position));
        if (best > tolerance * tolerance) ++moved;
    }
    return a.arrows.empty() ? 0.0 : static_cast<double>(moved) / static_cast<double>(a.arrows.size());
}

// Fraction of a's line points with no point of b within `tolerance`.
double linesMoved(const StreamlineSet& a, const StreamlineSet& b, double tolerance) {
    std::size_t moved = 0;
    std::size_t total = 0;
    for (const auto& line : a.lines) {
        for (std::size_t i = 0; i < line.size(); i += 3) {
            ++total;
            if (nearestSquared(b, line[i]) > tolerance * tolerance) ++moved;
        }
    }
    return total == 0 ? 0.0 : static_cast<double>(moved) / static_cast<double>(total);
}

FlowField liftingCylinder(double circulation) {
    return {.elements = {UniformFlow{.speed = 1.0}, Doublet{.strength = 1.0}, Vortex{.circulation = circulation}}};
}

StreamlineSettings plotLikeSettings() {
    StreamlineSettings settings;
    settings.bounds = {{-4.05, -2.55}, {4.05, 2.55}};
    settings.separation = 5.0 / 24.0;
    settings.singularRadius = 0.02;
    settings.arrowSpacing = 0.7;
    return settings;
}

} // namespace

TEST_CASE("streamlines from a set's own seeds, in the same flow, draw the same picture") {
    const FlowField field = liftingCylinder(5.0);
    const StreamlineSettings settings = plotLikeSettings();
    const StreamlineSet first = computeStreamlines(field, settings);
    REQUIRE(first.seeds.size() == first.lines.size());
    const StreamlineSet again = computeStreamlines(field, settings, first.seeds);
    // Every line and arrow comes back, traced from near the middle of the line now, so the points
    // differ by RK4 truncation at most.
    CHECK(linesMoved(first, again, 1e-3 * settings.separation) < 0.005);
    CHECK(arrowsMoved(first, again, 1e-3 * settings.separation) < 0.01);
    // Seeds need a little less room when following lines, so a few gaps may get a line.
    CHECK(linesMoved(again, first, 1e-3 * settings.separation) < 0.02);
    CHECK(arrowsMoved(again, first, 1e-3 * settings.separation) < 0.03);
}

TEST_CASE("with the previous seeds, streamlines and arrows follow a sweep smoothly") {
    // A lifting-cylinder sweep as the plot draws it, about one frame's change of circulation at a time.
    const StreamlineSettings settings = plotLikeSettings();
    const double tolerance = 0.5 * settings.separation;
    StreamlineSet kept = computeStreamlines(liftingCylinder(3.0), settings);
    StreamlineSet fresh = kept;
    double keptArrows = 0.0, freshArrows = 0.0, keptLines = 0.0;
    constexpr int frames = 30;
    for (int k = 1; k <= frames; ++k) {
        const FlowField field = liftingCylinder(3.0 + 0.042 * k);
        StreamlineSet nextKept = computeStreamlines(field, settings, kept.seeds);
        StreamlineSet nextFresh = computeStreamlines(field, settings);
        keptArrows += arrowsMoved(nextKept, kept, tolerance) / frames;
        freshArrows += arrowsMoved(nextFresh, fresh, tolerance) / frames;
        keptLines += linesMoved(nextKept, kept, 0.5 * tolerance) / frames;
        kept = std::move(nextKept);
        fresh = std::move(nextFresh);
    }
    INFO("arrows jumping per frame: " << keptArrows << " with the seeds, " << freshArrows << " seeded afresh");
    CHECK(keptArrows < 0.05);
    CHECK(keptArrows < 0.25 * freshArrows);
    CHECK(keptLines < 0.03);
}

TEST_CASE("arrows stay at the same place on a line when its ends move") {
    // Panning a uniform flow: the lines through the kept seeds are the same, only longer or shorter.
    const FlowField field{.elements = {UniformFlow{.speed = 1.0}}};
    StreamlineSettings settings;
    settings.bounds = {{-2.0, -2.0}, {2.0, 2.0}};
    settings.separation = 0.25;
    settings.arrowSpacing = 0.6;
    const StreamlineSet before = computeStreamlines(field, settings);
    settings.bounds = {{-1.63, -2.0}, {2.37, 2.0}};
    const StreamlineSet after = computeStreamlines(field, settings, before.seeds);

    StreamlineSet inside; // the arrows of `before` away from both views' ends
    for (const Arrow& a : before.arrows) {
        if (a.position.x > -1.4 && a.position.x < 1.8) inside.arrows.push_back(a);
    }
    REQUIRE(inside.arrows.size() > 20);
    CHECK(arrowsMoved(inside, after, 1e-9) == 0.0);
}

TEST_CASE("streamlines stay inside a sector's wedge") {
    const FlowField field{.elements = {SectorFlow{.strength = 1.0, .wall1 = 30.0 * deg, .wall2 = 120.0 * deg}}};
    StreamlineSettings settings;
    settings.bounds = {{-2.0, -2.0}, {2.0, 2.0}};
    settings.separation = 0.2;
    const StreamlineSet set = computeStreamlines(field, settings);
    REQUIRE_FALSE(set.lines.empty());
    for (const auto& line : set.lines) {
        for (const Vec2 p : line) CHECK(field.inDomain(p));
    }
}

TEST_CASE("equipotential lines cross the streamlines at right angles and keep phi constant") {
    // Cylinder of radius 1 in a uniform flow at 20 degrees.
    const double alpha = 20.0 * deg;
    const FlowField field{.elements = {UniformFlow{.speed = 1.0, .angle = alpha},
                                       Doublet{.strength = 1.0, .angle = alpha}}};
    StreamlineSettings settings;
    settings.bounds = {{-3.0, -3.0}, {3.0, 3.0}};
    settings.separation = 0.25;
    settings.arrowSpacing = 0.0;
    settings.orthogonal = true;
    const StreamlineSet set = computeStreamlines(field, settings);

    REQUIRE(set.lines.size() > 20);
    CHECK(set.arrows.empty());
    // Checked in the physical flow (outside the cylinder), away from the stagnation points on its
    // surface, where the direction field turns quickly and RK4 with a fixed step is least accurate.
    const auto checked = [](Vec2 p) { return length(p) > 1.2; };
    // phi changes by about U * separation between neighbouring lines; along a line it must stay put.
    const double spacing = 1.0 * settings.separation;
    int checkedLines = 0;
    for (const auto& line : set.lines) {
        if (!std::ranges::all_of(line, checked)) continue;
        ++checkedLines;
        const double phi0 = field.complexPotential(line.front()).real();
        double drift = 0.0;
        for (std::size_t i = 0; i + 1 < line.size(); ++i) {
            const Vec2 d = line[i + 1] - line[i];
            const Vec2 v = field.velocity((line[i] + line[i + 1]) * 0.5);
            if (length(d) <= 0.0) continue;
            CHECK(std::abs(dot(d / length(d), v / length(v))) < 1e-3); // each step crosses the flow at 90 deg
            drift = std::max(drift, std::abs(field.complexPotential(line[i + 1]).real() - phi0));
        }
        INFO("line from (" << line.front().x << ", " << line.front().y << "), " << line.size() << " points");
        CHECK(drift < 1e-4 * spacing);
    }
    CHECK(checkedLines > 10);
}

TEST_CASE("equipotential lines of a vortex are radial, despite phi's branch cut") {
    const Vec2 centre{0.5, -0.25};
    const FlowField field{.elements = {Vortex{.circulation = 3.0, .position = centre}}};
    StreamlineSettings settings;
    settings.bounds = {{-2.0, -2.0}, {2.0, 2.0}};
    settings.separation = 0.3;
    settings.arrowSpacing = 0.0;
    settings.orthogonal = true;
    const StreamlineSet set = computeStreamlines(field, settings);

    REQUIRE(set.lines.size() >= 8);
    for (const auto& line : set.lines) {
        const Vec2 ray = (line.back() - centre) / length(line.back() - centre);
        for (const Vec2 p : line) CHECK(std::abs(cross(p - centre, ray)) < 1e-9);
    }
}

TEST_CASE("speed grid: |V| inside the domain, NaN outside, infinite at singular points") {
    const ScalarGrid uniform = sampleSpeed({.elements = {UniformFlow{.speed = 2.0, .angle = 1.0}}},
                                           {{-1.0, -1.0}, {1.0, 1.0}}, 4, 3, 0.01);
    REQUIRE(uniform.values.size() == 12);
    for (const double v : uniform.values) CHECK(v == doctest::Approx(2.0));

    // 2 x 2 cells over [-1, 1]^2: centres at (+-0.5, +-0.5). Only (0.5, 0.5) is inside a 90 degree wedge.
    const ScalarGrid sector = sampleSpeed({.elements = {SectorFlow{.strength = 1.0, .wall1 = 0.0, .wall2 = 90.0 * deg}}},
                                          {{-1.0, -1.0}, {1.0, 1.0}}, 2, 2, 0.01);
    CHECK(std::isnan(sector.at(0, 0)));
    CHECK(std::isnan(sector.at(1, 0)));
    CHECK(std::isnan(sector.at(0, 1)));
    CHECK(sector.at(1, 1) == doctest::Approx(2.0 * length(Vec2{0.5, 0.5}))); // |W| = 2 V r for n = 2

    const ScalarGrid source = sampleSpeed({.elements = {SourceSink{.strength = 1.0}}}, {{-1.0, -1.0}, {1.0, 1.0}}, 1, 1, 0.01);
    CHECK(std::isinf(source.at(0, 0)));
}

TEST_CASE("field and line CSV exports") {
    const Scene scene{.elements = {UniformFlow{.speed = 2.0, .angle = 90.0 * deg},
                                   SectorFlow{.strength = 1.0, .vertex = {0.0, 0.0}, .wall1 = 0.0, .wall2 = 180.0 * deg}}};
    const FlowModel model = buildFlowModel(scene);
    const std::string csv = fieldCsv(model, {{-1.0, -1.0}, {1.0, 1.0}}, 3, 2);

    std::vector<std::string> rows;
    std::istringstream in(csv);
    for (std::string row; std::getline(in, row);) rows.push_back(row);
    REQUIRE(rows.size() == 7);
    CHECK(rows[0] == "x,y,u,v,speed,phi,psi,cp");
    // Bottom row (y = -1) is outside the upper half-plane wedge.
    CHECK(rows[1] == "-1,-1,NaN,NaN,NaN,NaN,NaN,NaN");
    // (1, 1): uniform (0, 2) plus sector z^1 (1, 0): u = 1, v = 2, phi = 2 y + x, psi = -2 x + y, Cp = 1 - 5/4.
    std::vector<double> values;
    std::istringstream row(rows[6]);
    for (std::string cell; std::getline(row, cell, ',');) values.push_back(std::stod(cell));
    const std::vector<double> expected{1.0, 1.0, 1.0, 2.0, std::sqrt(5.0), 3.0, -1.0, -0.25};
    REQUIRE(values.size() == expected.size());
    for (std::size_t i = 0; i < values.size(); ++i) CHECK(values[i] == doctest::Approx(expected[i]).epsilon(1e-15));

    StreamlineSet lines;
    lines.lines = {{{0.0, 0.0}, {1.0, 0.5}}};
    StreamlineSet equipotentials;
    equipotentials.lines = {{{2.0, 2.0}}, {{3.0, -1.0}}};
    CHECK(linesCsv(lines, equipotentials) ==
          "kind,line,x,y\n"
          "streamline,1,0,0\n"
          "streamline,1,1,0.5\n"
          "equipotential,1,2,2\n"
          "equipotential,2,3,-1\n");
}

TEST_CASE("scenes survive a save/load round trip") {
    Scene scene;
    scene.elements = {UniformFlow{.speed = 1.5, .angle = 30.0 * deg},
                      SourceSink{.strength = -2.0, .position = {1.0, -0.5}},
                      Vortex{.circulation = 3.0, .position = {0.0, 2.0}},
                      Doublet{.strength = 0.5, .position = {-1.0, 0.0}, .angle = 45.0 * deg},
                      SectorFlow{.strength = 2.0, .vertex = {0.5, 0.5}, .wall1 = 10.0 * deg, .wall2 = 100.0 * deg},
                      SymmetryWall{.point = {0.0, -2.0}, .angle = 5.0 * deg}};
    scene.plot.view = {{-3.0, -1.0}, {5.0, 4.0}};
    scene.plot.streamlineDensity = 31;
    scene.plot.showArrows = false;
    scene.plot.showParticles = false;
    scene.plot.showEquipotentials = true;
    scene.plot.background = Background::PressureCoefficient;
    scene.plot.showSurfacePressure = true;

    const SceneLoadResult loaded = sceneFromJson(sceneToJson(scene));
    REQUIRE_MESSAGE(loaded.scene.has_value(), loaded.error);
    const Scene& copy = *loaded.scene;
    REQUIRE(copy.elements.size() == scene.elements.size());
    for (std::size_t i = 0; i < scene.elements.size(); ++i) {
        INFO("element " << i);
        CHECK(copy.elements[i].index() == scene.elements[i].index());
        CHECK(describe(copy.elements[i]) == describe(scene.elements[i]));
    }
    // Same velocity field.
    const FlowModel a = buildFlowModel(scene);
    const FlowModel b = buildFlowModel(copy);
    for (const Vec2 p : {Vec2{0.7, 1.1}, Vec2{-0.3, 0.9}}) {
        CHECK(a.field.velocity(p).x == doctest::Approx(b.field.velocity(p).x));
        CHECK(a.field.velocity(p).y == doctest::Approx(b.field.velocity(p).y));
    }
    CHECK(copy.plot.view == scene.plot.view);
    CHECK(copy.plot.streamlineDensity == 31);
    CHECK_FALSE(copy.plot.showArrows);
    CHECK_FALSE(copy.plot.showParticles);
    CHECK(copy.plot.showEquipotentials);
    CHECK(copy.plot.background == Background::PressureCoefficient);
    CHECK(copy.plot.showSurfacePressure);
}

TEST_CASE("invalid scene files are rejected with a message") {
    CHECK_FALSE(sceneFromJson("not json").scene);
    CHECK_FALSE(sceneFromJson(R"({"format": "something-else", "version": 1})").scene);
    const SceneLoadResult badSector = sceneFromJson(
        R"({"format": "pfl-scene", "version": 1, "elements": [
            {"type": "sector", "strength": 1, "x": 0, "y": 0, "theta1_deg": 90, "theta2_deg": 10}]})");
    CHECK_FALSE(badSector.scene);
    CHECK(badSector.error.find("element #1") != std::string::npos);
    const SceneLoadResult missing =
        sceneFromJson(R"({"format": "pfl-scene", "version": 1, "elements": [{"type": "vortex", "x": 0, "y": 0}]})");
    CHECK_FALSE(missing.scene);
    CHECK(missing.error.find("circulation") != std::string::npos);
}

TEST_CASE("every example scene loads") {
    int count = 0;
    for (const auto& entry : std::filesystem::directory_iterator(PFL_EXAMPLES_DIR)) {
        if (entry.path().extension() != ".json") continue;
        INFO(entry.path().filename().string());
        const SceneLoadResult result = loadSceneFile(entry.path());
        CHECK_MESSAGE(result.scene.has_value(), result.error);
        ++count;
    }
    CHECK(count > 0);
}

TEST_CASE("built-in examples: every one loads and matches its file in examples/") {
    for (const ExampleScene& example : exampleScenes()) {
        INFO(example.file);
        const SceneLoadResult embedded = loadExampleScene(example);
        REQUIRE_MESSAGE(embedded.scene.has_value(), embedded.error);
        CHECK_FALSE(embedded.scene->elements.empty());
        const SceneLoadResult file = loadSceneFile(std::filesystem::path(PFL_EXAMPLES_DIR) / example.file);
        REQUIRE(file.scene);
        CHECK(sceneToJson(*embedded.scene) == sceneToJson(*file.scene));
    }
}

TEST_CASE("built-in examples: every file in examples/ is in the menu") {
    for (const auto& entry : std::filesystem::directory_iterator(PFL_EXAMPLES_DIR)) {
        if (entry.path().extension() != ".json") continue;
        const std::string name = entry.path().filename().string();
        INFO(name);
        const auto examples = exampleScenes();
        CHECK(std::any_of(examples.begin(), examples.end(), [&name](const ExampleScene& e) { return name == e.file; }));
    }
}

TEST_CASE("plot background: older files with speed_contours still load") {
    const SceneLoadResult old =
        sceneFromJson(R"({"format": "pfl-scene", "version": 1, "elements": [], "plot": {"speed_contours": true}})");
    REQUIRE(old.scene);
    CHECK(old.scene->plot.background == Background::Speed);
    CHECK_FALSE(old.scene->plot.showSurfacePressure);
    const SceneLoadResult bad =
        sceneFromJson(R"({"format": "pfl-scene", "version": 1, "elements": [], "plot": {"background": "pressure"}})");
    CHECK_FALSE(bad.scene);
}
