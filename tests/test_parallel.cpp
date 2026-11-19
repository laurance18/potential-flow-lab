#include "core/dye.hpp"
#include "core/parallel.hpp"
#include "core/particles.hpp"
#include "core/scene.hpp"
#include "core/streamline.hpp"

#include <doctest/doctest.h>

#include <atomic>
#include <cmath>
#include <numbers>
#include <vector>

using namespace pfl;

TEST_CASE("parallelFor covers the range exactly once") {
    for (const std::size_t count : {std::size_t{0}, std::size_t{1}, std::size_t{7}, std::size_t{100}, std::size_t{10007}}) {
        for (const std::size_t grain : {std::size_t{1}, std::size_t{16}, std::size_t{1000}}) {
            std::vector<std::atomic<int>> hits(count);
            parallelFor(count, grain, [&](std::size_t begin, std::size_t end) {
                CHECK(begin < end);
                for (std::size_t i = begin; i < end; ++i) hits[i].fetch_add(1);
            });
            for (std::size_t i = 0; i < count; ++i) CHECK(hits[i].load() == 1);
        }
    }
}

namespace {

// A source-panel body and a thin airfoil in a stream: the evaluations the threads share.
FlowModel panelScene() {
    return buildFlowModel(Scene{.elements = {UniformFlow{.speed = 1.0, .angle = 0.1},
                                             PanelBody{.position = {0.5, -0.8}, .semiAxisA = 0.8, .semiAxisB = 0.3, .panels = 60},
                                             ThinAirfoil{.position = {-0.5, 0.6}, .panels = 50}}});
}

} // namespace

TEST_CASE("particles and dye on many threads give the same result every time") {
    const FlowModel model = panelScene();
    const FlowField& field = model.field;
    const auto velocity = [&field](Vec2 p) { return field.velocity(p); };
    const auto valid = [&field](Vec2 p) { return field.inDomain(p) && !field.nearSingularity(p, 0.05); };
    const auto run = [&] {
        ParticleSystem particles(3000, ParticleSettings{.bounds = {{-4.0, -2.5}, {4.0, 2.5}}, .maxAge = 2.0}, 77);
        particles.respawnAll(valid);
        for (int i = 0; i < 40; ++i) particles.step(1.0 / 30.0, velocity, valid);
        return std::vector<Particle>(particles.particles().begin(), particles.particles().end());
    };
    const std::vector<Particle> a = run();
    const std::vector<Particle> b = run();
    REQUIRE(a.size() == b.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        CHECK(a[i].pos == b[i].pos);
        CHECK(a[i].age == b[i].age);
    }

    const auto dye = [&] {
        DyeSystem system;
        const std::vector<DyeEmitter> emitters = {{{-3.5, 0.5}, 0}, {{-3.5, -0.2}, 0}, {{-3.5, 1.1}, 1}};
        for (int i = 0; i < 120; ++i) system.step(1.0 / 30.0, emitters, velocity, valid);
        std::vector<Vec2> points;
        for (const DyeStreak& s : system.streaks()) {
            for (const DyePoint& p : s.points) points.push_back(p.pos);
        }
        return points;
    };
    CHECK(dye() == dye());
}

TEST_CASE("streamlines from kept seeds, traced ahead on many threads, are the same every time") {
    const FlowModel model = panelScene();
    StreamlineSettings settings;
    settings.bounds = {{-4.2, -2.7}, {4.2, 2.7}};
    settings.separation = 5.0 / 24.0;
    const StreamlineSet first = computeStreamlines(model.field, settings);
    REQUIRE(first.lines.size() > 10);
    const StreamlineSet a = computeStreamlines(model.field, settings, first.seeds);
    const StreamlineSet b = computeStreamlines(model.field, settings, first.seeds);
    REQUIRE(a.lines.size() == b.lines.size());
    for (std::size_t i = 0; i < a.lines.size(); ++i) CHECK(a.lines[i] == b.lines[i]);
    // In the same flow the kept seeds bring back the same picture.
    CHECK(a.lines.size() >= first.lines.size() - 2);
}
