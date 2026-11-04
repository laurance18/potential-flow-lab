// dye.hpp: dye rakes and the streaklines they release.
//
// A rake is a straight segment carrying evenly spaced emitters. Every emitter releases a dye point
// at a fixed interval, and the points are advected with the flow (midpoint RK2, as the tracer
// particles). Joined in order of release they form a streakline: the curve through every fluid
// particle that has passed the emitter. In a steady flow it lies on the streamline through the
// emitter, so the dye shows the same curves as the streamlines; how far apart successive releases
// end up shows the speed. Each release has a serial number, so the dye can be drawn in pulses
// (alternating bands of equal release time) that bunch up where the flow is slow and stretch where
// it is fast.
//
// A point dies when it leaves the bounds, reaches a position the caller rejects (outside the flow
// domain, too close to a singular point), exceeds maxAge, or moves too fast to be followed (a step
// that would need more than maxSubsteps sub-steps, next to a singular point); a streakline is
// broken where its points have died. A streak stays fed by the same emitter index from step to step, so moving a rake drags
// its dye along; detach() cuts every streak loose (after rakes are added or removed, or their
// emitter count changes), and detached streaks drift away until their points die.
#pragma once

#include "core/particles.hpp"
#include "core/vec2.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace pfl {

struct DyeRake {
    Vec2 start{-3.0, -1.0};
    Vec2 end{-3.0, 1.0};
    int emitters = 9;

    friend bool operator==(const DyeRake&, const DyeRake&) = default;
};

inline constexpr int k_maxRakeEmitters = 64;

// Emitter positions: evenly spaced from start to end inclusive (the midpoint for one emitter).
std::vector<Vec2> emitterPositions(const DyeRake& rake);

// Why the rake is invalid, or nullopt if it is fine.
std::optional<std::string> dyeRakeError(const DyeRake& rake);

struct DyeEmitter {
    Vec2 position;
    std::size_t colour = 0; // the rake's index, for its colour
};

struct DyePoint {
    Vec2 pos;
    double age = 0.0;
    std::uint32_t serial = 0; // release number within its streak
    bool alive = true;
};

struct DyeStreak {
    std::size_t colour = 0;
    bool attached = true;          // still fed by an emitter
    std::deque<DyePoint> points;   // oldest first
};

struct DyeSettings {
    Rect bounds{{-4.0, -2.5}, {4.0, 2.5}};
    double maxAge = 12.0;
    double releaseInterval = 1.0 / 30.0;
    double maxStepLength = 0.01;
    int maxSubsteps = 16;
    int maxReleasesPerStep = 8; // a long step releases at most this many points per emitter
};

class DyeSystem {
public:
    // Advances the dye by dt and releases new points from `emitters` (in the same order every step
    // while they are attached). `velocity` is any callable Vec2 -> Vec2 and `isValid` any callable
    // Vec2 -> bool; an emitter at an invalid position releases nothing.
    template <typename VelocityFn, typename ValidFn>
    void step(double dt, std::span<const DyeEmitter> emitters, const VelocityFn& velocity, const ValidFn& isValid);

    // Stops feeding every streak; the next step starts new ones.
    void detach();
    void clear();

    DyeSettings& settings() { return m_settings; }
    const DyeSettings& settings() const { return m_settings; }
    std::span<const DyeStreak> streaks() const { return m_streaks; }

private:
    DyeSettings m_settings;
    std::vector<DyeStreak> m_streaks; // attached ones first, in emitter order
    std::size_t m_attached = 0;
    double m_sinceRelease = 0.0;
};

template <typename VelocityFn, typename ValidFn>
void DyeSystem::step(double dt, std::span<const DyeEmitter> emitters, const VelocityFn& velocity,
                     const ValidFn& isValid) {
    const DyeSettings& s = m_settings;
    // Advects p by `time`; false if it has to die.
    const auto move = [&](DyePoint& p, double time) {
        bool resolved = true;
        if (time > 0.0) advect(p.pos, time, velocity, s.maxStepLength, s.maxSubsteps, &resolved);
        return resolved && std::isfinite(p.pos.x) && std::isfinite(p.pos.y) && s.bounds.contains(p.pos) &&
               p.age <= s.maxAge && isValid(p.pos);
    };

    // Existing points.
    for (DyeStreak& streak : m_streaks) {
        for (DyePoint& p : streak.points) {
            if (!p.alive) continue;
            p.age += dt;
            p.alive = move(p, dt);
        }
    }

    // One attached streak per emitter, in order. Extra attached streaks (fewer emitters than last
    // step) are cut loose.
    if (emitters.size() < m_attached) {
        for (std::size_t i = emitters.size(); i < m_attached; ++i) m_streaks[i].attached = false;
        m_attached = emitters.size();
    }
    while (m_attached < emitters.size()) {
        m_streaks.insert(m_streaks.begin() + static_cast<std::ptrdiff_t>(m_attached), DyeStreak{});
        ++m_attached;
    }
    for (std::size_t i = 0; i < m_attached; ++i) m_streaks[i].colour = emitters[i].colour;

    // Releases during this step, each advected for the time since its release.
    m_sinceRelease += dt;
    const int releases = s.releaseInterval > 0.0 ? static_cast<int>(m_sinceRelease / s.releaseInterval) : 0;
    if (releases > 0) {
        m_sinceRelease -= releases * s.releaseInterval;
        const int released = std::min(releases, s.maxReleasesPerStep);
        for (std::size_t i = 0; i < m_attached; ++i) {
            const Vec2 origin = emitters[i].position;
            if (!isValid(origin)) continue;
            DyeStreak& streak = m_streaks[i];
            for (int k = released - 1; k >= 0; --k) { // oldest release first
                DyePoint p;
                p.pos = origin;
                p.age = m_sinceRelease + k * s.releaseInterval;
                p.serial = streak.points.empty() ? 0 : streak.points.back().serial + 1;
                p.alive = move(p, p.age);
                streak.points.push_back(p);
            }
        }
    }

    // Drop dead points from the old end, and detached streaks with nothing left.
    for (DyeStreak& streak : m_streaks) {
        while (!streak.points.empty() && !streak.points.front().alive) streak.points.pop_front();
    }
    std::erase_if(m_streaks, [](const DyeStreak& streak) { return !streak.attached && streak.points.empty(); });
}

} // namespace pfl
