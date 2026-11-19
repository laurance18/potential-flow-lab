#include "core/streamline.hpp"

#include "core/parallel.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <numbers>
#include <optional>
#include <span>
#include <utility>

namespace pfl {
namespace {

constexpr double k_stopFactor = 0.5;    // a line stops this fraction of a separation from another line
constexpr double k_seedFactor = 0.95;   // a seed needs this fraction of a separation of free space
constexpr double k_keepFactor = 0.8;    // ... or this fraction while following previous lines
constexpr double k_slideSeparations = 3.0; // how far along its streamline a crowded kept seed may move
constexpr double k_stepsPerSeparation = 6.0;
constexpr double k_stagnationFraction = 1e-6;
constexpr double k_nearStepFactor = 0.15; // step <= this fraction of the distance to a singular point
constexpr int k_ringSeeds = 8; // seeds placed around each singular point

// Uniform grid of cells (one separation wide) holding the points of finished streamlines.
class OccupancyGrid {
public:
    OccupancyGrid(const Rect& bounds, double cell) : m_origin(bounds.min), m_cell(cell) {
        const Vec2 size = bounds.size();
        m_nx = std::max(1, static_cast<int>(std::ceil(size.x / cell)));
        m_ny = std::max(1, static_cast<int>(std::ceil(size.y / cell)));
        m_cells.resize(static_cast<std::size_t>(m_nx) * static_cast<std::size_t>(m_ny));
    }

    void insert(Vec2 p) {
        const auto [i, j] = cellOf(p);
        if (i >= 0 && i < m_nx && j >= 0 && j < m_ny) m_cells[index(i, j)].push_back(p);
    }

    // True if any stored point lies within `radius` (at most one cell) of p.
    bool near(Vec2 p, double radius) const {
        const auto [ci, cj] = cellOf(p);
        const double r2 = radius * radius;
        for (int j = std::max(cj - 1, 0); j <= std::min(cj + 1, m_ny - 1); ++j) {
            for (int i = std::max(ci - 1, 0); i <= std::min(ci + 1, m_nx - 1); ++i) {
                for (const Vec2 q : m_cells[index(i, j)]) {
                    if (lengthSquared(q - p) < r2) return true;
                }
            }
        }
        return false;
    }

private:
    std::pair<int, int> cellOf(Vec2 p) const {
        return {static_cast<int>(std::floor((p.x - m_origin.x) / m_cell)),
                static_cast<int>(std::floor((p.y - m_origin.y) / m_cell))};
    }
    std::size_t index(int i, int j) const {
        return static_cast<std::size_t>(j) * static_cast<std::size_t>(m_nx) + static_cast<std::size_t>(i);
    }

    Vec2 m_origin;
    double m_cell;
    int m_nx = 1;
    int m_ny = 1;
    std::vector<std::vector<Vec2>> m_cells;
};

// A line traced from a seed with no other lines around; the last point is the seed again if it
// closed on itself, and `closing` is where that last step really ended.
struct Trajectory {
    std::vector<Vec2> points;
    Vec2 closing;
    bool closed = false;
};

// Both directions from a seed, traced ahead of time (the backward one only if the forward one does
// not close on itself).
struct Speculation {
    Trajectory forward;
    Trajectory backward;
    bool hasBackward = false;
};

class Tracer {
public:
    Tracer(const FlowField& field, const StreamlineSettings& settings, double minSpeed)
        : m_field(field), m_settings(settings), m_minSpeed(minSpeed) {
        for (const FlowElement& e : field.elements) {
            if (const std::optional<Vec2> centre = singularPoint(e)) m_centres.push_back(*centre);
        }
        // The flow turns round a thin airfoil's leading edge in a tiny space.
        for (const VortexSheet& sheet : field.vortexSheets) m_centres.push_back(sheet.leadingEdge());
    }

    // Step of at most |h|, shortened near singular points where streamlines curve tightly.
    double localStep(Vec2 p, double h) const {
        double distance = std::abs(h) / k_nearStepFactor;
        for (const Vec2 c : m_centres) distance = std::min(distance, length(p - c));
        for (const PeriodicRow& row : m_field.rows) distance = std::min(distance, length(p - nearestCopy(row, p)));
        return std::copysign(k_nearStepFactor * distance, h);
    }

    // Unit flow direction at p, or nullopt where a line must stop.
    std::optional<Vec2> direction(Vec2 p) const {
        if (!m_settings.bounds.contains(p) || !m_field.inDomain(p) ||
            m_field.nearSingularity(p, m_settings.singularRadius)) {
            return std::nullopt;
        }
        const Vec2 v = m_field.velocity(p);
        const double speed = length(v);
        if (!(speed > m_minSpeed) || !std::isfinite(speed)) return std::nullopt;
        return m_settings.orthogonal ? perp(v / speed) : v / speed;
    }

    // One RK4 step of arc length h (negative h integrates upstream).
    std::optional<Vec2> step(Vec2 p, double h) const {
        const auto k1 = direction(p);
        if (!k1) return std::nullopt;
        const auto k2 = direction(p + *k1 * (0.5 * h));
        if (!k2) return std::nullopt;
        const auto k3 = direction(p + *k2 * (0.5 * h));
        if (!k3) return std::nullopt;
        const auto k4 = direction(p + *k3 * h);
        if (!k4 || dot(*k1, *k4) < 0.0) return std::nullopt; // stepped across a stagnation point
        const Vec2 next = p + (*k1 + *k2 * 2.0 + *k3 * 2.0 + *k4) * (h / 6.0);
        if (!m_field.vortexSheets.empty() && m_field.crossesSheet(p, next)) return std::nullopt; // cut through a thin airfoil
        return next;
    }

    // The line from `seed` (included) in one direction, stopping near the lines in `grid`: other
    // lines only decide where a line stops, not where it goes. With the grid as it is before more
    // lines are added, the result holds every point the line keeps later (cut() with the fuller
    // grid). Safe to run on several threads at once while nobody changes the grid.
    Trajectory run(Vec2 seed, double h, const OccupancyGrid& grid) const {
        Trajectory t;
        t.points = {seed};
        const double stopDistance = k_stopFactor * m_settings.separation;
        double travelled = 0.0;
        while (static_cast<int>(t.points.size()) < m_settings.maxPointsPerLine) {
            const double hLocal = localStep(t.points.back(), h);
            const std::optional<Vec2> next = step(t.points.back(), hLocal);
            if (!next || grid.near(*next, stopDistance)) break;
            travelled += std::abs(hLocal);
            if (travelled > 3.0 * m_settings.separation && length(*next - seed) < std::abs(hLocal)) {
                t.points.push_back(seed); // closed on itself
                t.closing = *next;
                t.closed = true;
                break;
            }
            t.points.push_back(*next);
        }
        return t;
    }

    // The part of `t` that is kept with the lines in `grid`: up to its first point within half a
    // separation of one of them. Sets `closed` if the line still returns to its seed.
    std::vector<Vec2> cut(const Trajectory& t, const OccupancyGrid& grid, bool& closed) const {
        const double stopDistance = k_stopFactor * m_settings.separation;
        closed = false;
        for (std::size_t k = 1; k < t.points.size(); ++k) {
            const bool closing = t.closed && k + 1 == t.points.size();
            if (grid.near(closing ? t.closing : t.points[k], stopDistance)) {
                return {t.points.begin(), t.points.begin() + static_cast<std::ptrdiff_t>(k)};
            }
        }
        closed = t.closed;
        return t.points;
    }

    // Points from `seed` (included) in one direction, stopping near the lines in `grid`. Sets
    // `closed` if the line returned to its seed.
    std::vector<Vec2> trace(Vec2 seed, double h, const OccupancyGrid& grid, bool& closed) const {
        const Trajectory t = run(seed, h, grid);
        closed = t.closed;
        return t.points;
    }

    // The nearest point along the streamline through `seed`, within `maxArc` either way, that has
    // `room` of free space, with its signed arc length from the seed (> 0 downstream).
    std::optional<std::pair<Vec2, double>> slideToRoom(Vec2 seed, double h, const OccupancyGrid& grid, double room,
                                                       double maxArc) const {
        std::optional<Vec2> down = seed;
        std::optional<Vec2> up = seed;
        double downArc = 0.0;
        double upArc = 0.0;
        while ((down && downArc < maxArc) || (up && upArc < maxArc)) {
            if (down && downArc < maxArc) {
                const double step = localStep(*down, h);
                down = this->step(*down, step);
                downArc += step;
                if (down && !grid.near(*down, room)) return std::pair{*down, downArc};
            }
            if (up && upArc < maxArc) {
                const double step = localStep(*up, -h);
                up = this->step(*up, step);
                upArc -= step;
                if (up && !grid.near(*up, room)) return std::pair{*up, -upArc};
            }
        }
        return std::nullopt;
    }

private:
    const FlowField& m_field;
    const StreamlineSettings& m_settings;
    double m_minSpeed;
    std::vector<Vec2> m_centres;
};

// Arrows every `spacing` of arc length, counted from the seed (line[seedIndex]) both ways and
// shifted by `phase` in [0, 1) so neighbouring lines don't line them up. Anchored at the seed, they
// stay put when the line's ends grow or shrink. None within `margin` of an end; a line too short
// for that gets one arrow in its middle.
void addArrows(const std::vector<Vec2>& line, std::size_t seedIndex, double spacing, double minLength, double margin,
               double phase, std::vector<Arrow>& arrows) {
    std::vector<double> arc(line.size(), 0.0);
    for (std::size_t i = 1; i < line.size(); ++i) arc[i] = arc[i - 1] + length(line[i] - line[i - 1]);
    const double total = arc.back();
    if (total < minLength || spacing <= 0.0) return;

    const double origin = arc[seedIndex] + phase * spacing;
    std::vector<double> targets;
    for (double k = std::ceil((margin - origin) / spacing); origin + k * spacing <= total - margin; k += 1.0) {
        targets.push_back(origin + k * spacing);
    }
    if (targets.empty()) targets.push_back(0.5 * total);

    std::size_t i = 1;
    for (const double target : targets) {
        while (i + 1 < line.size() && arc[i] < target) ++i;
        const Vec2 segment = line[i] - line[i - 1];
        const double l = arc[i] - arc[i - 1];
        if (l <= 0.0) continue;
        const double t = std::clamp((target - arc[i - 1]) / l, 0.0, 1.0);
        arrows.push_back({line[i - 1] + segment * t, segment / l});
    }
}

// The seed recorded for an open line: the point nearest its middle that has a new seed's room from
// the lines before it (`hasRoom`), or else the seed it was started from, with the arrow phase
// shifted so the arrows stay where they are. A seed near the middle stays in the view when it pans
// and away from the line's ends, where the lines around it change most; the room means it is not
// crowded out next time by the same lines in a slightly different place.
template <class HasRoom>
StreamlineSeed recordedSeed(const std::vector<Vec2>& line, std::size_t seedIndex, double spacing, double phase,
                            const HasRoom& hasRoom) {
    std::vector<double> arc(line.size(), 0.0);
    for (std::size_t i = 1; i < line.size(); ++i) arc[i] = arc[i - 1] + length(line[i] - line[i - 1]);
    const double half = 0.5 * arc.back();
    std::size_t chosen = seedIndex;
    double best = std::abs(arc[seedIndex] - half);
    for (std::size_t i = 0; i < line.size(); ++i) {
        const double distance = std::abs(arc[i] - half);
        if (distance < best && hasRoom(line[i])) {
            chosen = i;
            best = distance;
        }
    }
    if (spacing <= 0.0) return {line[chosen], phase};
    const double shifted = phase - (arc[chosen] - arc[seedIndex]) / spacing;
    return {line[chosen], shifted - std::floor(shifted)};
}

} // namespace

StreamlineSet computeStreamlines(const FlowField& field, const StreamlineSettings& settings,
                                 std::span<const StreamlineSeed> previousSeeds) {
    StreamlineSet result;
    const Vec2 size = settings.bounds.size();
    if (!(settings.separation > 0.0) || !(size.x > 0.0) || !(size.y > 0.0)) return result;

    const double h = settings.separation / k_stepsPerSeparation;
    const double minSpeed = k_stagnationFraction * typicalSpeed(field, settings.bounds);
    const Tracer tracer(field, settings, minSpeed);
    OccupancyGrid grid(settings.bounds, settings.separation);

    // Fallback seeds: rings around singular points first (so sources and vortices get lines),
    // then a regular grid over the bounds.
    std::vector<Vec2> seeds;
    for (const Vec2 centre : field.singularPoints(settings.bounds)) {
        const double r = std::max(2.0 * settings.singularRadius, 0.5 * settings.separation);
        for (int k = 0; k < k_ringSeeds; ++k) {
            const double a = 2.0 * std::numbers::pi * (k + 0.5) / k_ringSeeds;
            seeds.push_back(centre + Vec2{r * std::cos(a), r * std::sin(a)});
        }
    }
    for (double y = settings.bounds.min.y + 0.5 * settings.separation; y < settings.bounds.max.y;
         y += settings.separation) {
        for (double x = settings.bounds.min.x + 0.5 * settings.separation; x < settings.bounds.max.x;
             x += settings.separation) {
            seeds.push_back({x, y});
        }
    }

    std::deque<Vec2> candidates; // seeds beside finished lines, tried before the fallback seeds
    const int seedStride = std::max(1, static_cast<int>(k_stepsPerSeparation));
    // Following previous lines, a seed needs less room: a kept line is not crowded out by a small
    // move of its neighbours, and the gaps opening where lines drift apart are filled (with the full
    // room the picture would end up about a fifth sparser than seeded afresh).
    const double room = (previousSeeds.empty() ? k_seedFactor : k_keepFactor) * settings.separation;

    // `phase`: the arrow offset of a kept seed; nullopt for a new seed. `ahead`: the line traced from
    // `seed` ahead of time with no other lines around, which the lines in the grid now cut short.
    auto tryLine = [&](Vec2 seed, std::optional<double> phase, const Speculation* ahead) {
        if (!tracer.direction(seed)) return;
        if (grid.near(seed, room)) {
            // A kept seed crowded by the lines before it: the same streamline may still have room a
            // little further along. Starting there keeps the line; the phase keeps its arrows in place.
            if (!phase) return;
            const auto slid = tracer.slideToRoom(seed, h, grid, room, k_slideSeparations * settings.separation);
            if (!slid) return;
            seed = slid->first;
            ahead = nullptr; // a different seed
            if (settings.arrowSpacing > 0.0) {
                const double shifted = *phase - slid->second / settings.arrowSpacing;
                phase = shifted - std::floor(shifted);
            }
        }
        bool closed = false;
        std::vector<Vec2> forward = ahead ? tracer.cut(ahead->forward, grid, closed) : tracer.trace(seed, h, grid, closed);
        std::vector<Vec2> line;
        if (!closed) {
            bool ignored = false;
            std::vector<Vec2> backward = ahead && ahead->hasBackward ? tracer.cut(ahead->backward, grid, ignored)
                                                                     : tracer.trace(seed, -h, grid, ignored);
            line.assign(backward.rbegin(), backward.rend() - 1); // upstream part, without the seed
        }
        const std::size_t seedIndex = line.size();
        line.insert(line.end(), forward.begin(), forward.end());
        if (line.size() < 3) return;

        const double arrowPhase =
            phase.value_or(std::fmod(0.618034 * static_cast<double>(result.lines.size()), 1.0)); // golden ratio
        // Recorded before the line itself is in the grid. A closed line keeps its seed, where the
        // spacing of its arrows has its seam.
        const auto hasRoom = [&](Vec2 p) { return !grid.near(p, k_seedFactor * settings.separation); };
        result.seeds.push_back(closed ? StreamlineSeed{seed, arrowPhase}
                                      : recordedSeed(line, seedIndex, settings.arrowSpacing, arrowPhase, hasRoom));

        for (const Vec2 p : line) grid.insert(p);
        for (std::size_t i = 0; i + 1 < line.size(); i += static_cast<std::size_t>(seedStride)) {
            const Vec2 d = line[i + 1] - line[i];
            const double l = length(d);
            if (l <= 0.0) continue;
            const Vec2 offset = perp(d / l) * settings.separation;
            candidates.push_back(line[i] + offset);
            candidates.push_back(line[i] - offset);
        }
        addArrows(line, seedIndex, settings.arrowSpacing, 0.5 * settings.separation, 0.25 * settings.separation,
                  arrowPhase, result.arrows);
        result.lines.push_back(std::move(line));
    };

    // The previous lines first, in their order, so the same lines come back. Each is traced ahead of
    // time on every core, stopping only at the edges, the bodies and the singular points; then they
    // are tried one after another as before, each cut short by the lines found before it, which gives
    // exactly the lines tried one at a time would. Nearly every kept seed keeps its line, so little
    // of that work is lost, and this is what runs every frame while something is dragged.
    std::vector<Speculation> traced(previousSeeds.size());
    parallelFor(previousSeeds.size(), 1, [&](std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
            const Vec2 seed = previousSeeds[i].position;
            traced[i].forward = tracer.run(seed, h, grid);
            traced[i].hasBackward = !traced[i].forward.closed;
            if (traced[i].hasBackward) traced[i].backward = tracer.run(seed, -h, grid);
        }
    });
    for (std::size_t i = 0; i < previousSeeds.size(); ++i) tryLine(previousSeeds[i].position, previousSeeds[i].arrowPhase, &traced[i]);
    traced.clear();

    // Then the gaps, one seed at a time: the candidates beside the lines, oldest first, and the
    // fallback seeds whenever there are none. (Traced ahead in parallel they would mostly be wasted:
    // the candidates beside a line lie along one streamline, which the first of them takes.)
    std::size_t nextSeed = 0;
    while (true) {
        if (!candidates.empty()) {
            const Vec2 seed = candidates.front();
            candidates.pop_front();
            tryLine(seed, std::nullopt, nullptr);
        } else if (nextSeed < seeds.size()) {
            tryLine(seeds[nextSeed++], std::nullopt, nullptr);
        } else {
            break;
        }
    }
    return result;
}

} // namespace pfl
