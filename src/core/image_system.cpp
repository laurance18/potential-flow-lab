#include "core/image_system.hpp"

#include <algorithm>
#include <numbers>

namespace pfl {
namespace {

constexpr double k_tolerance = 1e-9;
constexpr int k_maxGenerations = 16;
constexpr std::size_t k_maxImagesPerElement = 32;

bool close(double a, double b) { return std::abs(a - b) <= k_tolerance * (1.0 + std::max(std::abs(a), std::abs(b))); }
bool close(Vec2 a, Vec2 b) { return close(a.x, b.x) && close(a.y, b.y); }
bool sameAngle(double a, double b) { return std::abs(std::remainder(a - b, 2.0 * std::numbers::pi)) <= k_tolerance; }

// T = 2 h n when the walls are exactly two distinct parallel lines (walls along the same line count
// once) a distance h apart, with n the normal of the first wall; otherwise nullopt.
std::optional<Vec2> channelPeriod(std::span<const SymmetryWall> walls) {
    const SymmetryWall& first = walls.front();
    std::optional<double> gap;
    for (const SymmetryWall& wall : walls) {
        if (std::abs(std::sin(wall.angle - first.angle)) > k_tolerance) return std::nullopt;
        const double offset = dot(wall.point - first.point, first.normal());
        if (close(offset, 0.0)) continue;
        if (gap && !close(offset, *gap)) return std::nullopt;
        gap = offset;
    }
    if (!gap) return std::nullopt;
    return first.normal() * (2.0 * std::abs(*gap));
}

// True if the two point elements are copies in one row of period T.
bool sameRow(const FlowElement& a, const FlowElement& b, Vec2 period) {
    const Vec2 gap = *singularPoint(a) - *singularPoint(b);
    return sameElement(a, translated(b, period * std::round(dot(gap, period) / lengthSquared(period))));
}

// The element's own row and its mirror row. No mirror row for an element on a wall that is its own
// image, and none for an element whose image the scene already holds as another element: the two
// are then a pair, as for a single wall.
std::vector<PeriodicRow> channelRows(std::span<const FlowElement> elements, const SymmetryWall& wall, Vec2 period) {
    std::vector<PeriodicRow> rows;
    std::vector<bool> paired(elements.size(), false);
    for (std::size_t i = 0; i < elements.size(); ++i) {
        if (paired[i]) continue;
        rows.push_back({.element = elements[i], .period = period, .original = true});
        const FlowElement image = *mirror(elements[i], wall);
        if (sameRow(image, elements[i], period)) continue;
        bool found = false;
        for (std::size_t j = i + 1; j < elements.size() && !found; ++j) {
            if (paired[j] || !sameRow(image, elements[j], period)) continue;
            paired[j] = true;
            found = true;
            rows.push_back({.element = elements[j], .period = period, .original = true});
        }
        if (!found) rows.push_back({.element = image, .period = period});
    }
    return rows;
}

} // namespace

Vec2 reflect(Vec2 p, const SymmetryWall& wall) {
    const Vec2 t = wall.direction();
    const Vec2 d = p - wall.point;
    return wall.point + t * (2.0 * dot(d, t)) - d;
}

std::optional<FlowElement> mirror(const FlowElement& element, const SymmetryWall& wall) {
    if (const auto* e = std::get_if<SourceSink>(&element)) {
        return SourceSink{.strength = e->strength, .position = reflect(e->position, wall)};
    }
    if (const auto* e = std::get_if<Vortex>(&element)) {
        return Vortex{.circulation = -e->circulation, .position = reflect(e->position, wall)};
    }
    if (const auto* e = std::get_if<Doublet>(&element)) {
        return Doublet{.strength = e->strength, .position = reflect(e->position, wall), .angle = 2.0 * wall.angle - e->angle};
    }
    return std::nullopt;
}

bool sameElement(const FlowElement& a, const FlowElement& b) {
    if (a.index() != b.index()) return false;
    if (const auto* x = std::get_if<UniformFlow>(&a)) {
        const auto& y = std::get<UniformFlow>(b);
        return close(x->speed, y.speed) && sameAngle(x->angle, y.angle);
    }
    if (const auto* x = std::get_if<SourceSink>(&a)) {
        const auto& y = std::get<SourceSink>(b);
        return close(x->strength, y.strength) && close(x->position, y.position);
    }
    if (const auto* x = std::get_if<Vortex>(&a)) {
        const auto& y = std::get<Vortex>(b);
        return close(x->circulation, y.circulation) && close(x->position, y.position);
    }
    if (const auto* x = std::get_if<Doublet>(&a)) {
        const auto& y = std::get<Doublet>(b);
        return close(x->strength, y.strength) && close(x->position, y.position) && sameAngle(x->angle, y.angle);
    }
    const auto& x = std::get<SectorFlow>(a);
    const auto& y = std::get<SectorFlow>(b);
    return close(x.strength, y.strength) && close(x.vertex, y.vertex) && close(x.wall1, y.wall1) &&
           close(x.wall2, y.wall2);
}

ImageSystem buildImageSystem(std::span<const FlowElement> elements, std::span<const SymmetryWall> walls) {
    ImageSystem result;
    if (walls.empty()) return result;

    // Every distinct copy so far (originals first), and the copies added in the last generation.
    std::vector<FlowElement> all;
    for (const FlowElement& e : elements) {
        if (mirror(e, walls.front())) all.push_back(e);
    }
    if (const std::optional<Vec2> period = channelPeriod(walls)) {
        result.rows = channelRows(all, walls.front(), *period);
        return result;
    }
    const std::size_t originals = all.size();
    const std::size_t cap = originals * (k_maxImagesPerElement + 1);
    std::vector<FlowElement> frontier = all;

    for (int generation = 0; generation < k_maxGenerations && !frontier.empty(); ++generation) {
        std::vector<FlowElement> next;
        for (const FlowElement& e : frontier) {
            for (const SymmetryWall& wall : walls) {
                FlowElement image = *mirror(e, wall);
                const bool known = std::ranges::any_of(all, [&](const FlowElement& o) { return sameElement(o, image); });
                if (known || all.size() >= cap) continue;
                all.push_back(image);
                next.push_back(std::move(image));
            }
        }
        frontier = std::move(next);
        if (all.size() >= cap) break;
    }
    // Leftover copies that were never reflected (or reflections dropped at the cap) mean the
    // reflection group did not close.
    result.exact = frontier.empty() && all.size() < cap;
    result.images.assign(all.begin() + static_cast<std::ptrdiff_t>(originals), all.end());
    return result;
}

UniformFlow constrainToWalls(const UniformFlow& flow, std::span<const SymmetryWall> walls) {
    if (walls.empty()) return flow;
    const double phi = walls.front().angle;
    for (const SymmetryWall& wall : walls) {
        if (std::abs(std::sin(wall.angle - phi)) > k_tolerance) return UniformFlow{.speed = 0.0, .angle = flow.angle};
    }
    const double along = flow.speed * std::cos(flow.angle - phi);
    if (along >= 0.0) return UniformFlow{.speed = along, .angle = phi};
    return UniformFlow{.speed = -along, .angle = phi + std::numbers::pi};
}

} // namespace pfl
