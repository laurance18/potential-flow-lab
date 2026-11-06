#include "core/flow_field.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace pfl {
namespace {

bool within(const FlowElement& element, Vec2 p, double radius) {
    const std::optional<Vec2> centre = singularPoint(element);
    return centre && lengthSquared(p - *centre) < radius * radius;
}

bool within(const PeriodicRow& row, Vec2 p, double radius) {
    return lengthSquared(p - nearestCopy(row, p)) < radius * radius;
}

} // namespace

Complex FlowField::complexVelocity(Vec2 p) const {
    Complex w{};
    for (const FlowElement& element : elements) {
        if (!within(element, p, singularRadius)) w += pfl::complexVelocity(element, p);
    }
    for (const PeriodicRow& row : rows) {
        if (!within(row, p, singularRadius)) w += pfl::complexVelocity(row, p);
    }
    return w;
}

Complex FlowField::complexPotential(Vec2 p) const {
    Complex f{};
    for (const FlowElement& element : elements) {
        if (!within(element, p, singularRadius)) f += pfl::complexPotential(element, p);
    }
    for (const PeriodicRow& row : rows) {
        if (!within(row, p, singularRadius)) f += pfl::complexPotential(row, p);
    }
    return f;
}

bool FlowField::inDomain(Vec2 p) const {
    return std::ranges::all_of(elements, [p](const FlowElement& e) { return inElementDomain(e, p); });
}

bool FlowField::nearSingularity(Vec2 p, double radius) const {
    return std::ranges::any_of(elements, [p, radius](const FlowElement& e) { return within(e, p, radius); }) ||
           std::ranges::any_of(rows, [p, radius](const PeriodicRow& r) { return within(r, p, radius); });
}

std::vector<Vec2> FlowField::singularPoints(const Rect& region) const {
    std::vector<Vec2> points;
    for (const FlowElement& element : elements) {
        const std::optional<Vec2> centre = singularPoint(element);
        if (centre && region.contains(*centre)) points.push_back(*centre);
    }
    for (const PeriodicRow& row : rows) {
        for (const FlowElement& copy : rowCopies(row, region)) points.push_back(*singularPoint(copy));
    }
    return points;
}

double typicalSpeed(const FlowField& field, const Rect& region, int gridSize) {
    const Vec2 size = region.size();
    const double guard = 0.02 * std::max(size.x, size.y);
    std::vector<double> speeds;
    speeds.reserve(static_cast<std::size_t>(gridSize) * static_cast<std::size_t>(gridSize));
    for (int j = 0; j < gridSize; ++j) {
        for (int i = 0; i < gridSize; ++i) {
            const Vec2 p = region.min + Vec2{size.x * (i + 0.5) / gridSize, size.y * (j + 0.5) / gridSize};
            if (!field.inDomain(p) || field.nearSingularity(p, guard)) continue;
            const double s = length(field.velocity(p));
            if (std::isfinite(s)) speeds.push_back(s);
        }
    }
    if (speeds.empty()) return 0.0;
    auto middle = speeds.begin() + static_cast<std::ptrdiff_t>(speeds.size() / 2);
    std::ranges::nth_element(speeds, middle);
    return *middle;
}

ScalarGrid sampleSpeed(const FlowField& field, const Rect& bounds, int columns, int rows, double singularRadius) {
    ScalarGrid grid{.bounds = bounds, .columns = std::max(columns, 0), .rows = std::max(rows, 0)};
    grid.values.reserve(static_cast<std::size_t>(grid.columns) * static_cast<std::size_t>(grid.rows));
    const Vec2 cell{bounds.size().x / std::max(columns, 1), bounds.size().y / std::max(rows, 1)};
    for (int j = 0; j < grid.rows; ++j) {
        for (int i = 0; i < grid.columns; ++i) {
            const Vec2 p = bounds.min + Vec2{cell.x * (i + 0.5), cell.y * (j + 0.5)};
            if (!field.inDomain(p)) {
                grid.values.push_back(std::numeric_limits<double>::quiet_NaN());
            } else if (field.nearSingularity(p, singularRadius)) {
                grid.values.push_back(std::numeric_limits<double>::infinity());
            } else {
                grid.values.push_back(length(field.velocity(p)));
            }
        }
    }
    return grid;
}

} // namespace pfl
