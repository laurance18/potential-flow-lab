#include "core/flow_field.hpp"

#include "core/parallel.hpp"

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
    for (const SourcePanels& body : panelBodies) w += body.complexVelocity(p);
    for (const VortexSheet& sheet : vortexSheets) w += sheet.complexVelocity(p);
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
    for (const SourcePanels& body : panelBodies) f += body.complexPotential(p);
    for (const VortexSheet& sheet : vortexSheets) f += sheet.complexPotential(p);
    return f;
}

bool FlowField::inDomain(Vec2 p) const {
    return std::ranges::all_of(elements, [p](const FlowElement& e) { return inElementDomain(e, p); }) &&
           panelBodyAt(p) == nullptr;
}

const SourcePanels* FlowField::panelBodyAt(Vec2 p) const {
    for (const SourcePanels& body : panelBodies) {
        if (body.contains(p)) return &body;
    }
    return nullptr;
}

bool FlowField::nearSingularity(Vec2 p, double radius) const {
    return std::ranges::any_of(elements, [p, radius](const FlowElement& e) { return within(e, p, radius); }) ||
           std::ranges::any_of(rows, [p, radius](const PeriodicRow& r) { return within(r, p, radius); }) ||
           std::ranges::any_of(vortexSheets, [p, radius](const VortexSheet& s) {
               return lengthSquared(p - s.leadingEdge()) < radius * radius;
           });
}

bool FlowField::crossesSheet(Vec2 a, Vec2 b) const {
    return std::ranges::any_of(vortexSheets, [a, b](const VortexSheet& s) { return s.crossedBy(a, b); });
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
    const std::size_t n = static_cast<std::size_t>(std::max(gridSize, 0));
    std::vector<double> samples(n * n, std::numeric_limits<double>::quiet_NaN()); // NaN: skipped
    parallelFor(n, 4, [&](std::size_t begin, std::size_t end) {
        for (std::size_t j = begin; j < end; ++j) {
            for (std::size_t i = 0; i < n; ++i) {
                const Vec2 p = region.min + Vec2{size.x * (static_cast<double>(i) + 0.5) / gridSize,
                                                 size.y * (static_cast<double>(j) + 0.5) / gridSize};
                if (!field.inDomain(p) || field.nearSingularity(p, guard)) continue;
                samples[j * n + i] = length(field.velocity(p));
            }
        }
    });
    std::vector<double> speeds;
    speeds.reserve(samples.size());
    for (const double s : samples) {
        if (std::isfinite(s)) speeds.push_back(s);
    }
    if (speeds.empty()) return 0.0;
    auto middle = speeds.begin() + static_cast<std::ptrdiff_t>(speeds.size() / 2);
    std::ranges::nth_element(speeds, middle);
    return *middle;
}

ScalarGrid sampleSpeed(const FlowField& field, const Rect& bounds, int columns, int rows, double singularRadius) {
    ScalarGrid grid{.bounds = bounds, .columns = std::max(columns, 0), .rows = std::max(rows, 0)};
    const std::size_t width = static_cast<std::size_t>(grid.columns);
    grid.values.resize(width * static_cast<std::size_t>(grid.rows));
    const Vec2 cell{bounds.size().x / std::max(columns, 1), bounds.size().y / std::max(rows, 1)};
    // Row by row on every core.
    parallelFor(static_cast<std::size_t>(grid.rows), 4, [&](std::size_t begin, std::size_t end) {
        for (std::size_t j = begin; j < end; ++j) {
            for (std::size_t i = 0; i < width; ++i) {
                const Vec2 p = bounds.min + Vec2{cell.x * (static_cast<double>(i) + 0.5), cell.y * (static_cast<double>(j) + 0.5)};
                double& value = grid.values[j * width + i];
                if (!field.inDomain(p)) {
                    value = std::numeric_limits<double>::quiet_NaN();
                } else if (field.nearSingularity(p, singularRadius)) {
                    value = std::numeric_limits<double>::infinity();
                } else {
                    value = length(field.velocity(p));
                }
            }
        }
    });
    return grid;
}

} // namespace pfl
