#include "core/scene.hpp"

#include <cmath>
#include <format>
#include <numbers>

namespace pfl {
namespace {

using std::numbers::pi;
constexpr double k_degrees = 180.0 / pi;
constexpr double k_coincidenceTolerance = 1e-9;

template <typename... Ts>
struct Overloaded : Ts... {
    using Ts::operator()...;
};

std::string formatPoint(Vec2 p) { return std::format("({:.4g}, {:.4g})", p.x, p.y); }

bool coincide(Vec2 a, Vec2 b) {
    const double scale = 1.0 + std::max(length(a), length(b));
    return length(a - b) <= k_coincidenceTolerance * scale;
}

double distanceToWall(Vec2 p, const SymmetryWall& wall) { return std::abs(dot(p - wall.point, wall.normal())); }

// The flow element of a scene element, or nullopt for a wall.
std::optional<FlowElement> asFlowElement(const SceneElement& element) {
    return std::visit(Overloaded{
                          [](const SymmetryWall&) -> std::optional<FlowElement> { return std::nullopt; },
                          [](const auto& e) -> std::optional<FlowElement> { return FlowElement{e}; },
                      },
                      element);
}

} // namespace

const char* kindName(ElementKind kind) {
    switch (kind) {
    case ElementKind::Uniform: return "Uniform flow";
    case ElementKind::SourceSink: return "Source / sink";
    case ElementKind::Vortex: return "Point vortex";
    case ElementKind::Doublet: return "Doublet";
    case ElementKind::Sector: return "Sector (wedge) flow";
    case ElementKind::Wall: return "Symmetry wall";
    }
    return "?";
}

SceneElement defaultElement(ElementKind kind, Vec2 position) {
    switch (kind) {
    case ElementKind::Uniform: return UniformFlow{.speed = 1.0, .angle = 0.0};
    case ElementKind::SourceSink: return SourceSink{.strength = 2.0 * pi, .position = position};
    case ElementKind::Vortex: return Vortex{.circulation = 2.0 * pi, .position = position};
    case ElementKind::Doublet: return Doublet{.strength = 1.0, .position = position, .angle = 0.0};
    case ElementKind::Sector: return SectorFlow{.strength = 1.0, .vertex = position, .wall1 = 0.0, .wall2 = pi / 2.0};
    case ElementKind::Wall: return SymmetryWall{.point = position, .angle = 0.0};
    }
    return UniformFlow{};
}

std::string describe(const SceneElement& element) {
    return std::visit(
        Overloaded{
            [](const UniformFlow& e) { return std::format("U = {:.4g}, α = {:.4g}°", e.speed, e.angle * k_degrees); },
            [](const SourceSink& e) {
                return std::format("Λ = {:.4g} at {}{}", e.strength, formatPoint(e.position),
                                   e.strength < 0.0 ? " (sink)" : "");
            },
            [](const Vortex& e) { return std::format("Γ = {:.4g} at {}", e.circulation, formatPoint(e.position)); },
            [](const Doublet& e) {
                return std::format("μ = {:.4g} at {}, β = {:.4g}°", e.strength, formatPoint(e.position),
                                   e.angle * k_degrees);
            },
            [](const SectorFlow& e) {
                return std::format("V = {:.4g} at {}, θ₁ = {:.4g}°, θ₂ = {:.4g}° (n = {:.4g})", e.strength,
                                   formatPoint(e.vertex), e.wall1 * k_degrees, e.wall2 * k_degrees, e.exponent());
            },
            [](const SymmetryWall& e) {
                return std::format("through {} at {:.4g}°", formatPoint(e.point), e.angle * k_degrees);
            },
        },
        element);
}

std::optional<std::string> sceneElementError(const SceneElement& element) {
    if (const auto* wall = std::get_if<SymmetryWall>(&element)) {
        if (!std::isfinite(wall->point.x) || !std::isfinite(wall->point.y) || !std::isfinite(wall->angle)) {
            return "All parameters must be finite numbers.";
        }
        return std::nullopt;
    }
    return validationError(*asFlowElement(element));
}

std::optional<Vec2> anchorPoint(const SceneElement& element) {
    return std::visit(Overloaded{
                          [](const UniformFlow&) -> std::optional<Vec2> { return std::nullopt; },
                          [](const SectorFlow& e) -> std::optional<Vec2> { return e.vertex; },
                          [](const SymmetryWall& e) -> std::optional<Vec2> { return e.point; },
                          [](const auto& e) -> std::optional<Vec2> { return e.position; },
                      },
                      element);
}

SceneElement movedTo(SceneElement element, Vec2 position) {
    std::visit(Overloaded{
                   [](UniformFlow&) {},
                   [position](SectorFlow& e) { e.vertex = position; },
                   [position](SymmetryWall& e) { e.point = position; },
                   [position](auto& e) { e.position = position; },
               },
               element);
    return element;
}

FlowModel buildFlowModel(const Scene& scene) {
    FlowModel model;
    std::vector<FlowElement> flows;
    for (const SceneElement& element : scene.elements) {
        if (const auto* wall = std::get_if<SymmetryWall>(&element)) {
            model.walls.push_back(*wall);
        } else {
            flows.push_back(*asFlowElement(element));
        }
    }

    ImageSystem images = buildImageSystem(flows, model.walls);
    model.imagesExact = images.exact;
    model.images = std::move(images.images);

    Complex freeStream{};
    for (const FlowElement& flow : flows) {
        if (const auto* uniform = std::get_if<UniformFlow>(&flow)) {
            const UniformFlow constrained = constrainToWalls(*uniform, model.walls);
            freeStream += complexVelocity(constrained, {});
            model.field.elements.push_back(constrained);
        } else if (images.rows.empty() || !mirror(flow, model.walls.front())) {
            model.field.elements.push_back(flow); // otherwise it is in its row
        }
    }
    model.field.rows = std::move(images.rows);
    model.field.elements.insert(model.field.elements.end(), model.images.begin(), model.images.end());
    model.freeStream = velocityFromW(freeStream);
    model.freeStreamSpeed = std::abs(freeStream);
    return model;
}

std::vector<std::string> sceneWarnings(const Scene& scene, const FlowModel& model) {
    std::vector<std::string> warnings;
    const auto& elements = scene.elements;

    // Singular points (and sector vertices) of the original elements, with their 1-based numbers.
    struct Centre {
        std::size_t number;
        Vec2 point;
    };
    std::vector<Centre> centres;
    std::vector<std::pair<std::size_t, SymmetryWall>> walls;
    for (std::size_t i = 0; i < elements.size(); ++i) {
        const SceneElement& e = elements[i];
        if (const auto* wall = std::get_if<SymmetryWall>(&e)) {
            walls.emplace_back(i + 1, *wall);
        } else if (const auto* sector = std::get_if<SectorFlow>(&e)) {
            centres.push_back({i + 1, sector->vertex});
        } else if (const std::optional<Vec2> point = singularPoint(*asFlowElement(e))) {
            centres.push_back({i + 1, *point});
        }
    }

    for (std::size_t a = 0; a < centres.size(); ++a) {
        for (std::size_t b = a + 1; b < centres.size(); ++b) {
            if (!coincide(centres[a].point, centres[b].point)) continue;
            warnings.push_back(std::format("Elements #{} and #{} share the point {}; their fields add there "
                                           "(e.g. source + vortex is a spiral vortex).",
                                           centres[a].number, centres[b].number, formatPoint(centres[a].point)));
        }
    }

    for (const auto& [wallNumber, wall] : walls) {
        for (const Centre& c : centres) {
            if (distanceToWall(c.point, wall) > k_coincidenceTolerance * (1.0 + length(c.point))) continue;
            if (std::holds_alternative<SectorFlow>(elements[c.number - 1])) continue; // handled below
            warnings.push_back(std::format("Element #{} lies on wall #{} and coincides with its own image "
                                           "(a vortex there is cancelled by it).",
                                           c.number, wallNumber));
        }
    }

    for (std::size_t i = 0; i < elements.size(); ++i) {
        if (const auto* uniform = std::get_if<UniformFlow>(&elements[i]); uniform != nullptr && !walls.empty()) {
            const UniformFlow constrained = constrainToWalls(*uniform, model.walls);
            if (std::abs(constrained.speed - uniform->speed) > 1e-9 * (1.0 + uniform->speed)) {
                warnings.push_back(std::format("Uniform flow #{}: only the component along the walls is kept "
                                               "(U = {:.4g} instead of {:.4g}); an infinite wall admits no flow "
                                               "across it.",
                                               i + 1, constrained.speed, uniform->speed));
            }
        }
        if (const auto* sector = std::get_if<SectorFlow>(&elements[i])) {
            for (const auto& [wallNumber, wall] : walls) {
                // A symmetry wall is a streamline of the sector term only if it runs along one of the
                // sector's own walls (through the vertex, parallel to theta1 or theta2).
                const bool throughVertex = distanceToWall(sector->vertex, wall) <= 1e-9 * (1.0 + length(sector->vertex));
                const bool alongWall = std::abs(std::sin(wall.angle - sector->wall1)) <= 1e-9 ||
                                       std::abs(std::sin(wall.angle - sector->wall2)) <= 1e-9;
                if (throughVertex && alongWall) continue;
                warnings.push_back(std::format("Sector #{} is not mirrored by wall #{}, so that wall is not a "
                                               "streamline of the sector flow.",
                                               i + 1, wallNumber));
            }
        }
    }

    const auto sectorCount = std::ranges::count_if(elements, [](const SceneElement& e) {
        return std::holds_alternative<SectorFlow>(e);
    });
    if (sectorCount > 1) {
        warnings.push_back("Several sector flows: the flow domain is the intersection of their wedges.");
    }
    if (!model.field.rows.empty()) {
        const Vec2 normal = model.walls.front().normal();
        const auto offset = [&](Vec2 p) { return dot(p - model.walls.front().point, normal); };
        double low = 0.0;
        double high = 0.0;
        for (const SymmetryWall& wall : model.walls) {
            low = std::min(low, offset(wall.point));
            high = std::max(high, offset(wall.point));
        }
        for (const Centre& c : centres) {
            if (std::holds_alternative<SectorFlow>(elements[c.number - 1])) continue;
            const double tolerance = k_coincidenceTolerance * (1.0 + length(c.point));
            if (offset(c.point) >= low - tolerance && offset(c.point) <= high + tolerance) continue;
            warnings.push_back(std::format("Element #{} lies outside the channel, so one of its images lies inside it "
                                           "and acts as a real element there.",
                                           c.number));
        }
    }
    if (!model.imagesExact) {
        warnings.push_back("The walls do not form a closed image system (e.g. three parallel walls, or walls not "
                           "meeting at 180°/k); the images were truncated, so the wall condition is approximate.");
    }
    return warnings;
}

WallCheck checkWall(const FlowField& field, const SymmetryWall& wall, const Rect& region, int samples) {
    WallCheck check;
    const Vec2 t = wall.direction();
    double t0 = -1e300;
    double t1 = 1e300;
    if (samples < 1 || !clipLine(wall.point, t, region, t0, t1)) return check;

    const Vec2 n = wall.normal();
    const double guard = 1e-3 * length(region.size());
    for (int i = 0; i < samples; ++i) {
        const double s = t0 + (t1 - t0) * (i + 0.5) / samples;
        const Vec2 p = wall.point + t * s;
        if (!field.inDomain(p) || field.nearSingularity(p, guard)) continue;
        const Vec2 v = field.velocity(p);
        check.maxNormalSpeed = std::max(check.maxNormalSpeed, std::abs(dot(v, n)));
        check.maxSpeed = std::max(check.maxSpeed, length(v));
        ++check.samples;
    }
    return check;
}

} // namespace pfl
