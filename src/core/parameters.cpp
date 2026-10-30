#include "core/parameters.hpp"

#include <array>
#include <format>

namespace pfl {
namespace {

template <typename T, double T::*Field>
constexpr ElementParameter field(const char* label, const char* help, ParameterKind kind) {
    return {label, help, kind, [](const SceneElement& e) { return std::get<T>(e).*Field; },
            [](SceneElement& e, double v) { std::get<T>(e).*Field = v; }};
}

template <typename T, Vec2 T::*Point, double Vec2::*Coord>
constexpr ElementParameter coordinate(const char* label, const char* help, ParameterKind kind) {
    return {label, help, kind, [](const SceneElement& e) { return std::get<T>(e).*Point.*Coord; },
            [](SceneElement& e, double v) { std::get<T>(e).*Point.*Coord = v; }};
}

template <typename T, Vec2 T::*Point>
constexpr ElementParameter xOf(const char* help) {
    return coordinate<T, Point, &Vec2::x>("x₀", help, ParameterKind::X);
}

template <typename T, Vec2 T::*Point>
constexpr ElementParameter yOf(const char* help) {
    return coordinate<T, Point, &Vec2::y>("y₀", help, ParameterKind::Y);
}

constexpr std::array k_uniform = {
    field<UniformFlow, &UniformFlow::speed>("U∞", "Free-stream speed (not negative)", ParameterKind::Speed),
    field<UniformFlow, &UniformFlow::angle>("α", "Direction, counter-clockwise from +x", ParameterKind::Angle),
};
constexpr std::array k_source = {
    field<SourceSink, &SourceSink::strength>("Λ", "Strength: > 0 source, < 0 sink", ParameterKind::Strength),
    xOf<SourceSink, &SourceSink::position>("x coordinate of the source / sink"),
    yOf<SourceSink, &SourceSink::position>("y coordinate of the source / sink"),
};
constexpr std::array k_vortex = {
    field<Vortex, &Vortex::circulation>("Γ", "Circulation: > 0 clockwise, < 0 counter-clockwise", ParameterKind::Strength),
    xOf<Vortex, &Vortex::position>("x coordinate of the vortex"),
    yOf<Vortex, &Vortex::position>("y coordinate of the vortex"),
};
constexpr std::array k_doublet = {
    field<Doublet, &Doublet::strength>("μ", "Strength, as in F = μ e^(iβ) / (z − z₀)", ParameterKind::Strength),
    xOf<Doublet, &Doublet::position>("x coordinate of the doublet"),
    yOf<Doublet, &Doublet::position>("y coordinate of the doublet"),
    field<Doublet, &Doublet::angle>("β", "Axis angle, counter-clockwise", ParameterKind::Angle),
};
constexpr std::array k_sector = {
    field<SectorFlow, &SectorFlow::strength>("V", "Scaling constant in F = V (ζ e^(−iθ₁))ⁿ", ParameterKind::Strength),
    xOf<SectorFlow, &SectorFlow::vertex>("x coordinate of the wedge vertex"),
    yOf<SectorFlow, &SectorFlow::vertex>("y coordinate of the wedge vertex"),
    field<SectorFlow, &SectorFlow::wall1>("θ₁", "Angle of the first wall; 0 < θ₂ − θ₁ ≤ 360°", ParameterKind::Angle),
    field<SectorFlow, &SectorFlow::wall2>("θ₂", "Angle of the second wall; 0 < θ₂ − θ₁ ≤ 360°", ParameterKind::Angle),
};
constexpr std::array k_wall = {
    xOf<SymmetryWall, &SymmetryWall::point>("x coordinate of a point on the wall"),
    yOf<SymmetryWall, &SymmetryWall::point>("y coordinate of a point on the wall"),
    field<SymmetryWall, &SymmetryWall::angle>("φ", "Direction of the wall line, counter-clockwise from +x",
                                              ParameterKind::Angle),
};

} // namespace

std::span<const ElementParameter> elementParameters(ElementKind kind) {
    switch (kind) {
    case ElementKind::Uniform: return k_uniform;
    case ElementKind::SourceSink: return k_source;
    case ElementKind::Vortex: return k_vortex;
    case ElementKind::Doublet: return k_doublet;
    case ElementKind::Sector: return k_sector;
    case ElementKind::Wall: return k_wall;
    }
    return {};
}

const ElementParameter* findParameter(const Scene& scene, ParameterRef ref) {
    if (ref.element >= scene.elements.size()) return nullptr;
    const std::span<const ElementParameter> parameters = elementParameters(kindOf(scene.elements[ref.element]));
    return ref.parameter < parameters.size() ? &parameters[ref.parameter] : nullptr;
}

std::string parameterName(const Scene& scene, ParameterRef ref) {
    const ElementParameter* parameter = findParameter(scene, ref);
    if (parameter == nullptr) return "-";
    return std::format("#{} {}: {}", ref.element + 1, kindName(kindOf(scene.elements[ref.element])), parameter->label);
}

std::optional<SceneElement> withParameter(const Scene& scene, ParameterRef ref, double value) {
    const ElementParameter* parameter = findParameter(scene, ref);
    if (parameter == nullptr) return std::nullopt;
    SceneElement element = scene.elements[ref.element];
    parameter->set(element, value);
    if (sceneElementError(element)) return std::nullopt;
    return element;
}

} // namespace pfl
