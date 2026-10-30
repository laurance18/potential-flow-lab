// parameters.hpp: the numeric parameters of each kind of scene element, read and written one at a
// time, for the live sliders and the parameter sweep.
//
// Each kind lists its parameters in the order the parameter dialog shows them. Angles are read and
// written in radians (the UI shows degrees).
#pragma once

#include "core/scene.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <string>

namespace pfl {

// What a parameter measures; the UI picks the slider range and units from it.
enum class ParameterKind { Speed, Strength, Angle, X, Y };

struct ElementParameter {
    const char* label; // UTF-8 symbol, e.g. "Γ", "x₀", "α"
    const char* help;
    ParameterKind kind;
    double (*get)(const SceneElement&);
    void (*set)(SceneElement&, double);
};

std::span<const ElementParameter> elementParameters(ElementKind kind);

// One parameter of one element of a scene.
struct ParameterRef {
    std::size_t element = 0;
    std::size_t parameter = 0;

    bool operator==(const ParameterRef&) const = default;
};

// The parameter `ref` points to, or nullptr if the scene has no such element or parameter.
const ElementParameter* findParameter(const Scene& scene, ParameterRef ref);

// E.g. "#3 Point vortex: Γ" (element numbers are 1-based, as in the Elements menu).
std::string parameterName(const Scene& scene, ParameterRef ref);

// The element `ref` points to with the parameter set to `value`, or nullopt if `ref` is invalid or
// the value makes the element invalid (see sceneElementError()).
std::optional<SceneElement> withParameter(const Scene& scene, ParameterRef ref, double value);

} // namespace pfl
