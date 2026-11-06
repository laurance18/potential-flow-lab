// scene.hpp: what the user builds, edits and saves: the element list (flows and symmetry walls),
// the dye rakes and the plot settings, and the flow model derived from it.
//
// buildFlowModel() turns a scene into the field that is actually evaluated: every flow element,
// the wall images of its sources, vortices and doublets, and its uniform flows reduced to what the
// walls allow (see image_system.hpp). The UI shows warnings for configurations that are legal but
// easy to misread (coincident singularities, elements on a wall, truncated images, ...).
#pragma once

#include "core/dye.hpp"
#include "core/flow_element.hpp"
#include "core/flow_field.hpp"
#include "core/image_system.hpp"
#include "core/vec2.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <variant>
#include <vector>

namespace pfl {

using SceneElement = std::variant<UniformFlow, SourceSink, Vortex, Doublet, SectorFlow, SymmetryWall>;

// One per SceneElement alternative, in the same order.
enum class ElementKind { Uniform, SourceSink, Vortex, Doublet, Sector, Wall };
inline constexpr std::array<ElementKind, 6> k_elementKinds = {ElementKind::Uniform, ElementKind::SourceSink,
                                                               ElementKind::Vortex,  ElementKind::Doublet,
                                                               ElementKind::Sector,  ElementKind::Wall};

inline ElementKind kindOf(const SceneElement& element) { return static_cast<ElementKind>(element.index()); }
const char* kindName(ElementKind kind);

// A new element with sensible default parameters, placed at `position` where that applies.
SceneElement defaultElement(ElementKind kind, Vec2 position);

// One-line parameter summary in UTF-8, e.g. "Λ = 6.283 at (1, 0)". Angles in degrees.
std::string describe(const SceneElement& element);

// Why the parameters are invalid, or nullopt if they are fine.
std::optional<std::string> sceneElementError(const SceneElement& element);

// The point an element is drawn at and dragged by: the position of a source, vortex or doublet, the
// vertex of a sector, the reference point of a wall. nullopt for a uniform flow, which has none.
std::optional<Vec2> anchorPoint(const SceneElement& element);

// `element` translated so that its anchor point is at `position` (a uniform flow is unchanged).
SceneElement movedTo(SceneElement element, Vec2 position);

// Filled contour bands drawn under the streamlines.
enum class Background { None, Speed, PressureCoefficient };

struct PlotSettings {
    Rect view{{-4.0, -2.5}, {4.0, 2.5}};
    int streamlineDensity = 24; // streamlines across the view height
    bool showArrows = true;
    double arrowSpacingPx = 110.0; // screen distance between direction arrows along a streamline
    bool showImages = true;        // draw the wall images of elements
    bool showEquipotentials = false;
    Background background = Background::None;
    bool showSurfacePressure = false; // body surfaces coloured by Cp, pressure arrows and the Cp panel
    bool showStagnationPoints = true;
    bool showParticles = true;
    bool particleTrails = true;
};

struct Scene {
    std::vector<SceneElement> elements;
    std::vector<DyeRake> rakes; // flow visualization only; they do not change the flow
    PlotSettings plot;
};

struct FlowModel {
    FlowField field;                  // what is evaluated: flows (uniform ones constrained) + images
    std::vector<FlowElement> images;  // the wall images on their own, for drawing (a channel has
                                      // none: its images are in field.rows, see rowCopies())
    std::vector<SymmetryWall> walls;
    bool imagesExact = true;
    Vec2 freeStream;                  // sum of the constrained uniform flows
    double freeStreamSpeed = 0.0;     // |freeStream|; the Cp reference
};

FlowModel buildFlowModel(const Scene& scene);

// Human-readable notes about the scene (element numbers are 1-based, as in the Elements menu).
std::vector<std::string> sceneWarnings(const Scene& scene, const FlowModel& model);

// Numerical check that a wall is a streamline: velocity sampled along the part of the wall inside
// `region`, skipping points near singularities and outside the flow domain.
struct WallCheck {
    double maxNormalSpeed = 0.0; // largest |v . n|
    double maxSpeed = 0.0;       // largest |v|, for scale
    int samples = 0;
};
WallCheck checkWall(const FlowField& field, const SymmetryWall& wall, const Rect& region, int samples = 64);

} // namespace pfl
