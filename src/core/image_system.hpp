// image_system.hpp: method of images for straight symmetry walls.
//
// A symmetry wall is the infinite straight line through `point` at angle phi. Reflecting a point
// element across it gives an image whose normal velocity on the line cancels the original's:
//     source / sink: same strength Lambda      vortex: circulation -Gamma
//     doublet: same mu, axis angle beta' = 2 phi - beta
// (each image is F_image(z) = conj(F(reflect(z))), so the line is a streamline of the pair).
//
// Several walls: images are reflected again across every wall until no new copies appear. That
// closes exactly for one wall or walls meeting at pi/k; otherwise (e.g. two parallel walls, which
// need an infinite row of images) the series is truncated and ImageSystem::exact is false.
//
// Uniform flow is not mirrored, since its image would double the component along the wall. It is
// instead reduced to the part every wall allows: the component along the walls when they are all
// parallel, nothing when two walls cross. Sector flows are not mirrored either; a wall only keeps
// zero normal velocity for a sector if it lies along one of the sector's own walls.
#pragma once

#include "core/flow_element.hpp"
#include "core/vec2.hpp"

#include <cmath>
#include <optional>
#include <span>
#include <vector>

namespace pfl {

struct SymmetryWall {
    Vec2 point;         // any point on the wall
    double angle = 0.0; // phi: direction of the line, counter-clockwise from +x

    Vec2 direction() const { return {std::cos(angle), std::sin(angle)}; }
    Vec2 normal() const { return perp(direction()); }

    friend bool operator==(const SymmetryWall&, const SymmetryWall&) = default;
};

Vec2 reflect(Vec2 p, const SymmetryWall& wall);

// Image of a source, vortex or doublet; nullopt for uniform and sector flows (see above).
std::optional<FlowElement> mirror(const FlowElement& element, const SymmetryWall& wall);

struct ImageSystem {
    std::vector<FlowElement> images; // mirrored copies only, not the originals
    bool exact = true;               // false if the series of images had to be truncated
};

ImageSystem buildImageSystem(std::span<const FlowElement> elements, std::span<const SymmetryWall> walls);

// The uniform flow left after imposing the walls (see above). Unchanged when there are no walls.
UniformFlow constrainToWalls(const UniformFlow& flow, std::span<const SymmetryWall> walls);

// True if two elements of the same type have the same parameters (to a relative 1e-9).
bool sameElement(const FlowElement& a, const FlowElement& b);

} // namespace pfl
