// streamline.hpp: evenly spaced streamlines with direction arrows.
//
// Streamlines are integrated with classic RK4 along the unit direction field v/|v| (arc-length
// steps, so slow and fast regions get the same resolution). Seeding follows Jobard & Lefer (1997):
// new seeds are placed one separation distance beside existing lines, and a line stops when it
// comes within half a separation of another one, so the picture has an even density.
//
// A line also stops when it leaves the bounds or the flow domain, comes within singularRadius of a
// singular point, reaches a stagnation point (speed below 1e-6 of the typical speed, or the
// direction reversing within one step), closes on itself, or reaches maxPointsPerLine.
//
// With `orthogonal` set, the same tracer follows the rotated field perp(v)/|v| instead. Those curves
// cross every streamline at right angles, and since grad phi = v is normal to them, phi is constant
// along each: they are the equipotential lines. Tracing them (rather than contouring phi) avoids the
// branch cut of phi around a vortex, where phi is multivalued.
//
// Live changes (a sweep, a slider, a drag, a pan) recompute the lines every frame. Seeded from
// scratch, a tiny change in the field changes which seeds find room, and that cascades through the
// picture: lines pop in and out from frame to frame. So the set records each line's seed, and the
// next computation can start from those seeds first, in the same order. A line through the same
// seed in a slightly different field moves only slightly. Seeds then need a little less free space,
// and a kept seed crowded by the lines before it slides along its own streamline to where there is
// room, so the same line stays; lines appear or vanish mostly where the flow stretches or really
// changes shape. Each seed is recorded near the middle of its line, away from the ends. Arrows are placed at fixed arc-length steps from the line's seed, with
// an offset kept with the seed, so they move with their line instead of jumping along it.
#pragma once

#include "core/flow_field.hpp"
#include "core/vec2.hpp"

#include <span>
#include <vector>

namespace pfl {

struct StreamlineSettings {
    Rect bounds;                  // trace only inside this rectangle
    double separation = 0.25;     // distance between neighbouring streamlines
    double singularRadius = 0.02; // stop this close to a singular point
    double arrowSpacing = 1.0;    // arc length between direction arrows on a line; <= 0 for none
    int maxPointsPerLine = 5000;
    bool orthogonal = false;      // trace equipotential lines instead of streamlines
};

struct Arrow {
    Vec2 position;
    Vec2 direction; // unit vector along the flow
};

struct StreamlineSeed {
    Vec2 position;
    double arrowPhase = 0.0; // offset of the line's arrows from the seed, a fraction of arrowSpacing
};

struct StreamlineSet {
    std::vector<std::vector<Vec2>> lines; // each ordered in the flow direction
    std::vector<Arrow> arrows;
    std::vector<StreamlineSeed> seeds; // seeds[i] started lines[i]
};

// `previousSeeds`: the seeds of the lines shown before (with the same separation), tried first so
// the picture changes smoothly; empty to seed from scratch.
StreamlineSet computeStreamlines(const FlowField& field, const StreamlineSettings& settings,
                                 std::span<const StreamlineSeed> previousSeeds = {});

} // namespace pfl
