// flow_field.hpp: superposition of flow elements.
//
// Holds the list of active elements and evaluates the total complex velocity and potential by
// summing each element's W and F. This is the single source of truth that the UI edits and the
// renderers and streamline integrator read.
//
// Singularities: an element's term is skipped within singularRadius of its singular point (see
// singularPoint()), so the total stays finite there. Callers that should stop or respawn near a
// singularity (streamlines, particles) ask nearSingularity() with their own radius.
//
// Domain: sector flows are defined only inside their wedge, so the admissible domain is the
// intersection of all sector wedges. Outside it (inDomain() false) the field is undefined and must
// be masked. Other elements are superposed unrestricted and do not keep the sector walls as
// streamlines. Wall images are added as ordinary elements by buildFlowModel() (scene.hpp), or, in a
// channel, as periodic rows that take the place of the elements they contain.
#pragma once

#include "core/flow_element.hpp"
#include "core/vec2.hpp"

#include <vector>

namespace pfl {

struct FlowField {
    std::vector<FlowElement> elements;
    std::vector<PeriodicRow> rows;
    double singularRadius = 1e-9;

    bool empty() const { return elements.empty() && rows.empty(); }

    Complex complexVelocity(Vec2 p) const;
    Complex complexPotential(Vec2 p) const;
    Vec2 velocity(Vec2 p) const { return velocityFromW(complexVelocity(p)); }

    bool inDomain(Vec2 p) const;
    bool nearSingularity(Vec2 p, double radius) const;

    // Every singular point in `region`, those of the periodic rows included.
    std::vector<Vec2> singularPoints(const Rect& region) const;
};

// Median speed over a grid of points in `region` (skipping points outside the domain or close to a
// singularity): a scale for colour maps, animation speed and stagnation thresholds. 0 if no flow.
double typicalSpeed(const FlowField& field, const Rect& region, int gridSize = 24);

// A scalar sampled at the cell centres of a columns x rows grid over `bounds`; row 0 is at
// bounds.min.y.
struct ScalarGrid {
    Rect bounds;
    int columns = 0;
    int rows = 0;
    std::vector<double> values; // row by row

    double at(int i, int j) const {
        return values[static_cast<std::size_t>(j) * static_cast<std::size_t>(columns) + static_cast<std::size_t>(i)];
    }
};

// |V| on a grid: NaN outside the flow domain, +infinity within singularRadius of a singular point.
ScalarGrid sampleSpeed(const FlowField& field, const Rect& bounds, int columns, int rows, double singularRadius);

} // namespace pfl
