#include "core/sweep.hpp"

#include "core/surface_pressure.hpp"

#include <algorithm>
#include <limits>

namespace pfl {

std::vector<double> sweepValues(double from, double to, int steps) {
    if (steps < 2) return {from};
    std::vector<double> values(static_cast<std::size_t>(steps));
    for (int i = 0; i < steps; ++i) {
        const double t = static_cast<double>(i) / (steps - 1);
        values[static_cast<std::size_t>(i)] = i == steps - 1 ? to : from + t * (to - from);
    }
    return values;
}

std::optional<SweepSample> sampleSweep(const Scene& scene, ParameterRef ref, double value, const Rect& view) {
    std::optional<SceneElement> element = withParameter(scene, ref, value);
    if (!element) return std::nullopt;
    Scene swept = scene;
    swept.elements[ref.element] = std::move(*element);
    const FlowModel model = buildFlowModel(swept);

    SweepSample sample;
    sample.value = value;
    if (const double speed = typicalSpeed(model.field, view); speed > 0.0) {
        sample.stagnationPoints = findStagnationPoints(model.field, view, speed);
    }

    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    for (const BodySurface& surface : findBodySurfaces(model, view)) {
        SweepBody body;
        body.closed = surface.closed;
        body.centroid = surface.centroid;
        body.cpMin = surface.cp.empty() ? nan : std::ranges::min(surface.cp);
        body.chord = body.circulation = body.lift = body.drag = body.exactLift = nan;
        if (surface.closed) {
            const BodyForces forces = bodyForces(surface, model);
            body.chord = forces.chord;
            body.circulation = forces.circulation;
            body.lift = forces.lift;
            body.drag = forces.drag;
            if (forces.isolated) body.exactLift = 2.0 * forces.circulation / model.freeStreamSpeed;
        }
        sample.bodies.push_back(body);
    }
    return sample;
}

} // namespace pfl
