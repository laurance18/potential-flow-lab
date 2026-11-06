#include "core/surface_pressure.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace pfl {
namespace {

constexpr double pi = std::numbers::pi;

constexpr int k_newtonIterations = 100;
constexpr double k_zeroSpeed = 1e-9;        // |W| / speedScale accepted as a stagnation point
constexpr double k_mergeDistance = 1e-4;    // stagnation points closer than this (x region size) are one
constexpr int k_directionSamples = 720;     // samples on the circle around a stagnation point
constexpr double k_directionOffset = 0.382; // of the first sample, in sample spacings
constexpr double k_circleRadius = 1e-4;     // its radius, x body region size (smaller near other points)
constexpr double k_snapDistance = 0.75;     // a traced line ends at a stagnation point this close (x step)
constexpr double k_cylinderTolerance = 1e-3; // surface points within this fraction of a of r = a
constexpr double k_closureTolerance = 0.1;  // a body closes if its surfaces pass the rear stagnation point
                                             // this close (x distance between the stagnation points)

Vec2 toVec(Complex z) { return {z.real(), z.imag()}; }

bool finite(Complex z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }

// A streamline traced from a stagnation point, starting with the stagnation point itself.
struct Branch {
    enum class End { Stagnation, Exit, Other };
    std::vector<Vec2> points;
    bool outgoing = true;
    End end = End::Other;
    std::size_t target = 0; // the stagnation point it reached, when end == Stagnation
    // Per stagnation point: the closest approach to it (after leaving the start) and that point's index.
    std::vector<std::pair<double, std::size_t>> closest;
};

// Shoelace formula; the polygon's last point repeats its first.
double signedArea(const std::vector<Vec2>& polygon) {
    double twice = 0.0;
    for (std::size_t k = 0; k + 1 < polygon.size(); ++k) twice += cross(polygon[k], polygon[k + 1]);
    return 0.5 * twice;
}

Vec2 areaCentroid(const std::vector<Vec2>& polygon) {
    const Vec2 origin = polygon.front(); // shifted for accuracy
    double twice = 0.0;
    Vec2 sum;
    for (std::size_t k = 0; k + 1 < polygon.size(); ++k) {
        const Vec2 a = polygon[k] - origin;
        const Vec2 b = polygon[k + 1] - origin;
        const double c = cross(a, b);
        twice += c;
        sum += (a + b) * c;
    }
    return twice != 0.0 ? origin + sum / (3.0 * twice) : origin;
}

// Even-odd rule.
bool insidePolygon(const std::vector<Vec2>& polygon, Vec2 p) {
    bool inside = false;
    for (std::size_t k = 0; k + 1 < polygon.size(); ++k) {
        const Vec2 a = polygon[k];
        const Vec2 b = polygon[k + 1];
        if ((a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y)) inside = !inside;
    }
    return inside;
}

// Where the singular points are, grown by the length scale of each element relative to the free
// stream (the size of the body it makes): every body lies inside. Nullopt without singular points.
std::optional<Rect> bodyRegion(const FlowField& field, double freeStreamSpeed) {
    std::optional<Rect> box;
    double scale = 0.0;
    const auto include = [&](const FlowElement& element) {
        const std::optional<Vec2> centre = singularPoint(element);
        if (!centre) return;
        box = box ? Rect{{std::min(box->min.x, centre->x), std::min(box->min.y, centre->y)},
                         {std::max(box->max.x, centre->x), std::max(box->max.y, centre->y)}}
                  : Rect{*centre, *centre};
        if (const auto* source = std::get_if<SourceSink>(&element)) {
            scale = std::max(scale, std::abs(source->strength) / freeStreamSpeed); // half-body width
        } else if (const auto* vortex = std::get_if<Vortex>(&element)) {
            scale = std::max(scale, std::abs(vortex->circulation) / (2.0 * pi * freeStreamSpeed));
        } else if (const auto* doublet = std::get_if<Doublet>(&element)) {
            scale = std::max(scale, std::sqrt(std::abs(doublet->strength) / freeStreamSpeed)); // cylinder radius
        }
    };
    for (const FlowElement& element : field.elements) include(element);
    for (const PeriodicRow& row : field.rows) include(row.element); // the other copies repeat the same bodies
    if (!box) return std::nullopt;
    const Vec2 size = box->size();
    const double margin = 2.0 * scale + 0.1 * std::max(size.x, size.y);
    if (margin <= 0.0) return std::nullopt;
    return Rect{box->min - Vec2{margin, margin}, box->max + Vec2{margin, margin}};
}

// Directions (angles) in which streamlines leave or reach the stagnation point s: where v is radial
// on a small circle around it, i.e. where cross(r, v) changes sign. Outgoing where dot(r, v) > 0.
std::vector<std::pair<double, bool>> branchDirections(const FlowField& field, Vec2 s, double radius) {
    const auto sample = [&](double angle) {
        const Vec2 r{std::cos(angle), std::sin(angle)};
        const Vec2 v = field.velocity(s + r * radius);
        return std::pair{cross(r, v), dot(r, v)};
    };
    std::vector<std::pair<double, bool>> directions;
    // Offset so that no sample falls on a symmetry axis of the flow: on a branch cross(r, v) is exactly
    // 0, which shows no sign change.
    const auto sampleAngle = [](int k) { return 2.0 * pi * (k + k_directionOffset) / k_directionSamples; };
    double t0 = sampleAngle(0);
    double g0 = sample(t0).first;
    for (int k = 1; k <= k_directionSamples; ++k) {
        double t1 = sampleAngle(k);
        const double g1 = sample(t1).first;
        if ((g0 > 0.0) != (g1 > 0.0)) {
            double a = t0;
            double b = t1;
            const bool aPositive = g0 > 0.0;
            for (int i = 0; i < 40; ++i) {
                const double m = 0.5 * (a + b);
                ((sample(m).first > 0.0) == aPositive ? a : b) = m;
            }
            const double angle = 0.5 * (a + b);
            directions.emplace_back(angle, sample(angle).second > 0.0);
        }
        t0 = t1;
        g0 = g1;
    }
    return directions;
}

struct TraceSettings {
    Rect bounds;
    double step = 0.0;
    double singularGuard = 0.0;
    int maxPoints = 0;
};

// Radius of the circle sampled around stagnation[i], small enough to hold no other stagnation or
// singular point.
double circleRadius(const FlowField& field, const std::vector<Vec2>& stagnation, std::size_t i, double extent) {
    double radius = k_circleRadius * extent;
    for (std::size_t k = 0; k < stagnation.size(); ++k) {
        if (k != i) radius = std::min(radius, 0.25 * length(stagnation[k] - stagnation[i]));
    }
    for (const FlowElement& element : field.elements) {
        if (const std::optional<Vec2> c = singularPoint(element)) radius = std::min(radius, 0.25 * length(*c - stagnation[i]));
    }
    for (const PeriodicRow& row : field.rows) {
        radius = std::min(radius, 0.25 * length(nearestCopy(row, stagnation[i]) - stagnation[i]));
    }
    return radius;
}

Branch traceBranch(const FlowField& field, const std::vector<Vec2>& stagnation, std::size_t start, double angle,
                   bool outgoing, double radius, const TraceSettings& settings) {
    Branch branch;
    branch.outgoing = outgoing;
    const Vec2 s = stagnation[start];
    Vec2 p = s + Vec2{std::cos(angle), std::sin(angle)} * radius;
    branch.points = {s, p};
    branch.closest.assign(stagnation.size(), {std::numeric_limits<double>::infinity(), 0});

    const double sign = outgoing ? 1.0 : -1.0;
    bool valid = true;
    const auto direction = [&](Vec2 q) {
        const Vec2 v = field.velocity(q);
        const double l = length(v);
        if (!(l > 0.0) || !std::isfinite(l)) {
            valid = false;
            return Vec2{};
        }
        return v * (sign / l);
    };

    const double snap = k_snapDistance * settings.step;
    double travelled = radius;
    for (int i = 0; i < settings.maxPoints; ++i) {
        double nearest = 1e300;
        for (std::size_t k = 0; k < stagnation.size(); ++k) {
            const double d = length(p - stagnation[k]);
            if (d < snap && (k != start || travelled > 10.0 * settings.step)) {
                branch.points.push_back(stagnation[k]);
                branch.end = Branch::End::Stagnation;
                branch.target = k;
                return branch;
            }
            nearest = std::min(nearest, d);
            if (k != start && d < branch.closest[k].first) branch.closest[k] = {d, branch.points.size() - 1};
        }
        // v/|v| turns quickly near a stagnation point, so the step shrinks there.
        const double h = std::clamp(0.1 * nearest, 0.01 * settings.step, settings.step);
        const Vec2 k1 = direction(p);
        const Vec2 k2 = direction(p + k1 * (0.5 * h));
        const Vec2 k3 = direction(p + k2 * (0.5 * h));
        const Vec2 k4 = direction(p + k3 * h);
        if (!valid) return branch;
        p += (k1 + 2.0 * k2 + 2.0 * k3 + k4) * (h / 6.0);
        travelled += h;
        branch.points.push_back(p);
        if (!settings.bounds.contains(p)) {
            branch.end = Branch::End::Exit;
            return branch;
        }
        if (!field.inDomain(p) || field.nearSingularity(p, settings.singularGuard)) return branch;
    }
    branch.end = Branch::End::Exit; // ran out of points: treated like leaving the region
    return branch;
}

// Fills in Cp and the coordinates along the surface from its points.
void evaluateSurface(BodySurface& body, const FlowModel& model, std::size_t stagnationIndex) {
    const FlowField& field = model.field;
    const Vec2 e = model.freeStream / model.freeStreamSpeed;
    const std::size_t n = body.points.size();
    body.cp.resize(n);
    body.arcLength.resize(n);
    body.streamwise.resize(n);
    double s = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        if (k > 0) s += length(body.points[k] - body.points[k - 1]);
        body.arcLength[k] = s;
        body.cp[k] = pressureCoefficient(length(field.velocity(body.points[k])), model.freeStreamSpeed);
    }
    const double s0 = body.arcLength[stagnationIndex];
    for (double& value : body.arcLength) value -= s0;

    if (body.closed) {
        body.centroid = areaCentroid(body.points);
        double lo = 1e300;
        double hi = -1e300;
        for (const Vec2 p : body.points) {
            lo = std::min(lo, dot(p, e));
            hi = std::max(hi, dot(p, e));
        }
        body.angle.resize(n);
        for (std::size_t k = 0; k < n; ++k) {
            const Vec2 r = body.points[k] - body.centroid;
            body.streamwise[k] = hi > lo ? (dot(body.points[k], e) - lo) / (hi - lo) : 0.0;
            const double degrees = std::atan2(cross(e, r), dot(e, r)) * 180.0 / pi;
            body.angle[k] = degrees < 0.0 ? degrees + 360.0 : degrees;
        }
    } else {
        body.centroid = body.stagnationPoints.front();
        for (std::size_t k = 0; k < n; ++k) body.streamwise[k] = dot(body.points[k] - body.centroid, e);
    }
}

BodySurface closedBody(std::vector<Vec2> polygon, std::vector<Vec2> stagnationPoints, double gap,
                       const FlowModel& model) {
    BodySurface body;
    body.closed = true;
    body.closureGap = gap;
    if (signedArea(polygon) > 0.0) std::ranges::reverse(polygon); // clockwise; still starts at the front
    body.points = std::move(polygon);
    body.stagnationPoints = std::move(stagnationPoints);
    evaluateSurface(body, model, 0);
    return body;
}

// The part of `branch` up to stagnation point j (which is `sj`) and how far it misses it: 0 if the
// branch ends there, else its closest approach to sj, accepted up to `tolerance`.
std::optional<std::pair<std::vector<Vec2>, double>> pathTo(const Branch& branch, std::size_t j, Vec2 sj, double tolerance) {
    if (branch.end == Branch::End::Stagnation && branch.target == j) return std::pair{branch.points, 0.0};
    const auto [distance, index] = branch.closest[j];
    if (!(distance <= tolerance)) return std::nullopt;
    std::vector<Vec2> points(branch.points.begin(), branch.points.begin() + static_cast<std::ptrdiff_t>(index) + 1);
    points.push_back(sj);
    return std::pair{std::move(points), distance};
}

// W is analytic, so its derivative along x is dW/dz. Damped far from a root; stops at a root or on
// leaving `searchArea`.
Complex newtonOnVelocity(const FlowField& field, Complex z, const Rect& searchArea, double extent) {
    const double h = 1e-6 * extent;
    const auto w = [&field](Complex q) { return field.complexVelocity(toVec(q)); };
    for (int iteration = 0; iteration < k_newtonIterations; ++iteration) {
        const Complex value = w(z);
        const Complex derivative = (w(z + h) - w(z - h)) / (2.0 * h);
        if (!finite(value) || !finite(derivative) || std::abs(derivative) == 0.0) break;
        Complex step = value / derivative;
        const double stepLength = std::abs(step);
        if (stepLength > 0.1 * extent) step *= 0.1 * extent / stepLength;
        z -= step;
        if (!searchArea.contains(toVec(z)) || stepLength < 1e-14 * extent) break;
    }
    return z;
}

bool isStagnationPoint(const FlowField& field, Vec2 p, double speedScale, double extent) {
    return field.inDomain(p) && !field.nearSingularity(p, 1e-6 * extent) &&
           std::abs(field.complexVelocity(p)) <= k_zeroSpeed * speedScale;
}

BodySurface openBody(const Branch& a, const Branch& b, const FlowModel& model) {
    // Body on the right, like the clockwise closed surfaces. A source's half-body extends downstream,
    // so that runs in along the lower branch and out along the upper one; a sink's the other way.
    const Vec2 s = a.points.front();
    const bool aUpper = cross(model.freeStream, a.points.back() - s) > 0.0;
    const bool downstream = dot(model.freeStream, a.points.back() - s) > 0.0;
    const Branch& in = aUpper == downstream ? b : a;
    const Branch& out = aUpper == downstream ? a : b;
    BodySurface body;
    body.points.assign(in.points.rbegin(), in.points.rend());
    body.points.insert(body.points.end(), out.points.begin() + 1, out.points.end());
    body.stagnationPoints = {s};
    evaluateSurface(body, model, in.points.size() - 1);
    return body;
}

} // namespace

std::vector<Vec2> findStagnationPoints(const FlowField& field, const Rect& region, double speedScale, int seedGrid) {
    std::vector<Vec2> found;
    const Vec2 size = region.size();
    const double extent = std::max(size.x, size.y);
    if (!(extent > 0.0) || !(speedScale > 0.0)) return found;
    const Rect searchArea = region.expanded(0.1);

    for (int j = 0; j < seedGrid; ++j) {
        for (int i = 0; i < seedGrid; ++i) {
            const Vec2 seed = region.min + Vec2{size.x * (i + 0.5) / seedGrid, size.y * (j + 0.5) / seedGrid};
            const Vec2 p = toVec(newtonOnVelocity(field, toComplex(seed), searchArea, extent));
            if (!region.contains(p) || !isStagnationPoint(field, p, speedScale, extent)) continue;
            const bool known = std::ranges::any_of(found, [&](Vec2 q) {
                return lengthSquared(p - q) < k_mergeDistance * k_mergeDistance * extent * extent;
            });
            if (!known) found.push_back(p);
        }
    }
    return found;
}

std::optional<Vec2> followStagnationPoint(const FlowField& field, Vec2 from, double reach, double speedScale) {
    if (!(reach > 0.0) || !(speedScale > 0.0)) return std::nullopt;
    const Vec2 corner{reach, reach};
    const Vec2 p = toVec(newtonOnVelocity(field, toComplex(from), {from - corner, from + corner}, reach));
    if (!(length(p - from) <= reach) || !isStagnationPoint(field, p, speedScale, reach)) return std::nullopt;
    return p;
}

StreamlineSet dividingStreamlines(const FlowField& field, Vec2 point, std::vector<Vec2> stagnation, const Rect& bounds,
                                  const SurfaceSettings& settings) {
    StreamlineSet set;
    const double extent = std::max(bounds.size().x, bounds.size().y);
    if (!(extent > 0.0)) return set;

    std::erase_if(stagnation, [&](Vec2 q) { return length(q - point) < k_mergeDistance * extent; });
    stagnation.insert(stagnation.begin(), point);
    TraceSettings trace;
    trace.bounds = {{std::min(bounds.min.x, point.x), std::min(bounds.min.y, point.y)},
                    {std::max(bounds.max.x, point.x), std::max(bounds.max.y, point.y)}};
    trace.step = settings.stepFraction * extent;
    trace.singularGuard = 2.0 * trace.step;
    trace.maxPoints = settings.maxPointsPerBranch;

    const double radius = circleRadius(field, stagnation, 0, extent);
    for (const auto& [angle, outgoing] : branchDirections(field, point, radius)) {
        // At a sector vertex the velocity formula also has radial directions outside the wedge.
        if (!field.inDomain(point + Vec2{std::cos(angle), std::sin(angle)} * radius)) continue;
        std::vector<Vec2> line = traceBranch(field, stagnation, 0, angle, outgoing, radius, trace).points;
        if (!outgoing) std::ranges::reverse(line);

        double total = 0.0;
        for (std::size_t k = 0; k + 1 < line.size(); ++k) total += length(line[k + 1] - line[k]);
        double along = 0.0;
        for (std::size_t k = 0; k + 1 < line.size(); ++k) {
            const Vec2 d = line[k + 1] - line[k];
            const double l = length(d);
            along += l;
            if (along >= 0.5 * total && l > 0.0) {
                set.arrows.push_back({line[k], d / l});
                break;
            }
        }
        set.lines.push_back(std::move(line));
    }
    return set;
}

std::vector<BodySurface> findBodySurfaces(const FlowModel& model, const Rect& view, const SurfaceSettings& settings) {
    std::vector<BodySurface> bodies;
    if (!(model.freeStreamSpeed > 0.0)) return bodies;
    const FlowField& field = model.field;
    const std::optional<Rect> region = bodyRegion(field, model.freeStreamSpeed);
    if (!region) return bodies;

    const std::vector<Vec2> stagnation =
        findStagnationPoints(field, *region, model.freeStreamSpeed, settings.seedGrid);
    const double extent = std::max(region->size().x, region->size().y);
    const Rect shown = view.expanded(0.05);
    TraceSettings trace;
    trace.bounds = {{std::min(region->min.x, shown.min.x), std::min(region->min.y, shown.min.y)},
                    {std::max(region->max.x, shown.max.x), std::max(region->max.y, shown.max.y)}};
    trace.step = settings.stepFraction * extent;
    trace.singularGuard = 2.0 * trace.step;
    trace.maxPoints = settings.maxPointsPerBranch;

    std::vector<std::vector<Branch>> branches(stagnation.size());
    for (std::size_t i = 0; i < stagnation.size(); ++i) {
        const double radius = circleRadius(field, stagnation, i, extent);
        for (const auto& [angle, outgoing] : branchDirections(field, stagnation[i], radius)) {
            branches[i].push_back(traceBranch(field, stagnation, i, angle, outgoing, radius, trace));
        }
    }

    for (std::size_t i = 0; i < stagnation.size(); ++i) {
        const std::vector<Branch>& from = branches[i];
        std::vector<const Branch*> exitsOut;
        std::vector<const Branch*> exitsIn;
        const std::size_t bodiesBefore = bodies.size();
        for (std::size_t a = 0; a < from.size(); ++a) {
            const Branch& first = from[a];
            if (first.end == Branch::End::Exit) (first.outgoing ? exitsOut : exitsIn).push_back(&first);
            if (first.outgoing && first.end == Branch::End::Stagnation && first.target == i) {
                // A streamline that leaves the stagnation point and returns to it.
                bodies.push_back(closedBody(first.points, {stagnation[i]}, 0.0, model));
            }
        }
        // Two streamlines from stagnation point i to the same point j: upper and lower surface. Unless the
        // flow is symmetric they miss j (other bodies nearby), which is accepted up to k_closureTolerance
        // and reported as the body's closure gap.
        for (std::size_t a = 0; a < from.size(); ++a) {
            for (std::size_t b = a + 1; b < from.size(); ++b) {
                if (!from[a].outgoing || !from[b].outgoing) continue;
                std::optional<std::pair<std::vector<Vec2>, double>> best[2];
                std::size_t bestTarget = 0;
                for (std::size_t j = 0; j < stagnation.size(); ++j) {
                    if (j == i) continue;
                    const double tolerance = k_closureTolerance * length(stagnation[j] - stagnation[i]);
                    auto upper = pathTo(from[a], j, stagnation[j], tolerance);
                    auto lower = pathTo(from[b], j, stagnation[j], tolerance);
                    if (!upper || !lower) continue;
                    const double gap = std::max(upper->second, lower->second);
                    if (best[0] && gap >= std::max(best[0]->second, best[1]->second)) continue;
                    best[0] = std::move(upper);
                    best[1] = std::move(lower);
                    bestTarget = j;
                }
                if (!best[0]) continue;
                std::vector<Vec2> polygon = std::move(best[0]->first);
                polygon.insert(polygon.end(), best[1]->first.rbegin() + 1, best[1]->first.rend());
                bodies.push_back(closedBody(std::move(polygon), {stagnation[i], stagnation[bestTarget]},
                                            std::max(best[0]->second, best[1]->second), model));
            }
        }
        if (bodies.size() > bodiesBefore) continue; // the branches that left the region belong to a closed body
        if (exitsOut.size() == 2) bodies.push_back(openBody(*exitsOut[0], *exitsOut[1], model));
        if (exitsIn.size() == 2) bodies.push_back(openBody(*exitsIn[0], *exitsIn[1], model));
    }
    return bodies;
}

BodyForces bodyForces(const BodySurface& body, const FlowModel& model) {
    BodyForces forces;
    if (!body.closed || body.points.size() < 3 || !(model.freeStreamSpeed > 0.0)) return forces;
    const Vec2 e = model.freeStream / model.freeStreamSpeed;

    // Clockwise polygon: perp(d) is the outward normal times the segment length. The constant part of
    // the pressure integrates to zero over the closed polygon, so only Cp matters.
    Vec2 force;
    double lo = 1e300;
    double hi = -1e300;
    for (std::size_t k = 0; k + 1 < body.points.size(); ++k) {
        const Vec2 d = body.points[k + 1] - body.points[k];
        force -= perp(d) * (0.5 * (body.cp[k] + body.cp[k + 1]));
        lo = std::min(lo, dot(body.points[k], e));
        hi = std::max(hi, dot(body.points[k], e));
    }
    forces.drag = dot(force, e);
    forces.lift = dot(force, perp(e));
    forces.chord = hi - lo;

    forces.isolated = true;
    for (const FlowElement& element : model.field.elements) {
        if (std::holds_alternative<SectorFlow>(element)) forces.isolated = false;
        const std::optional<Vec2> centre = singularPoint(element);
        if (!centre) continue;
        const bool inside = insidePolygon(body.points, *centre);
        if (!inside) forces.isolated = false;
        if (const auto* vortex = std::get_if<Vortex>(&element); vortex != nullptr && inside) {
            forces.circulation += vortex->circulation;
        }
    }
    // In a channel the body feels every other copy, and holds at most the nearest one of a row.
    for (const PeriodicRow& row : model.field.rows) {
        forces.isolated = false;
        const auto* vortex = std::get_if<Vortex>(&row.element);
        if (vortex != nullptr && insidePolygon(body.points, nearestCopy(row, body.points.front()))) {
            forces.circulation += vortex->circulation;
        }
    }
    return forces;
}

double CylinderReference::cp(double theta) const {
    const double tangential = 2.0 * std::sin(theta - freeStreamAngle) + circulation / (2.0 * pi * radius * freeStreamSpeed);
    return 1.0 - tangential * tangential;
}

std::optional<CylinderReference> cylinderReference(const BodySurface& body, const FlowModel& model) {
    if (!body.closed || !(model.freeStreamSpeed > 0.0) || !model.field.rows.empty()) return std::nullopt;
    const Doublet* doublet = nullptr;
    const Vortex* vortex = nullptr;
    for (const FlowElement& element : model.field.elements) {
        if (std::holds_alternative<UniformFlow>(element)) continue;
        if (const auto* d = std::get_if<Doublet>(&element); d != nullptr && doublet == nullptr) {
            doublet = d;
        } else if (const auto* v = std::get_if<Vortex>(&element); v != nullptr && vortex == nullptr) {
            vortex = v;
        } else {
            return std::nullopt;
        }
    }
    if (doublet == nullptr || !(doublet->strength > 0.0)) return std::nullopt;

    CylinderReference reference;
    reference.centre = doublet->position;
    reference.freeStreamSpeed = model.freeStreamSpeed;
    reference.freeStreamAngle = std::atan2(model.freeStream.y, model.freeStream.x);
    reference.radius = std::sqrt(doublet->strength / model.freeStreamSpeed);
    if (std::cos(doublet->angle - reference.freeStreamAngle) < 1.0 - 1e-9) return std::nullopt;
    if (vortex != nullptr) {
        if (length(vortex->position - doublet->position) > 1e-9 * reference.radius) return std::nullopt;
        reference.circulation = vortex->circulation;
    }
    const bool onCircle = std::ranges::all_of(body.points, [&](Vec2 p) {
        return std::abs(length(p - reference.centre) - reference.radius) <= k_cylinderTolerance * reference.radius;
    });
    if (!onCircle) return std::nullopt;
    return reference;
}

} // namespace pfl
