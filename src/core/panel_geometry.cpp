#include "core/panel_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <numeric>

namespace pfl {
namespace {

constexpr double pi = std::numbers::pi;
constexpr double k_smoothCorner = 100.0 * pi / 180.0;  // a smooth outline keeps sharper turns as corners
constexpr double k_straightCorner = 0.5 * pi / 180.0;  // a straight outline: every kink but rounding
constexpr double k_clusterCorner = 30.0 * pi / 180.0;  // panels close up towards corners sharper than this
constexpr int k_splineSamples = 24;                    // per stretch between two points of a smooth outline
constexpr double k_curvatureWeight = 0.5;
constexpr double k_curvatureCap = 40.0; // of R kappa, so a near-cusp does not take every panel

double signedArea(std::span<const Vec2> v) {
    double twice = 0.0;
    for (std::size_t k = 0; k < v.size(); ++k) twice += cross(v[k], v[(k + 1) % v.size()]);
    return 0.5 * twice;
}

// Turning angle at b on the way a -> b -> c, in [0, pi].
double turn(Vec2 a, Vec2 b, Vec2 c) { return std::abs(std::atan2(cross(b - a, c - b), dot(b - a, c - b))); }

double orientation(Vec2 a, Vec2 b, Vec2 c) { return cross(b - a, c - a); }

bool onSegment(Vec2 a, Vec2 b, Vec2 p) {
    return std::min(a.x, b.x) <= p.x && p.x <= std::max(a.x, b.x) && std::min(a.y, b.y) <= p.y && p.y <= std::max(a.y, b.y);
}

bool segmentsMeet(Vec2 a, Vec2 b, Vec2 c, Vec2 d) {
    const double o1 = orientation(a, b, c);
    const double o2 = orientation(a, b, d);
    const double o3 = orientation(c, d, a);
    const double o4 = orientation(c, d, b);
    if (((o1 > 0.0 && o2 < 0.0) || (o1 < 0.0 && o2 > 0.0)) && ((o3 > 0.0 && o4 < 0.0) || (o3 < 0.0 && o4 > 0.0))) return true;
    return (o1 == 0.0 && onSegment(a, b, c)) || (o2 == 0.0 && onSegment(a, b, d)) || (o3 == 0.0 && onSegment(c, d, a)) ||
           (o4 == 0.0 && onSegment(c, d, b));
}

// Solves the tridiagonal system sub[i] x[i-1] + diag[i] x[i] + super[i] x[i+1] = rhs[i] (Thomas).
std::vector<Vec2> solveTridiagonal(std::vector<double> sub, std::vector<double> diag, std::vector<double> super,
                                   std::vector<Vec2> rhs) {
    const std::size_t n = diag.size();
    for (std::size_t i = 1; i < n; ++i) {
        const double f = sub[i] / diag[i - 1];
        diag[i] -= f * super[i - 1];
        rhs[i] -= rhs[i - 1] * f;
    }
    std::vector<Vec2> x(n);
    x[n - 1] = rhs[n - 1] / diag[n - 1];
    for (std::size_t i = n - 1; i-- > 0;) x[i] = (rhs[i] - x[i + 1] * super[i]) / diag[i];
    return x;
}

// A cubic spline through points q (in x and y, against the chord length s between them, so the
// points need not be evenly spaced), C2 throughout: natural (no curvature) at the two ends of an open
// run, periodic around a closed one. Returns the second derivatives at the points.
std::vector<Vec2> splineCurvatures(std::span<const Vec2> q, std::span<const double> h, bool periodic) {
    const std::size_t n = q.size();
    if (!periodic) {
        // Unknowns M_1 .. M_(n-2); M_0 = M_(n-1) = 0.
        if (n < 3) return std::vector<Vec2>(n);
        const std::size_t m = n - 2;
        std::vector<double> sub(m), diag(m), super(m);
        std::vector<Vec2> rhs(m);
        for (std::size_t i = 1; i + 1 < n; ++i) {
            sub[i - 1] = h[i - 1];
            diag[i - 1] = 2.0 * (h[i - 1] + h[i]);
            super[i - 1] = h[i];
            rhs[i - 1] = ((q[i + 1] - q[i]) / h[i] - (q[i] - q[i - 1]) / h[i - 1]) * 6.0;
        }
        std::vector<Vec2> inner = solveTridiagonal(sub, diag, super, rhs);
        std::vector<Vec2> curvature(n);
        std::ranges::copy(inner, curvature.begin() + 1);
        return curvature;
    }
    // Cyclic: h[i] from q[i] to q[i + 1 mod n]. Sherman-Morrison on the two corner coefficients.
    std::vector<double> sub(n), diag(n), super(n);
    std::vector<Vec2> rhs(n);
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t before = (i + n - 1) % n;
        sub[i] = h[before];
        diag[i] = 2.0 * (h[before] + h[i]);
        super[i] = h[i];
        rhs[i] = ((q[(i + 1) % n] - q[i]) / h[i] - (q[i] - q[before]) / h[before]) * 6.0;
    }
    const double alpha = super[n - 1]; // row n-1, column 0
    const double beta = sub[0];        // row 0, column n-1
    const double gamma = -diag[0];
    std::vector<double> modified = diag;
    modified[0] -= gamma;
    modified[n - 1] -= alpha * beta / gamma;
    const std::vector<Vec2> x = solveTridiagonal(sub, modified, super, rhs);
    std::vector<Vec2> u(n);
    u[0] = {gamma, gamma};
    u[n - 1] = {alpha, alpha};
    const std::vector<Vec2> z = solveTridiagonal(sub, modified, super, u);
    std::vector<Vec2> curvature(n);
    const auto component = [&](double Vec2::*c) {
        const double fact = (x[0].*c + beta * (x[n - 1].*c) / gamma) / (1.0 + z[0].*c + beta * (z[n - 1].*c) / gamma);
        for (std::size_t i = 0; i < n; ++i) curvature[i].*c = x[i].*c - fact * (z[i].*c);
    };
    component(&Vec2::x);
    component(&Vec2::y);
    return curvature;
}

std::vector<bool> cornersOf(std::span<const Vec2> p, bool smooth) {
    const std::size_t n = p.size();
    const double threshold = smooth ? k_smoothCorner : k_straightCorner;
    std::vector<bool> corner(n);
    for (std::size_t i = 0; i < n; ++i) corner[i] = turn(p[(i + n - 1) % n], p[i], p[(i + 1) % n]) > threshold;
    return corner;
}

// The outline as a fine polyline (closed, not repeated at the end), with its corners marked.
struct DenseOutline {
    std::vector<Vec2> points;
    std::vector<bool> corner;
};

// Samples of the spline through `q` (h: chord lengths, curvature: from splineCurvatures) on every
// stretch from q[i] to q[i + 1], q[i] itself included, the end of the last stretch not.
void sampleSpline(std::span<const Vec2> q, std::span<const double> h, std::span<const Vec2> curvature, bool periodic,
                  std::vector<Vec2>& out) {
    const std::size_t n = q.size();
    const std::size_t stretches = periodic ? n : n - 1;
    for (std::size_t i = 0; i < stretches; ++i) {
        const std::size_t j = (i + 1) % n;
        for (int k = 0; k < k_splineSamples; ++k) {
            const double b = static_cast<double>(k) / k_splineSamples;
            const double a = 1.0 - b;
            out.push_back(q[i] * a + q[j] * b +
                          (curvature[i] * (a * a * a - a) + curvature[j] * (b * b * b - b)) * (h[i] * h[i] / 6.0));
        }
    }
}

DenseOutline densify(std::span<const Vec2> p, bool smooth) {
    const std::size_t n = p.size();
    const std::vector<bool> corner = cornersOf(p, smooth);
    DenseOutline dense;
    if (!smooth) {
        dense.points.assign(p.begin(), p.end());
        dense.corner = corner;
        return dense;
    }
    std::vector<std::size_t> corners;
    for (std::size_t i = 0; i < n; ++i) {
        if (corner[i]) corners.push_back(i);
    }
    const auto chords = [](std::span<const Vec2> q, bool periodic) {
        std::vector<double> h(q.size());
        for (std::size_t i = 0; i + (periodic ? 0 : 1) < q.size(); ++i) h[i] = length(q[(i + 1) % q.size()] - q[i]);
        return h;
    };
    if (corners.empty()) {
        const std::vector<double> h = chords(p, true);
        sampleSpline(p, h, splineCurvatures(p, h, true), true, dense.points);
        dense.corner.assign(dense.points.size(), false);
        return dense;
    }
    // One open spline per run from a corner to the next, natural at both (sharp) ends.
    for (std::size_t r = 0; r < corners.size(); ++r) {
        std::vector<Vec2> run;
        for (std::size_t i = corners[r];; i = (i + 1) % n) {
            run.push_back(p[i]);
            if (run.size() > 1 && i == corners[(r + 1) % corners.size()]) break;
        }
        const std::vector<double> h = chords(run, false);
        const std::size_t before = dense.points.size();
        sampleSpline(run, h, splineCurvatures(run, h, false), false, dense.points);
        dense.corner.resize(dense.points.size(), false);
        dense.corner[before] = true;
    }
    return dense;
}

} // namespace

double nacaHalfThickness(double t, double x) {
    return 5.0 * t * (0.2969 * std::sqrt(x) - 0.1260 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1036 * x * x * x * x);
}

double nacaCamber(double m, double p, double x) {
    if (m == 0.0 || p <= 0.0) return 0.0;
    return x < p ? m / (p * p) * (2.0 * p * x - x * x) : m / ((1.0 - p) * (1.0 - p)) * ((1.0 - 2.0 * p) + 2.0 * p * x - x * x);
}

double nacaCamberSlope(double m, double p, double x) {
    if (m == 0.0 || p <= 0.0) return 0.0;
    return x < p ? 2.0 * m / (p * p) * (p - x) : 2.0 * m / ((1.0 - p) * (1.0 - p)) * (p - x);
}

std::vector<Vec2> nacaVertices(double m, double p, double t, int panels) {
    const int half = panels / 2;
    const auto surface = [&](int i, bool upper) {
        const double x = 0.5 * (1.0 - std::cos(pi * i / half));
        const double yt = nacaHalfThickness(t, x);
        const double theta = std::atan(nacaCamberSlope(m, p, x));
        const double sign = upper ? 1.0 : -1.0;
        return Vec2{x - sign * yt * std::sin(theta), nacaCamber(m, p, x) + sign * yt * std::cos(theta)};
    };
    std::vector<Vec2> vertices;
    vertices.reserve(static_cast<std::size_t>(2 * half));
    vertices.push_back({1.0, 0.0}); // the closed trailing edge, shared by both surfaces
    for (int i = half - 1; i >= 1; --i) vertices.push_back(surface(i, false));
    vertices.push_back({0.0, 0.0}); // the leading edge
    for (int i = 1; i <= half - 1; ++i) vertices.push_back(surface(i, true));
    return vertices;
}

std::vector<Vec2> cleanOutline(std::span<const Vec2> points) {
    std::vector<Vec2> clean;
    if (points.empty()) return clean;
    Vec2 lo = points.front();
    Vec2 hi = points.front();
    for (const Vec2 p : points) {
        lo = {std::min(lo.x, p.x), std::min(lo.y, p.y)};
        hi = {std::max(hi.x, p.x), std::max(hi.y, p.y)};
    }
    const double tolerance = 1e-12 * std::max(length(hi - lo), 1e-300);
    for (const Vec2 p : points) {
        if (clean.empty() || length(p - clean.back()) > tolerance) clean.push_back(p);
    }
    while (clean.size() > 1 && length(clean.back() - clean.front()) <= tolerance) clean.pop_back();
    if (signedArea(clean) > 0.0) std::ranges::reverse(clean);
    return clean;
}

bool isSimplePolygon(std::span<const Vec2> v) {
    const std::size_t n = v.size();
    if (n < 3) return false;
    // A spike (an edge doubling straight back) is not simple either.
    for (std::size_t i = 0; i < n; ++i) {
        const Vec2 a = v[(i + n - 1) % n];
        const Vec2 b = v[i];
        const Vec2 c = v[(i + 1) % n];
        if (cross(b - a, c - b) == 0.0 && dot(b - a, c - b) < 0.0) return false;
    }
    // Sweep in x: only edges whose x ranges overlap are compared.
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    const auto minX = [&](std::size_t e) { return std::min(v[e].x, v[(e + 1) % n].x); };
    const auto maxX = [&](std::size_t e) { return std::max(v[e].x, v[(e + 1) % n].x); };
    std::ranges::sort(order, {}, minX);
    for (std::size_t a = 0; a < n; ++a) {
        const std::size_t e = order[a];
        for (std::size_t b = a + 1; b < n && minX(order[b]) <= maxX(e); ++b) {
            const std::size_t f = order[b];
            const bool neighbours = (e + 1) % n == f || (f + 1) % n == e;
            if (neighbours) continue;
            if (segmentsMeet(v[e], v[(e + 1) % n], v[f], v[(f + 1) % n])) return false;
        }
    }
    return true;
}

std::size_t outlineCorners(std::span<const Vec2> outline, bool smooth) {
    const std::vector<bool> corner = cornersOf(outline, smooth);
    return static_cast<std::size_t>(std::ranges::count(corner, true));
}

std::vector<Vec2> distributeAlongOutline(std::span<const Vec2> outline, bool smooth, int panels) {
    const DenseOutline dense = densify(outline, smooth);
    const std::vector<Vec2>& d = dense.points;
    const std::size_t m = d.size();

    // Edge lengths and curvature-weighted lengths (the curvature is 0 at the corners: their
    // clustering is the cosine spacing).
    std::vector<double> edge(m);
    double perimeter = 0.0;
    for (std::size_t e = 0; e < m; ++e) {
        edge[e] = length(d[(e + 1) % m] - d[e]);
        perimeter += edge[e];
    }
    const double radius = perimeter / (2.0 * pi);
    std::vector<double> curvature(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
        if (dense.corner[i]) continue;
        const double ds = 0.5 * (edge[(i + m - 1) % m] + edge[i]);
        if (ds > 0.0) curvature[i] = std::min(radius * turn(d[(i + m - 1) % m], d[i], d[(i + 1) % m]) / ds, k_curvatureCap);
    }
    std::vector<double> weight(m);
    for (std::size_t e = 0; e < m; ++e) weight[e] = edge[e] * (1.0 + k_curvatureWeight * 0.5 * (curvature[e] + curvature[(e + 1) % m]));

    // Weighted length from dense point 0 to each dense point, twice around so a run may wrap.
    std::vector<double> cumulative(2 * m + 1, 0.0);
    for (std::size_t e = 0; e < 2 * m; ++e) cumulative[e + 1] = cumulative[e] + weight[e % m];
    // The point at weighted length `target` along the outline from dense point `start`.
    const auto pointAlong = [&](std::size_t start, double target) {
        const double goal = cumulative[start] + target;
        const auto it = std::upper_bound(cumulative.begin() + static_cast<std::ptrdiff_t>(start),
                                         cumulative.begin() + static_cast<std::ptrdiff_t>(start + m), goal);
        const std::size_t e = static_cast<std::size_t>(it - cumulative.begin()) - 1;
        const double w = cumulative[e + 1] - cumulative[e];
        const double f = w > 0.0 ? std::clamp((goal - cumulative[e]) / w, 0.0, 1.0) : 0.0;
        return d[e % m] + (d[(e + 1) % m] - d[e % m]) * f;
    };

    std::vector<std::size_t> corners;
    for (std::size_t i = 0; i < m; ++i) {
        if (dense.corner[i]) corners.push_back(i);
    }
    std::vector<Vec2> vertices;
    vertices.reserve(static_cast<std::size_t>(panels));
    const double total = std::accumulate(weight.begin(), weight.end(), 0.0);
    if (corners.empty()) {
        for (int k = 0; k < panels; ++k) vertices.push_back(pointAlong(0, total * k / panels));
        return vertices;
    }

    // Panels per run between two corners: in proportion to its weighted length, at least one each.
    const std::size_t runs = corners.size();
    std::vector<double> runWeight(runs, 0.0);
    for (std::size_t r = 0; r < runs; ++r) {
        const std::size_t end = corners[(r + 1) % runs];
        std::size_t e = corners[r];
        do {
            runWeight[r] += weight[e];
            e = (e + 1) % m;
        } while (e != end);
    }
    std::vector<int> count(runs);
    int assigned = 0;
    for (std::size_t r = 0; r < runs; ++r) {
        count[r] = std::max(1, static_cast<int>(std::floor(panels * runWeight[r] / total)));
        assigned += count[r];
    }
    const auto share = [&](std::size_t r) { return panels * runWeight[r] / total - count[r]; };
    while (assigned < panels) {
        std::size_t best = 0;
        for (std::size_t r = 1; r < runs; ++r) {
            if (share(r) > share(best)) best = r;
        }
        ++count[best];
        ++assigned;
    }
    while (assigned > panels) {
        std::size_t best = runs;
        for (std::size_t r = 0; r < runs; ++r) {
            if (count[r] > 1 && (best == runs || share(r) < share(best))) best = r;
        }
        if (best == runs) break; // more corners than panels: outlineError() rejects it
        --count[best];
        --assigned;
    }
    // Panels get shorter towards a sharp corner (cosine spacing), not towards a mere kink of a
    // straight outline, where the spacing carries on smoothly.
    const auto sharp = [&](std::size_t i) { return turn(d[(i + m - 1) % m], d[i], d[(i + 1) % m]) > k_clusterCorner; };
    for (std::size_t r = 0; r < runs; ++r) {
        const bool first = sharp(corners[r]);
        const bool last = sharp(corners[(r + 1) % runs]);
        for (int k = 0; k < count[r]; ++k) {
            const double t = static_cast<double>(k) / count[r];
            const double u = first && last ? 0.5 * (1.0 - std::cos(pi * t))
                             : first       ? 1.0 - std::cos(0.5 * pi * t)
                             : last        ? std::sin(0.5 * pi * t)
                                           : t;
            vertices.push_back(pointAlong(corners[r], u * runWeight[r]));
        }
    }
    return vertices;
}

std::optional<std::string> outlineError(std::span<const Vec2> points, bool smooth, int panels) {
    if (!std::ranges::all_of(points, [](Vec2 p) { return std::isfinite(p.x) && std::isfinite(p.y); })) {
        return "The outline's points must be finite numbers.";
    }
    const std::vector<Vec2> outline = cleanOutline(points);
    if (outline.size() < 3) return "The outline needs at least 3 different points (draw one, or import airfoil coordinates).";
    if (outline.size() > k_maxOutlinePoints) {
        return "The outline has more than " + std::to_string(k_maxOutlinePoints) + " points.";
    }
    if (std::abs(signedArea(outline)) == 0.0) return "The outline encloses no area.";
    if (!isSimplePolygon(outline)) return "The outline crosses itself.";
    const std::size_t corners = outlineCorners(outline, smooth);
    if (corners > static_cast<std::size_t>(panels)) {
        return "The outline has " + std::to_string(corners) + (smooth ? " sharp corners" : " corners") +
               ", each a panel vertex, so it needs at least that many panels N" +
               (smooth ? "." : " (or a smooth outline, which keeps only the sharp ones).");
    }
    const std::vector<Vec2> vertices = distributeAlongOutline(outline, smooth, panels);
    if (!isSimplePolygon(vertices)) {
        return "With N = " + std::to_string(panels) + " the panels cut across each other; use more panels.";
    }
    return std::nullopt;
}

} // namespace pfl
