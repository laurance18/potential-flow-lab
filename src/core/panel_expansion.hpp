// panel_expansion.hpp: far-field (Laurent) series of groups of panels, to evaluate panel bodies and
// vortex sheets faster near them.
//
// The velocity of a set of panels lying within a radius r of a centre z_q is analytic outside that
// circle, so there it is a Laurent series
//     W(z) = sum_m a_m / (z - z_q)^m = sum_m b_m (r / (z - z_q))^m,     b_m = a_m / r^m,
// whose terms fall off like (r / |z - z_q|)^m. Beyond twice the radius, m up to 36.8 / ln(|z - z_q| / r)
// terms (at most 50) give full double precision, and cost a few multiplications each, against a
// logarithm and an arctangent per panel for the direct sum.
//
// A whole body already uses one such series beyond twice its radius. Closer in, its panels are
// grouped, k_clusterPanels consecutive ones at a time, each group with a series of its own: at a
// point near the body only the groups close to it are summed panel by panel. The result is the same
// to rounding (the Verify checks "Panel sums, near and far" and "Vortex panel sums" compare it with
// the panels summed directly). Near a 100-panel body it is about twice as fast; a body of 16 panels
// or fewer is one group and gains nothing.
#pragma once

#include "core/flow_element.hpp"
#include "core/vec2.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace pfl {

inline constexpr double k_farFieldRatio = 0.5; // a series is used where radius / |z - z_q| is below this
inline constexpr int k_laurentTerms = 50;      // 0.5^50 ~ 1e-15: full precision at the switch-over distance
inline constexpr std::size_t k_clusterPanels = 16; // fastest of 4 to 16 on the example scenes

// Panels [first, last) of a body or sheet, the circle that holds them, and their series.
struct PanelCluster {
    std::size_t first = 0;
    std::size_t last = 0;
    Vec2 centre;
    double radius = 0.0;
    std::vector<Complex> laurent; // b_1 .. b_M
};

// W from a series about `centre` at offset dz = z - centre (|dz| > radius / k_farFieldRatio), Horner
// in w = radius / dz with just enough terms for full precision.
inline Complex laurentVelocity(const std::vector<Complex>& laurent, double radius, Complex dz, double distance) {
    const Complex w = radius / dz;
    const double ratio = radius / distance;
    const int needed = static_cast<int>(std::ceil(-36.8 / std::log(ratio))) + 1;
    const int terms = needed < k_laurentTerms ? needed : k_laurentTerms;
    Complex sum{};
    for (int m = terms; m >= 1; --m) sum = (sum + laurent[static_cast<std::size_t>(m - 1)]) * w;
    return sum;
}

// True if p is far enough from the circle (centre, radius) for its series.
inline bool farFrom(Vec2 centre, double radius, Vec2 p) {
    return radius * radius < k_farFieldRatio * k_farFieldRatio * lengthSquared(p - centre);
}

// Whether a group within `radius` of its centre is summed by its series at squared distance d2 from
// that centre. A group of 16 panels needs at most 50 terms there, which cost less than its 16
// logarithms and arctangents.
inline bool useSeries(double radius, double d2) { return radius * radius < k_farFieldRatio * k_farFieldRatio * d2; }

// The clusters of `panels` panels, k_clusterPanels at a time, each around the panel ends
// `end(j)` = (start, finish) of its panels; the series are left empty.
template <typename EndsFn>
std::vector<PanelCluster> clusterPanels(std::size_t panels, const EndsFn& ends) {
    std::vector<PanelCluster> clusters;
    for (std::size_t first = 0; first < panels; first += k_clusterPanels) {
        PanelCluster c;
        c.first = first;
        c.last = first + k_clusterPanels < panels ? first + k_clusterPanels : panels;
        Rect box{ends(first).first, ends(first).first};
        for (std::size_t j = c.first; j < c.last; ++j) {
            for (const Vec2 p : {ends(j).first, ends(j).second}) {
                box.min = {std::min(box.min.x, p.x), std::min(box.min.y, p.y)};
                box.max = {std::max(box.max.x, p.x), std::max(box.max.y, p.y)};
            }
        }
        c.centre = box.center();
        for (std::size_t j = c.first; j < c.last; ++j) {
            c.radius = std::max({c.radius, length(ends(j).first - c.centre), length(ends(j).second - c.centre)});
        }
        clusters.push_back(std::move(c));
    }
    return clusters;
}

} // namespace pfl
