// panel_geometry.hpp: the shapes source-panel bodies are cut from, beyond the ellipse and the
// rectangle: NACA 4-digit airfoils and outlines given as points (drawn on the plot or read from an
// airfoil coordinate file), and how N panels are laid along them.
//
// NACA 4-digit "m p tt" (Abbott & von Doenhoff): maximum camber m % of the chord at p tenths of the
// chord, thickness tt % of the chord. With x along the chord (0 at the leading edge, 1 at the
// trailing edge):
//     y_t = 5 t (0.2969 sqrt(x) - 0.1260 x - 0.3516 x^2 + 0.2843 x^3 - 0.1036 x^4)
// (the last coefficient closes the trailing edge to a point, which panels need), the camber line
//     y_c = m / p^2 (2 p x - x^2) for x < p,  m / (1 - p)^2 ((1 - 2 p) + 2 p x - x^2) beyond,
// and the surfaces are offset from it along its normal: x -+ y_t sin(theta), y_c +- y_t cos(theta),
// tan(theta) = dy_c/dx. The panel vertices are cosine-spaced in x, N/2 on each surface, so the panels
// are short at the leading edge (where the curvature is) and at the trailing edge.
//
// An outline is a closed list of points, in either direction (it is turned clockwise), with no
// crossings. The points give the shape; the N panels are laid along it independently of them:
//   - smooth outline: a cubic spline through the points against the chord length between them (C2,
//     periodic around the outline), except at sharp corners (turning by more than 100 degrees, such
//     as a trailing edge), which stay sharp: there the spline is cut and ends with no curvature;
//   - straight outline: straight lines between the points, every kink a corner.
// The outline is cut at its corners, a corner always being a panel vertex, and each stretch between
// two corners gets panels in proportion to its length weighted by 1 + R kappa / 2 (kappa the
// curvature, R the perimeter / 2 pi, so a circle has equal panels), cosine-spaced towards those of
// its end corners that are sharp (turning by more than 30 degrees). Without corners the panels are
// equally spaced in that weighted length.
//
// Constant-strength panels converge at first order on such a spacing in general (the max Cp error
// halves when N doubles); the presets' spacings are special, with second order on an ellipse and
// exact Cp at the control points of a circle.
#pragma once

#include "core/vec2.hpp"

#include <optional>
#include <span>
#include <string>
#include <vector>

namespace pfl {

inline constexpr std::size_t k_maxOutlinePoints = 5000;

// Half the thickness y_t and the camber line y_c with its slope, at x along a unit chord; camber and
// position as fractions (m = 0.02, p = 0.4 for a NACA 24xx), t the thickness as a fraction.
double nacaHalfThickness(double t, double x);
double nacaCamber(double m, double p, double x);
double nacaCamberSlope(double m, double p, double x);

// The N (even) panel vertices of a NACA airfoil of unit chord from (0, 0) to (1, 0), clockwise: from
// the trailing edge along the lower surface to the leading edge, and back along the upper surface.
std::vector<Vec2> nacaVertices(double m, double p, double t, int panels);

// The outline without repeated points (nor its first point repeated at the end), clockwise.
std::vector<Vec2> cleanOutline(std::span<const Vec2> points);

// True if no two edges of the closed polygon cross or touch, apart from neighbours at their shared
// vertex.
bool isSimplePolygon(std::span<const Vec2> polygon);

// The number of corners of a clean outline (see above).
std::size_t outlineCorners(std::span<const Vec2> outline, bool smooth);

// The N panel vertices along a clean outline, clockwise.
std::vector<Vec2> distributeAlongOutline(std::span<const Vec2> outline, bool smooth, int panels);

// Why an outline cannot be cut into N panels (too few points, crossing edges, more corners than
// panels, panels that cross each other), or nullopt.
std::optional<std::string> outlineError(std::span<const Vec2> points, bool smooth, int panels);

} // namespace pfl
