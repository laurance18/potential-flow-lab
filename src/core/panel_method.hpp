// panel_method.hpp: the source panel method for non-lifting bodies of any shape (Anderson, §3.17).
//
// A body's surface is approximated by N straight panels, numbered clockwise, each carrying a
// source sheet of constant strength lambda_j (volume flow per unit length). The strengths are not
// given but solved for: at each panel's midpoint (its control point) the flow normal to the panel,
// from the rest of the flow and from every panel, must vanish,
//     lambda_i / 2 + sum_{j != i} lambda_j I_ij + V_ext . n_i = 0,       i = 1 .. N,
// where lambda_i / 2 is the panel's own normal velocity just outside it and I_ij the normal velocity
// at control point i of panel j with unit strength. All bodies of a scene are solved together, so
// they feel each other, and the rest of the flow (uniform flows, sources, vortices, doublets, their
// wall images) enters through V_ext. The surface is then a streamline (exactly at the control
// points), and the surface speed is the tangential velocity there.
//
// One panel from z1 to z2, of length L and unit direction t (a complex number), gives
//     W(z) = (lambda / 2 pi) conj(t) log((z - z1) / (z - z2)),
// which only jumps across the panel itself. Summed over a closed body, the panels' logarithms
// regroup by vertex v_k:
//     W(z) = (1 / 2 pi) sum_k c_k log((z - v_k) / (z - z_c)),   c_k = lambda_k conj(t_k) - lambda_{k-1} conj(t_{k-1}),
// (the c_k add up to zero) with z_c the body's centre. Each logarithm is cut only along the segment
// from v_k to the centre, inside the body, so W is continuous everywhere in the fluid and at the
// control points it is the limit from outside. The complex potential follows the same way:
//     F(z) = (1 / 2 pi) sum_k c_k (z - v_k) log((z - v_k) / (z - z_c)) + (Q / 2 pi) log(z - z_c),
// with Q = sum lambda_j L_j the net source strength of the body (zero for the exact solution of a
// closed body; the panels get close to it, which the Verify checks measure).
// Far from the body (more than twice its radius from z_c) W is summed as its Laurent series in
// 1 / (z - z_c), to full double precision but much faster than N logarithms.
//
// Inside the body the field is not physical (the panel sheets make their own interior flow), so
// the flow domain excludes it (FlowField::inDomain), and the plot fills it in.
//
// Without circulation a body has no lift; a vortex placed inside one gives it circulation
// (around a circle that is the lifting cylinder) but nothing fixes how much: that needs vortex
// panels and the Kutta condition.
#pragma once

#include "core/flow_element.hpp"
#include "core/vec2.hpp"

#include <optional>
#include <span>
#include <string>
#include <vector>

namespace pfl {

struct FlowField;

enum class PanelShape { Ellipse, Rectangle };

inline constexpr int k_minPanels = 4;
inline constexpr int k_maxPanels = 400;

// A body to be modelled with source panels, as the user describes it.
struct PanelBody {
    PanelShape shape = PanelShape::Ellipse;
    Vec2 position;           // centre
    double semiAxisA = 1.0;  // a: half the length along the body's own x axis (a circle's radius when a = b)
    double semiAxisB = 1.0;  // b: half the length across it
    double angle = 0.0;      // turns the body counter-clockwise about its centre
    int panels = 40;         // N

    friend bool operator==(const PanelBody&, const PanelBody&) = default;
};

const char* shapeName(PanelShape shape);

// Why the body is invalid (axes not positive, too few or too many panels, an odd number on a
// rectangle, which has the same number on opposite sides), or nullopt.
std::optional<std::string> panelBodyError(const PanelBody& body);

// The N panel vertices, clockwise (so the body is on the right of each panel and perp(t) points out
// of it). Panel j runs from vertex j to vertex j + 1 (mod N).
// Ellipse: x = a cos(eta), y = b sin(eta) at eta_k = pi + pi/N - 2 pi k / N, equal steps in eta, which
// puts a control point on each end of the a axis and closes the panels up at the pointed ends of a
// slender ellipse. Rectangle (2a by 2b): the panels shared out between the sides in proportion to
// their length, at least one per side, cosine-spaced along each side so they are short at the corners.
std::vector<Vec2> panelVertices(const PanelBody& body);

// Even-odd rule for the polygon with these vertices (not repeated at the end).
bool insidePolygon(std::span<const Vec2> vertices, Vec2 p);

// The eccentric-anomaly parameter eta of each panel midpoint of an ellipse (for the exact solution).
std::vector<double> ellipseMidpointParameters(const PanelBody& body);

// A body after the solve: its panels and their strengths.
struct SourcePanels {
    std::size_t element = 0; // the scene element it was made from
    PanelBody body;
    std::vector<Vec2> vertices;      // N, clockwise
    std::vector<Vec2> controlPoints; // panel midpoints
    std::vector<Vec2> normals;       // outward unit normals
    std::vector<double> lengths;
    std::vector<double> strengths; // lambda_j

    Vec2 centre;          // where the vertex logarithms are cut to (the centroid)
    double radius = 0.0;  // largest distance of a vertex from the centre
    Rect box;             // bounding box of the vertices
    std::vector<Complex> vertexWeights; // c_k
    std::vector<Complex> laurent;       // a_m, m = 1 .. M: W = sum a_m / (z - z_c)^m beyond 2 radius
    double netSource = 0.0;             // Q = sum lambda_j L_j

    std::size_t size() const { return vertices.size(); }
    Complex complexVelocity(Vec2 p) const;
    Complex complexPotential(Vec2 p) const;
    bool contains(Vec2 p) const; // inside the panel polygon
    double perimeter() const;
};

// Solves every body at once in the flow of `external` (which must not contain them). The result is in
// the order of `bodies`, each tagged with its scene index from `elements`.
std::vector<SourcePanels> solveSourcePanels(std::span<const PanelBody> bodies, std::span<const std::size_t> elements,
                                            const FlowField& external);

// Largest |v . n| at the control points over the largest |v| there, from `field` (which holds the
// bodies): 0 up to rounding for a solved body.
double normalFlowResidual(const SourcePanels& body, const FlowField& field);

// The exact flow past an ellipse with semi-axes a (along its own x axis) and b, centred at `centre`
// and turned by `angle`, in a uniform stream `freeStream`, with circulation `circulation` (> 0
// clockwise): the flow past the circle |zeta| = R = (a + b) / 2 mapped by z = zeta + k^2 / zeta, k^2 =
// (a^2 - b^2) / 4. Outside a body the flow is fixed by the free stream and the circulation alone, so
// this is also the flow with a vortex anywhere inside the ellipse.
// Velocity at a point outside the ellipse.
Complex exactEllipseVelocity(double a, double b, Vec2 centre, double angle, Vec2 freeStream, double circulation, Vec2 p);
// Surface speed at the point of parameter eta (x = a cos eta, y = b sin eta in the body's frame):
//     V = |2 V_inf sin(eta - alpha) + Gamma / (2 pi R)| R / sqrt(a^2 sin^2 eta + b^2 cos^2 eta),
// with alpha the stream direction measured from the body's own x axis. For a = b it is the cylinder's.
double exactEllipseSurfaceSpeed(double a, double b, double angle, Vec2 freeStream, double circulation, double eta);
// Strength of the continuous source sheet on the ellipse that makes this flow outside it (no
// circulation). The sheet's potential is continuous, so inside it is the harmonic function with
// the surface values of the outside potential, 2 V_inf R cos(eta - alpha): a uniform flow
// (2 R q_x / a, 2 R q_y / b) in the body's frame, q the stream. The sheet strength is the jump in
// normal velocity, 0 outside minus the inside one: lambda = -v_inside . n. For a circle,
// lambda = -2 V_inf cos(theta - alpha), 2 V_inf on the point facing the stream.
double exactEllipseSourceStrength(double a, double b, double angle, Vec2 freeStream, double eta);

} // namespace pfl
