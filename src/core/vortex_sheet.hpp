// vortex_sheet.hpp: thin airfoils, modelled as in classical thin-airfoil theory by a vortex sheet on
// the camber line (Anderson, §4.7 to 4.8; Katz & Plotkin, §5 and §11).
//
// The element. A camber line (a flat plate, a NACA 4-digit mean line or a circular arc) of chord c,
// placed by the middle of its chord and turned counter-clockwise by theta about it (so in a stream
// along +x the angle of attack is alpha = -theta, as for the NACA source-panel body). It has no
// thickness: the flow goes past both sides of the line.
//
// The numerical method: linear-strength vortex panels. The camber line is cut into N straight
// panels with cosine spacing along the chord, x_k = (c/2)(1 - cos(k pi / N)), k = 0 .. N, numbered
// from the leading edge (k = 0) to the trailing edge (k = N). gamma (circulation per unit length,
// > 0 clockwise like Gamma) varies linearly along each panel between its values gamma_k at the
// vertices, so the sheet strength is continuous. The unknowns are gamma_0 .. gamma_(N-1); the
// trailing-edge value is gamma_N = 0, which is the Kutta condition gamma(TE) = 0 stated directly.
// The N equations say that no fluid crosses the camber line at the panel midpoints (the control
// points). They are solved in one dense system with every source-panel body of the scene
// (panel_method.hpp), against the rest of the flow, wall images included, so airfoils and bodies
// feel each other.
//
// Why these panels: constant-strength vortex panels collocated at their midpoints are symmetric
// fore and aft and carry no Kutta condition (on a flat plate they converge to the solution without
// circulation), and the lumped-vortex 1/4-3/4 rule uses point vortices, whose velocity spikes
// near the line spoil the streamlines and particles there. Linear panels keep the velocity finite
// along the line except at the leading edge, where the exact gamma is infinite (like 1/sqrt(x))
// and the panels give a weak logarithmic singularity instead; the cosine spacing puts short panels
// there. Their solution approaches the exact flow past the infinitely thin line as N grows, with
// the boundary condition on the line itself (not linearized onto the chord), so for a flat plate
// c_l -> 2 pi sin(alpha), the exact value, and for a circular arc the conformal-map result.
//
// One panel from z1 to z2 (unit direction t, length L, gamma1 and gamma2 at its ends) gives, with
// zeta = conj(t) (z - z1) and Lambda = log((z - z1) / (z - z2)),
//     W(z) = (i conj(t) / 2 pi) [Lambda (gamma1 (1 - zeta / L) + gamma2 zeta / L) - (gamma2 - gamma1)],
// and above a panel along +x u = +gamma / 2, below u = -gamma / 2: the flow over the top is
// faster where gamma > 0. Lambda is cut only along the panel itself, so the velocity jumps only
// across the sheet, as it should. Exactly on the sheet (the control points) the velocity is the
// mean of the two sides (the principal value), and the two sides are that mean +-(gamma / 2) t.
// Far from the sheet (beyond twice its radius from the middle of the chord) W is summed as a
// Laurent series in 1 / (z - z_c), and nearer in each group of 16 panels uses a series of its own
// where the point is far enough from it, as for the source-panel bodies (panel_expansion.hpp).
//
// The complex potential F is multivalued: phi changes by Gamma around the airfoil. Its branch cuts
// are the sheet itself (phi really jumps across a vortex sheet) and a ray from the trailing edge
// straight on along the chord, downstream at small incidence: phi jumps by Gamma across that ray
// (the wake line of the Kutta condition). Written panel by panel, with the logarithms of each panel
// telescoped along the sheet to the trailing edge so no other cut enters the fluid:
//     F = (i / 2 pi) [ sum_m (A_m(z) + Gamma_<m) Lambda_m + P_m(z) + Gamma log_w(z - z_TE) ],
// A_m = gamma1 zeta + g zeta^2 / 2, P_m = -gamma1 L - g zeta L / 2 - g L^2 / 4, g = (gamma2 - gamma1) / L,
// Gamma_<m the circulation of the panels ahead of panel m and log_w the logarithm cut along the ray.
//
// Thin-airfoil theory (thinAirfoilTheory()) linearizes the same problem: the boundary condition
// moved onto the chord, dz/dx small, alpha small. With x = (c/2)(1 - cos theta0),
//     gamma(theta) = 2 V_inf [A0 (1 + cos theta) / sin theta + sum_n An sin n theta],
//     A0 = alpha - (1/pi) int_0^pi dz/dx dtheta0,   An = (2/pi) int_0^pi dz/dx cos n theta0 dtheta0,
//     c_l = 2 pi (A0 + A1 / 2) = 2 pi (alpha - alpha_L0),   c_m,c/4 = (pi/4)(A2 - A1),
//     alpha_L0 = -(1/pi) int_0^pi dz/dx (cos theta0 - 1) dtheta0,   Delta Cp = 2 gamma / V_inf.
// The integrals are in closed form for the flat plate and the NACA mean line, and summed by the
// trapezoidal rule in theta0 (exponentially accurate: dz/dx is smooth and periodic in theta0) for
// the circular arc. The panels and the theory agree to first order in camber and incidence; what is
// left is the theory's linearization (2 pi alpha against 2 pi sin alpha on a plate).
#pragma once

#include "core/flow_element.hpp"
#include "core/panel_expansion.hpp"
#include "core/vec2.hpp"

#include <optional>
#include <string>
#include <vector>

namespace pfl {

struct FlowField;

enum class CamberLine { Flat, Naca, Arc };

inline constexpr int k_minSheetPanels = 4;
inline constexpr int k_maxSheetPanels = 400;

// A thin airfoil as the user describes it.
struct ThinAirfoil {
    CamberLine camber = CamberLine::Naca;
    Vec2 position;      // the middle of the chord
    double chord = 2.0; // c
    double angle = 0.0; // theta: turns it counter-clockwise about the middle of the chord
    int panels = 40;    // N
    // NACA m p: the camber m (% of the chord, negative upside down) at p tenths of the chord.
    int nacaCamber = 2;
    int nacaPosition = 4;
    // Circular arc: its height h above the chord at mid-chord, as a fraction h / c (negative: sagging).
    double arcCamber = 0.05;

    friend bool operator==(const ThinAirfoil&, const ThinAirfoil&) = default;
};

// "flat", "naca", "arc": the names used in scene files.
const char* camberName(CamberLine camber);

// "flat plate", "NACA 24xx mean line", "circular arc h/c = 5%".
std::string airfoilName(const ThinAirfoil& airfoil);

// Why the airfoil is invalid (sizes, panel count, camber out of range), or nullopt.
std::optional<std::string> thinAirfoilError(const ThinAirfoil& airfoil);

// The camber line of a unit chord: its height z and slope dz/dx at x in [0, 1] from the leading
// edge, in the airfoil's own frame (chord along x).
double camberHeight(const ThinAirfoil& airfoil, double x);
double camberSlope(const ThinAirfoil& airfoil, double x);

// The N + 1 panel vertices in the plot, from the leading edge to the trailing edge.
std::vector<Vec2> sheetVertices(const ThinAirfoil& airfoil);

// Unit vector along the chord, from the leading edge to the trailing edge.
Vec2 chordDirection(const ThinAirfoil& airfoil);

// The angle of attack in a stream `freeStream`: its direction relative to the chord, in (-pi, pi].
double angleOfAttack(const ThinAirfoil& airfoil, Vec2 freeStream);

// --- Thin-airfoil theory ------------------------------------------------------------------

inline constexpr int k_theoryTerms = 100; // A0 .. A100, for the gamma(x) curve

struct ThinAirfoilTheory {
    double alpha = 0.0;           // angle of attack (rad)
    std::vector<double> A;        // A0 .. A_k_theoryTerms
    double alphaZeroLift = 0.0;   // alpha_L0 (rad)
    double cl = 0.0;              // 2 pi (A0 + A1 / 2)
    double cmQuarter = 0.0;       // (pi / 4)(A2 - A1)

    // Delta Cp = Cp_lower - Cp_upper = 2 gamma / V_inf at x/c (infinite at the leading edge unless A0 = 0).
    double deltaCp(double x) const;
};

ThinAirfoilTheory thinAirfoilTheory(const ThinAirfoil& airfoil, double alpha);

// --- The solved sheet ---------------------------------------------------------------------

// The velocity of one linear-strength vortex panel from a to b at p: W = gammaA * first + gammaB *
// second, with gammaA and gammaB its strengths at a and at b. On the panel itself: the mean of the
// two sides.
struct PanelInfluence {
    Complex first;
    Complex second;
};
PanelInfluence unitVortexPanelVelocity(Vec2 a, Vec2 b, Vec2 p);

struct VortexSheet {
    std::size_t element = 0; // the scene element it was made from
    ThinAirfoil airfoil;
    std::vector<Vec2> vertices;      // N + 1, leading edge to trailing edge
    std::vector<Vec2> controlPoints; // N panel midpoints
    std::vector<Vec2> tangents;      // N unit vectors, towards the trailing edge
    std::vector<Vec2> normals;       // N unit vectors perp(tangent): towards the upper side
    std::vector<double> lengths;     // N
    std::vector<double> chordwise;   // N + 1: x/c of each vertex, along the chord from the leading edge
    std::vector<double> gamma;       // N + 1 strengths at the vertices; gamma.back() = 0 (Kutta)

    double circulation = 0.0;            // Gamma = int gamma ds, > 0 clockwise
    std::vector<double> circulationAhead; // N: that of the panels ahead of each panel
    // Gamma for a unit stream along +x and along +y with the other panel bodies and airfoils there,
    // when nothing else is in the flow (only uniform flows): then Gamma(beta) = V (Gx cos beta + Gy
    // sin beta) for a stream V at beta, which gives alpha_L0 and the lift slope of the panels.
    std::optional<Vec2> unitCirculation;

    Vec2 centre;          // the middle of the chord
    double radius = 0.0;  // largest distance of a vertex from it
    Rect box;             // bounding box of the vertices
    std::vector<Complex> laurent; // as SourcePanels::laurent
    std::vector<PanelCluster> clusters; // groups of panels with series of their own (panel_expansion.hpp)

    std::size_t size() const { return lengths.size(); } // N
    Vec2 leadingEdge() const { return vertices.front(); }
    Vec2 trailingEdge() const { return vertices.back(); }
    double strengthAt(std::size_t panel) const { return 0.5 * (gamma[panel] + gamma[panel + 1]); } // at its control point

    Complex complexVelocity(Vec2 p) const;
    Complex complexPotential(Vec2 p) const;
    // The two sides of the sheet at control point `panel`: the mean velocity plus or minus (gamma/2) t.
    Vec2 upperVelocity(const FlowField& field, std::size_t panel) const;
    Vec2 lowerVelocity(const FlowField& field, std::size_t panel) const;
    // True if the segment from a to b crosses the sheet.
    bool crossedBy(Vec2 a, Vec2 b) const;
    double distanceTo(Vec2 p) const;
    double longestPanel() const;
};

// The geometry of a sheet (no strengths yet: all gamma 0).
VortexSheet sheetGeometry(const ThinAirfoil& airfoil, std::size_t element);

// Sets the strengths at the vertices (N + 1 values, the last 0) and everything that follows from them.
void setSheetStrengths(VortexSheet& sheet, std::vector<double> gamma);

// Largest |v . n| at the control points over the largest |v| there, from `field` (which holds the sheet).
double sheetNormalFlowResidual(const VortexSheet& sheet, const FlowField& field);

// --- What the sheet does: loading, lift and moment, against thin-airfoil theory -------------

struct AirfoilLoading {
    std::size_t sheet = 0;           // index in FlowField::vortexSheets
    double alpha = 0.0;              // angle of attack in the free stream
    double speed = 0.0;              // V_inf
    std::vector<double> x;           // x/c of the vertices
    std::vector<double> deltaCp;     // 2 gamma / V_inf at the vertices
    std::vector<Vec2> points;        // the vertices
    std::vector<double> cpUpper;     // Cp on each side at the control points (Bernoulli)
    std::vector<double> cpLower;
    double circulation = 0.0;        // Gamma
    double clKuttaJoukowski = 0.0;   // 2 Gamma / (V_inf c)
    double clPressure = 0.0;         // int_0^1 Delta Cp d(x/c), along the chord
    double cmQuarter = 0.0;          // -(2 / V_inf c^2) int gamma (x - c/4) dx, nose up > 0
    double trailingEdgeGamma = 0.0;  // gamma at the last control point / V_inf (gamma_N itself is 0)
    std::optional<double> alphaZeroLift; // of the panels, when the flow is uniform apart from the panels
    std::optional<double> liftSlope;     // d c_l / d alpha (per radian) of the panels, likewise
    ThinAirfoilTheory theory;
};

// Nullopt without a free stream (no reference speed).
std::optional<AirfoilLoading> airfoilLoading(const FlowField& field, std::size_t sheet, Vec2 freeStream);

} // namespace pfl
