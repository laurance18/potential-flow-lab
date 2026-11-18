#include "core/vortex_sheet.hpp"

#include "core/flow_field.hpp"
#include "core/panel_geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <limits>
#include <numbers>

namespace pfl {
namespace {

constexpr double pi = std::numbers::pi;
constexpr Complex k_i{0.0, 1.0};
constexpr double k_farFieldRatio = 0.5; // the Laurent series is used where radius / |z - z_c| is below this
constexpr int k_laurentTerms = 50;      // as for the source-panel bodies
constexpr int k_arcSamples = 1024;      // trapezoidal rule in theta0 for the circular arc's theory
constexpr double k_maxArcCamber = 0.25; // h/c: a quarter of the chord (a half circle would be 0.5)

Vec2 rotated(Vec2 v, double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    return {c * v.x - s * v.y, s * v.x + c * v.y};
}

// Gauss-Legendre nodes and weights on [0, 1], enough to integrate the Laurent moments (polynomials
// of degree k_laurentTerms along a panel) exactly.
constexpr int k_gaussPoints = k_laurentTerms / 2 + 1;
struct GaussRule {
    std::array<double, k_gaussPoints> node{};
    std::array<double, k_gaussPoints> weight{};
};

const GaussRule& gaussRule() {
    static const GaussRule rule = [] {
        GaussRule r;
        const int n = k_gaussPoints;
        for (int i = 0; i < n; ++i) {
            double x = std::cos(pi * (i + 0.75) / (n + 0.5)); // Newton on P_n from the usual first guess
            double derivative = 1.0;
            for (int iteration = 0; iteration < 100; ++iteration) {
                double p0 = 1.0;
                double p1 = x;
                for (int k = 2; k <= n; ++k) {
                    const double p2 = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k;
                    p0 = p1;
                    p1 = p2;
                }
                derivative = n * (x * p1 - p0) / (x * x - 1.0);
                const double step = p1 / derivative;
                x -= step;
                if (std::abs(step) < 1e-16) break;
            }
            r.node[static_cast<std::size_t>(i)] = 0.5 * (1.0 - x);
            r.weight[static_cast<std::size_t>(i)] = 1.0 / ((1.0 - x * x) * derivative * derivative); // 2 / (...) on [-1, 1], halved
        }
        return r;
    }();
    return rule;
}

// log((z - z1) / (z - z2)) of a panel of length L, from zeta = conj(t) (z - z1); on the panel itself
// (where the ratio is negative real) the principal value, the mean of the two sides.
Complex logRatio(Complex zeta, double length) {
    const Complex ratio = zeta / (zeta - length);
    if (ratio.real() < 0.0 && std::abs(ratio.imag()) <= 1e-12 * -ratio.real()) return {std::log(-ratio.real()), 0.0};
    return std::log(ratio);
}

bool crossesProperly(Vec2 a, Vec2 b, Vec2 c, Vec2 d) {
    const double o1 = cross(b - a, c - a);
    const double o2 = cross(b - a, d - a);
    const double o3 = cross(d - c, a - c);
    const double o4 = cross(d - c, b - c);
    return ((o1 > 0.0) != (o2 > 0.0)) && ((o3 > 0.0) != (o4 > 0.0));
}

double distanceToSegment(Vec2 p, Vec2 a, Vec2 b) {
    const Vec2 d = b - a;
    const double l2 = lengthSquared(d);
    const double t = l2 > 0.0 ? std::clamp(dot(p - a, d) / l2, 0.0, 1.0) : 0.0;
    return length(p - (a + d * t));
}

// int_a^b cos(k theta) d theta.
double cosIntegral(int k, double a, double b) {
    if (k == 0) return b - a;
    return (std::sin(k * b) - std::sin(k * a)) / k;
}

// int_0^pi dz/dx cos(n theta0) d theta0 for n = 0 .. terms.
std::vector<double> slopeMoments(const ThinAirfoil& airfoil, int terms) {
    std::vector<double> moments(static_cast<std::size_t>(terms) + 1, 0.0);
    const double m = airfoil.nacaCamber / 100.0;
    if (airfoil.camber == CamberLine::Naca && m != 0.0) {
        // dz/dx = K (cos theta0 - cos theta_p), K = m/p^2 ahead of the largest camber and m/(1-p)^2 behind it.
        const double p = airfoil.nacaPosition / 10.0;
        const double cp = 1.0 - 2.0 * p;
        const double thetaP = std::acos(cp);
        const double front = m / (p * p);
        const double back = m / ((1.0 - p) * (1.0 - p));
        const auto piece = [cp](int n, double a, double b) {
            // int (cos theta - cp) cos n theta = (1/2)(S_(n-1) + S_(n+1)) - cp S_n.
            return 0.5 * (cosIntegral(std::abs(n - 1), a, b) + cosIntegral(n + 1, a, b)) - cp * cosIntegral(n, a, b);
        };
        for (int n = 0; n <= terms; ++n) {
            moments[static_cast<std::size_t>(n)] = front * piece(n, 0.0, thetaP) + back * piece(n, thetaP, pi);
        }
    } else if (airfoil.camber == CamberLine::Arc && airfoil.arcCamber != 0.0) {
        // Trapezoidal rule: dz/dx(x(theta0)) is smooth, even and 2 pi periodic in theta0.
        std::vector<double> slope(k_arcSamples + 1);
        for (int k = 0; k <= k_arcSamples; ++k) {
            const double theta = pi * k / k_arcSamples;
            slope[static_cast<std::size_t>(k)] = camberSlope(airfoil, 0.5 * (1.0 - std::cos(theta)));
        }
        for (int n = 0; n <= terms; ++n) {
            double sum = 0.5 * (slope.front() + slope.back() * (n % 2 == 0 ? 1.0 : -1.0));
            for (int k = 1; k < k_arcSamples; ++k) sum += slope[static_cast<std::size_t>(k)] * std::cos(n * pi * k / k_arcSamples);
            moments[static_cast<std::size_t>(n)] = sum * pi / k_arcSamples;
        }
    }
    return moments;
}

} // namespace

const char* camberName(CamberLine camber) {
    switch (camber) {
    case CamberLine::Flat: return "flat";
    case CamberLine::Naca: return "naca";
    case CamberLine::Arc: return "arc";
    }
    return "?";
}

std::string airfoilName(const ThinAirfoil& airfoil) {
    switch (airfoil.camber) {
    case CamberLine::Flat: return "flat plate";
    case CamberLine::Naca:
        if (airfoil.nacaCamber == 0) return "flat plate (NACA 00xx mean line)";
        return std::format("{}NACA {}{}xx mean line", airfoil.nacaCamber < 0 ? "inverted " : "", std::abs(airfoil.nacaCamber),
                           airfoil.nacaPosition);
    case CamberLine::Arc:
        if (airfoil.arcCamber == 0.0) return "flat plate (circular arc h/c = 0)";
        return std::format("circular arc h/c = {:.3g}%", 100.0 * airfoil.arcCamber);
    }
    return "?";
}

std::optional<std::string> thinAirfoilError(const ThinAirfoil& airfoil) {
    if (!std::isfinite(airfoil.position.x) || !std::isfinite(airfoil.position.y) || !std::isfinite(airfoil.chord) ||
        !std::isfinite(airfoil.angle) || !std::isfinite(airfoil.arcCamber)) {
        return "All parameters must be finite numbers.";
    }
    if (!(airfoil.chord > 0.0)) return "The chord c must be positive.";
    if (airfoil.panels < k_minSheetPanels || airfoil.panels > k_maxSheetPanels) {
        return std::format("The number of panels N must be between {} and {}.", k_minSheetPanels, k_maxSheetPanels);
    }
    if (airfoil.camber == CamberLine::Naca) {
        if (airfoil.nacaCamber < -9 || airfoil.nacaCamber > 9) return "The camber m must be 0 to 9 % of the chord.";
        if (airfoil.nacaCamber != 0 && (airfoil.nacaPosition < 1 || airfoil.nacaPosition > 9)) {
            return "The position p of the camber must be 1 to 9 tenths of the chord.";
        }
    }
    if (airfoil.camber == CamberLine::Arc && std::abs(airfoil.arcCamber) > k_maxArcCamber) {
        return "The camber h/c of the arc must be at most 25%.";
    }
    return std::nullopt;
}

double camberHeight(const ThinAirfoil& airfoil, double x) {
    switch (airfoil.camber) {
    case CamberLine::Flat: return 0.0;
    case CamberLine::Naca:
        if (airfoil.nacaCamber == 0) return 0.0;
        return nacaCamber(airfoil.nacaCamber / 100.0, airfoil.nacaPosition / 10.0, x);
    case CamberLine::Arc: {
        const double h = airfoil.arcCamber;
        if (h == 0.0) return 0.0;
        // The circle through both ends and (1/2, h): radius R = (1/4 + h^2) / 2|h|. Written as
        // h - d^2 / (sqrt(R^2 - d^2) + R) so that a shallow arc loses no digits.
        const double radius = (0.25 + h * h) / (2.0 * std::abs(h));
        const double d = x - 0.5;
        return std::copysign(std::abs(h) - d * d / (std::sqrt(radius * radius - d * d) + radius), h);
    }
    }
    return 0.0;
}

double camberSlope(const ThinAirfoil& airfoil, double x) {
    switch (airfoil.camber) {
    case CamberLine::Flat: return 0.0;
    case CamberLine::Naca:
        if (airfoil.nacaCamber == 0) return 0.0;
        return nacaCamberSlope(airfoil.nacaCamber / 100.0, airfoil.nacaPosition / 10.0, x);
    case CamberLine::Arc: {
        const double h = airfoil.arcCamber;
        if (h == 0.0) return 0.0;
        const double radius = (0.25 + h * h) / (2.0 * std::abs(h));
        const double d = x - 0.5;
        return std::copysign(1.0, h) * -d / std::sqrt(radius * radius - d * d);
    }
    }
    return 0.0;
}

std::vector<Vec2> sheetVertices(const ThinAirfoil& airfoil) {
    const int n = airfoil.panels;
    std::vector<Vec2> vertices;
    vertices.reserve(static_cast<std::size_t>(n) + 1);
    for (int k = 0; k <= n; ++k) {
        const double x = k == n ? 1.0 : 0.5 * (1.0 - std::cos(pi * k / n));
        const Vec2 local{(x - 0.5) * airfoil.chord, camberHeight(airfoil, x) * airfoil.chord};
        vertices.push_back(airfoil.position + rotated(local, airfoil.angle));
    }
    return vertices;
}

Vec2 chordDirection(const ThinAirfoil& airfoil) { return {std::cos(airfoil.angle), std::sin(airfoil.angle)}; }

double angleOfAttack(const ThinAirfoil& airfoil, Vec2 freeStream) {
    const Vec2 e = chordDirection(airfoil);
    return std::atan2(cross(e, freeStream), dot(e, freeStream));
}

double ThinAirfoilTheory::deltaCp(double x) const {
    if (A.empty()) return 0.0;
    const double theta = std::acos(std::clamp(1.0 - 2.0 * x, -1.0, 1.0));
    const double s = std::sin(theta);
    if (!(s > 0.0)) {
        // The leading edge: infinite unless A0 = 0; the trailing edge: 0 (Kutta).
        if (theta < 1.0 && A[0] != 0.0) return std::copysign(std::numeric_limits<double>::infinity(), A[0]);
        return 0.0;
    }
    const double c = std::cos(theta);
    double sum = A[0] * (1.0 + c) / s;
    double previous = 0.0; // sin(0 theta)
    double current = s;    // sin(theta)
    for (std::size_t n = 1; n < A.size(); ++n) {
        sum += A[n] * current;
        const double next = 2.0 * c * current - previous;
        previous = current;
        current = next;
    }
    return 4.0 * sum;
}

ThinAirfoilTheory thinAirfoilTheory(const ThinAirfoil& airfoil, double alpha) {
    ThinAirfoilTheory theory;
    theory.alpha = alpha;
    const std::vector<double> moments = slopeMoments(airfoil, k_theoryTerms);
    theory.A.resize(moments.size());
    theory.A[0] = alpha - moments[0] / pi;
    for (std::size_t n = 1; n < moments.size(); ++n) theory.A[n] = 2.0 * moments[n] / pi;
    theory.alphaZeroLift = -(moments[1] - moments[0]) / pi;
    theory.cl = 2.0 * pi * (theory.A[0] + 0.5 * theory.A[1]);
    theory.cmQuarter = 0.25 * pi * (theory.A[2] - theory.A[1]);
    return theory;
}

PanelInfluence unitVortexPanelVelocity(Vec2 a, Vec2 b, Vec2 p) {
    const double l = length(b - a);
    const Complex t = toComplex(b - a) / l;
    const Complex zeta = std::conj(t) * toComplex(p - a);
    const Complex lambda = logRatio(zeta, l);
    const Complex factor = k_i * std::conj(t) / (2.0 * pi);
    return {factor * (lambda * (1.0 - zeta / l) + 1.0), factor * (lambda * zeta / l - 1.0)};
}

Complex VortexSheet::complexVelocity(Vec2 p) const {
    const Complex dz = toComplex(p - centre);
    const double distance = std::abs(dz);
    if (radius < k_farFieldRatio * distance) {
        const Complex w = radius / dz;
        const double ratio = radius / distance;
        const int terms = std::min(k_laurentTerms, static_cast<int>(std::ceil(-36.8 / std::log(ratio))) + 1);
        Complex sum{};
        for (int m = terms; m >= 1; --m) sum = (sum + laurent[static_cast<std::size_t>(m - 1)]) * w;
        return sum;
    }
    Complex sum{};
    for (std::size_t j = 0; j < lengths.size(); ++j) {
        const Complex t = toComplex(tangents[j]);
        const double l = lengths[j];
        const Complex zeta = std::conj(t) * toComplex(p - vertices[j]);
        const double g1 = gamma[j];
        const double g2 = gamma[j + 1];
        sum += std::conj(t) * (logRatio(zeta, l) * (g1 + (g2 - g1) * zeta / l) - (g2 - g1));
    }
    return k_i * sum / (2.0 * pi);
}

Complex VortexSheet::complexPotential(Vec2 p) const {
    Complex sum{};
    for (std::size_t j = 0; j < lengths.size(); ++j) {
        const Complex t = toComplex(tangents[j]);
        const double l = lengths[j];
        const Complex zeta = std::conj(t) * toComplex(p - vertices[j]);
        const double g1 = gamma[j];
        const double g = (gamma[j + 1] - g1) / l;
        const Complex a = g1 * zeta + 0.5 * g * zeta * zeta + circulationAhead[j];
        sum += a * logRatio(zeta, l) - g1 * l - 0.5 * g * zeta * l - 0.25 * g * l * l;
    }
    // The rest of each panel's logarithm, telescoped to the trailing edge: cut along the ray from it
    // straight on along the chord.
    sum += circulation * std::log(-std::polar(1.0, -airfoil.angle) * toComplex(p - trailingEdge()));
    return k_i * sum / (2.0 * pi);
}

Vec2 VortexSheet::upperVelocity(const FlowField& field, std::size_t panel) const {
    return field.velocity(controlPoints[panel]) + tangents[panel] * (0.5 * strengthAt(panel));
}

Vec2 VortexSheet::lowerVelocity(const FlowField& field, std::size_t panel) const {
    return field.velocity(controlPoints[panel]) - tangents[panel] * (0.5 * strengthAt(panel));
}

bool VortexSheet::crossedBy(Vec2 a, Vec2 b) const {
    if (std::max(a.x, b.x) < box.min.x || std::min(a.x, b.x) > box.max.x || std::max(a.y, b.y) < box.min.y ||
        std::min(a.y, b.y) > box.max.y) {
        return false;
    }
    for (std::size_t j = 0; j < lengths.size(); ++j) {
        if (crossesProperly(a, b, vertices[j], vertices[j + 1])) return true;
    }
    return false;
}

double VortexSheet::distanceTo(Vec2 p) const {
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t j = 0; j < lengths.size(); ++j) best = std::min(best, distanceToSegment(p, vertices[j], vertices[j + 1]));
    return best;
}

double VortexSheet::longestPanel() const { return lengths.empty() ? 0.0 : *std::ranges::max_element(lengths); }

VortexSheet sheetGeometry(const ThinAirfoil& airfoil, std::size_t element) {
    VortexSheet sheet;
    sheet.element = element;
    sheet.airfoil = airfoil;
    sheet.vertices = sheetVertices(airfoil);
    const std::size_t n = sheet.vertices.size() - 1;
    const Vec2 e = chordDirection(airfoil);
    const Vec2 leading = sheet.vertices.front();
    sheet.centre = airfoil.position;
    sheet.box = {leading, leading};
    for (std::size_t k = 0; k <= n; ++k) {
        const Vec2 v = sheet.vertices[k];
        sheet.chordwise.push_back(k == 0 ? 0.0 : k == n ? 1.0 : dot(v - leading, e) / airfoil.chord);
        sheet.radius = std::max(sheet.radius, length(v - sheet.centre));
        sheet.box.min = {std::min(sheet.box.min.x, v.x), std::min(sheet.box.min.y, v.y)};
        sheet.box.max = {std::max(sheet.box.max.x, v.x), std::max(sheet.box.max.y, v.y)};
        if (k == n) break;
        const Vec2 d = sheet.vertices[k + 1] - v;
        const double l = length(d);
        sheet.controlPoints.push_back(v + d * 0.5);
        sheet.lengths.push_back(l);
        sheet.tangents.push_back(d / l);
        sheet.normals.push_back(perp(d / l));
    }
    setSheetStrengths(sheet, std::vector<double>(n + 1, 0.0));
    return sheet;
}

void setSheetStrengths(VortexSheet& sheet, std::vector<double> gamma) {
    const std::size_t n = sheet.lengths.size();
    gamma.resize(n + 1, 0.0);
    gamma.back() = 0.0; // Kutta
    sheet.gamma = std::move(gamma);
    sheet.circulation = 0.0;
    sheet.circulationAhead.assign(n, 0.0);
    for (std::size_t j = 0; j < n; ++j) {
        sheet.circulationAhead[j] = sheet.circulation;
        sheet.circulation += sheet.lengths[j] * sheet.strengthAt(j);
    }

    // a_m / R^m = (i / 2 pi R) int gamma(s) u(s)^(m-1) ds, u = (z(s) - z_c) / R, by Gauss-Legendre on each panel.
    sheet.laurent.assign(k_laurentTerms, {});
    const GaussRule& rule = gaussRule();
    for (std::size_t j = 0; j < n; ++j) {
        const Complex u1 = toComplex(sheet.vertices[j] - sheet.centre) / sheet.radius;
        const Complex u2 = toComplex(sheet.vertices[j + 1] - sheet.centre) / sheet.radius;
        for (int q = 0; q < k_gaussPoints; ++q) {
            const double tau = rule.node[static_cast<std::size_t>(q)];
            const double strength = sheet.gamma[j] + (sheet.gamma[j + 1] - sheet.gamma[j]) * tau;
            const Complex u = u1 + (u2 - u1) * tau;
            Complex power = rule.weight[static_cast<std::size_t>(q)] * sheet.lengths[j] * strength;
            for (int m = 1; m <= k_laurentTerms; ++m) {
                sheet.laurent[static_cast<std::size_t>(m - 1)] += power;
                power *= u;
            }
        }
    }
    for (Complex& a : sheet.laurent) a *= k_i / (2.0 * pi * sheet.radius);
}

double sheetNormalFlowResidual(const VortexSheet& sheet, const FlowField& field) {
    double normal = 0.0;
    double scale = 0.0;
    for (std::size_t i = 0; i < sheet.size(); ++i) {
        const Vec2 v = field.velocity(sheet.controlPoints[i]);
        normal = std::max(normal, std::abs(dot(v, sheet.normals[i])));
        scale = std::max(scale, length(v));
    }
    return scale > 0.0 ? normal / scale : normal;
}

std::optional<AirfoilLoading> airfoilLoading(const FlowField& field, std::size_t index, Vec2 freeStream) {
    const double speed = length(freeStream);
    if (!(speed > 0.0) || index >= field.vortexSheets.size()) return std::nullopt;
    const VortexSheet& sheet = field.vortexSheets[index];
    const double c = sheet.airfoil.chord;
    const std::size_t n = sheet.size();

    AirfoilLoading loading;
    loading.sheet = index;
    loading.speed = speed;
    loading.alpha = angleOfAttack(sheet.airfoil, freeStream);
    loading.x = sheet.chordwise;
    loading.points = sheet.vertices;
    for (const double g : sheet.gamma) loading.deltaCp.push_back(2.0 * g / speed);
    for (std::size_t j = 0; j < n; ++j) {
        const double upper = length(sheet.upperVelocity(field, j)) / speed;
        const double lower = length(sheet.lowerVelocity(field, j)) / speed;
        loading.cpUpper.push_back(1.0 - upper * upper);
        loading.cpLower.push_back(1.0 - lower * lower);
    }

    loading.circulation = sheet.circulation;
    loading.clKuttaJoukowski = 2.0 * sheet.circulation / (speed * c);
    double moment = 0.0; // int gamma (x - c/4) dx
    for (std::size_t j = 0; j < n; ++j) {
        const double g1 = sheet.gamma[j];
        const double g2 = sheet.gamma[j + 1];
        const double dx = (loading.x[j + 1] - loading.x[j]) * c;
        loading.clPressure += 0.5 * (loading.deltaCp[j] + loading.deltaCp[j + 1]) * (loading.x[j + 1] - loading.x[j]);
        const double x1 = loading.x[j] * c - 0.25 * c;
        const double x2 = loading.x[j + 1] * c - 0.25 * c;
        moment += dx * (g1 * x1 + 0.5 * (g1 * (x2 - x1) + (g2 - g1) * x1) + (g2 - g1) * (x2 - x1) / 3.0);
    }
    loading.cmQuarter = -2.0 * moment / (speed * c * c);
    loading.trailingEdgeGamma = sheet.strengthAt(n - 1) / speed;

    if (sheet.unitCirculation) {
        // Gamma(beta) = V (Gx cos beta + Gy sin beta) = V R sin(beta + phi): no lift at beta = -phi.
        const Vec2 g = *sheet.unitCirculation;
        const double beta = std::atan2(freeStream.y, freeStream.x);
        double zero = std::remainder(-std::atan2(g.x, g.y) - sheet.airfoil.angle, 2.0 * pi);
        if (zero > 0.5 * pi) zero -= pi;
        if (zero <= -0.5 * pi) zero += pi;
        loading.alphaZeroLift = zero;
        loading.liftSlope = 2.0 * (-g.x * std::sin(beta) + g.y * std::cos(beta)) / c;
    }
    loading.theory = thinAirfoilTheory(sheet.airfoil, loading.alpha);
    return loading;
}

} // namespace pfl
