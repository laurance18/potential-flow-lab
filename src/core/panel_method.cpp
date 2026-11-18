#include "core/panel_method.hpp"

#include "core/flow_field.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace pfl {
namespace {

constexpr double pi = std::numbers::pi;
constexpr double k_farFieldRatio = 0.5; // the Laurent series is used where radius / |z - z_c| is below this
constexpr int k_laurentTerms = 50;      // 0.5^50 ~ 1e-15: full precision at the switch-over distance

Vec2 rotated(Vec2 v, double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    return {c * v.x - s * v.y, s * v.x + c * v.y};
}

// Unit direction of the panel from a to b, as a complex number.
Complex direction(Vec2 a, Vec2 b) {
    const Complex d = toComplex(b - a);
    return d / std::abs(d);
}

// W at z of panel a -> b with unit strength.
Complex unitPanelVelocity(Vec2 a, Vec2 b, Vec2 z) {
    const Complex t = direction(a, b);
    return std::conj(t) / (2.0 * pi) * std::log((toComplex(z) - toComplex(a)) / (toComplex(z) - toComplex(b)));
}

Vec2 polygonCentroid(const std::vector<Vec2>& v) {
    const Vec2 origin = v.front();
    double twice = 0.0;
    Vec2 sum;
    for (std::size_t k = 0; k < v.size(); ++k) {
        const Vec2 a = v[k] - origin;
        const Vec2 b = v[(k + 1) % v.size()] - origin;
        const double c = cross(a, b);
        twice += c;
        sum += (a + b) * c;
    }
    return twice != 0.0 ? origin + sum / (3.0 * twice) : origin;
}

// Solves the dense system a x = b (row-major n x n) by Gaussian elimination with partial pivoting.
std::vector<double> solveDense(std::vector<double> a, std::vector<double> b) {
    const std::size_t n = b.size();
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < n; ++row) {
            if (std::abs(a[row * n + col]) > std::abs(a[pivot * n + col])) pivot = row;
        }
        if (pivot != col) {
            for (std::size_t k = 0; k < n; ++k) std::swap(a[col * n + k], a[pivot * n + k]);
            std::swap(b[col], b[pivot]);
        }
        const double diagonal = a[col * n + col];
        if (diagonal == 0.0) continue; // singular (overlapping bodies): leave that unknown at 0
        for (std::size_t row = col + 1; row < n; ++row) {
            const double f = a[row * n + col] / diagonal;
            if (f == 0.0) continue;
            for (std::size_t k = col; k < n; ++k) a[row * n + k] -= f * a[col * n + k];
            b[row] -= f * b[col];
        }
    }
    std::vector<double> x(n, 0.0);
    for (std::size_t i = n; i-- > 0;) {
        double sum = b[i];
        for (std::size_t k = i + 1; k < n; ++k) sum -= a[i * n + k] * x[k];
        const double diagonal = a[i * n + i];
        x[i] = diagonal != 0.0 ? sum / diagonal : 0.0;
    }
    return x;
}

// Fills in everything that follows from the vertices and the strengths.
void finishBody(SourcePanels& body) {
    const std::size_t n = body.size();
    body.centre = polygonCentroid(body.vertices);
    body.radius = 0.0;
    body.box = {body.vertices.front(), body.vertices.front()};
    for (const Vec2 v : body.vertices) {
        body.radius = std::max(body.radius, length(v - body.centre));
        body.box.min = {std::min(body.box.min.x, v.x), std::min(body.box.min.y, v.y)};
        body.box.max = {std::max(body.box.max.x, v.x), std::max(body.box.max.y, v.y)};
    }

    body.netSource = 0.0;
    body.vertexWeights.assign(n, {});
    for (std::size_t j = 0; j < n; ++j) {
        const Complex t = direction(body.vertices[j], body.vertices[(j + 1) % n]);
        const Complex weight = body.strengths[j] * std::conj(t);
        body.vertexWeights[j] += weight;           // panel j starts at vertex j
        body.vertexWeights[(j + 1) % n] -= weight; // and ends at vertex j + 1
        body.netSource += body.strengths[j] * body.lengths[j];
    }

    // a_m / radius^m = -(1 / 2 pi m) sum_k c_k ((v_k - z_c) / radius)^m: scaled so no power overflows.
    body.laurent.assign(k_laurentTerms, {});
    std::vector<Complex> power(n, Complex{1.0, 0.0});
    for (int m = 1; m <= k_laurentTerms; ++m) {
        Complex sum{};
        for (std::size_t k = 0; k < n; ++k) {
            power[k] *= toComplex(body.vertices[k] - body.centre) / body.radius;
            sum += body.vertexWeights[k] * power[k];
        }
        body.laurent[static_cast<std::size_t>(m - 1)] = -sum / (2.0 * pi * m);
    }
}

} // namespace

const char* shapeName(PanelShape shape) {
    switch (shape) {
    case PanelShape::Ellipse: return "ellipse";
    case PanelShape::Rectangle: return "rectangle";
    }
    return "?";
}

std::optional<std::string> panelBodyError(const PanelBody& body) {
    if (!std::isfinite(body.position.x) || !std::isfinite(body.position.y) || !std::isfinite(body.semiAxisA) ||
        !std::isfinite(body.semiAxisB) || !std::isfinite(body.angle)) {
        return "All parameters must be finite numbers.";
    }
    if (!(body.semiAxisA > 0.0) || !(body.semiAxisB > 0.0)) return "The semi-axes a and b must be positive.";
    if (body.panels < k_minPanels || body.panels > k_maxPanels) {
        return "The number of panels N must be between " + std::to_string(k_minPanels) + " and " +
               std::to_string(k_maxPanels) + ".";
    }
    if (body.shape == PanelShape::Rectangle && body.panels % 2 != 0) {
        return "A rectangle needs an even number of panels N (the same number on opposite sides).";
    }
    return std::nullopt;
}

std::vector<Vec2> panelVertices(const PanelBody& body) {
    const int n = body.panels;
    const double a = body.semiAxisA;
    const double b = body.semiAxisB;
    std::vector<Vec2> local;
    local.reserve(static_cast<std::size_t>(n));
    if (body.shape == PanelShape::Ellipse) {
        for (int k = 0; k < n; ++k) {
            const double eta = pi + pi / n - 2.0 * pi * k / n;
            local.push_back({a * std::cos(eta), b * std::sin(eta)});
        }
    } else {
        // Up the left side, along the top, down the right side, back along the bottom.
        const int half = n / 2;
        const int along = std::clamp(static_cast<int>(std::lround(half * a / (a + b))), 1, half - 1);
        const int across = half - along;
        const Vec2 corners[] = {{-a, -b}, {-a, b}, {a, b}, {a, -b}};
        const int counts[] = {across, along, across, along};
        for (int side = 0; side < 4; ++side) {
            const Vec2 from = corners[side];
            const Vec2 to = corners[(side + 1) % 4];
            for (int i = 0; i < counts[side]; ++i) {
                const double s = 0.5 * (1.0 - std::cos(pi * i / counts[side]));
                local.push_back(from + (to - from) * s);
            }
        }
    }
    std::vector<Vec2> vertices;
    vertices.reserve(local.size());
    for (const Vec2 p : local) vertices.push_back(body.position + rotated(p, body.angle));
    return vertices;
}

std::vector<double> ellipseMidpointParameters(const PanelBody& body) {
    std::vector<double> eta;
    if (body.shape != PanelShape::Ellipse) return eta;
    for (int k = 0; k < body.panels; ++k) eta.push_back(pi - 2.0 * pi * k / body.panels);
    return eta;
}

Complex SourcePanels::complexVelocity(Vec2 p) const {
    const Complex dz = toComplex(p - centre);
    const double distance = std::abs(dz);
    if (radius < k_farFieldRatio * distance) {
        // Horner in w = radius / (z - z_c); terms fall off like (radius / |z - z_c|)^m.
        const Complex w = radius / dz;
        const double ratio = radius / distance;
        const int terms = std::min(k_laurentTerms, static_cast<int>(std::ceil(-36.8 / std::log(ratio))) + 1);
        Complex sum{};
        for (int m = terms; m >= 1; --m) sum = (sum + laurent[static_cast<std::size_t>(m - 1)]) * w;
        return sum;
    }
    // log of each ratio written out (ln |r| from |r|^2, so no square root): this loop is where the
    // particles and streamlines spend their time near a body.
    const Complex q = 1.0 / dz;
    double re = 0.0;
    double im = 0.0;
    for (std::size_t k = 0; k < vertices.size(); ++k) {
        const double ux = p.x - vertices[k].x;
        const double uy = p.y - vertices[k].y;
        const double rx = ux * q.real() - uy * q.imag();
        const double ry = ux * q.imag() + uy * q.real();
        const double lnr = 0.5 * std::log(rx * rx + ry * ry);
        const double arg = std::atan2(ry, rx);
        const Complex c = vertexWeights[k];
        re += c.real() * lnr - c.imag() * arg;
        im += c.imag() * lnr + c.real() * arg;
    }
    return Complex{re, im} / (2.0 * pi);
}

Complex SourcePanels::complexPotential(Vec2 p) const {
    const Complex dz = toComplex(p - centre);
    const Complex q = 1.0 / dz;
    Complex sum{};
    for (std::size_t k = 0; k < vertices.size(); ++k) {
        const Complex u = toComplex(p - vertices[k]);
        sum += vertexWeights[k] * u * std::log(u * q);
    }
    return (sum + netSource * std::log(dz)) / (2.0 * pi);
}

bool insidePolygon(std::span<const Vec2> vertices, Vec2 p) {
    bool inside = false;
    for (std::size_t k = 0, j = vertices.size() - 1; k < vertices.size(); j = k++) {
        const Vec2 a = vertices[k];
        const Vec2 b = vertices[j];
        if ((a.y > p.y) != (b.y > p.y) && p.x < a.x + (p.y - a.y) * (b.x - a.x) / (b.y - a.y)) inside = !inside;
    }
    return inside;
}

bool SourcePanels::contains(Vec2 p) const { return box.contains(p) && insidePolygon(vertices, p); }

double SourcePanels::perimeter() const {
    double sum = 0.0;
    for (const double l : lengths) sum += l;
    return sum;
}

std::vector<SourcePanels> solveSourcePanels(std::span<const PanelBody> bodies, std::span<const std::size_t> elements,
                                            const FlowField& external) {
    std::vector<SourcePanels> solved(bodies.size());
    std::size_t total = 0;
    for (std::size_t b = 0; b < bodies.size(); ++b) {
        SourcePanels& s = solved[b];
        s.element = b < elements.size() ? elements[b] : b;
        s.body = bodies[b];
        s.vertices = panelVertices(bodies[b]);
        const std::size_t n = s.size();
        for (std::size_t j = 0; j < n; ++j) {
            const Vec2 a = s.vertices[j];
            const Vec2 c = s.vertices[(j + 1) % n];
            s.controlPoints.push_back((a + c) * 0.5);
            s.lengths.push_back(length(c - a));
            s.normals.push_back(perp(c - a) / length(c - a));
        }
        total += n;
    }

    // Row i: normal velocity at control point i; column j: panel j with unit strength.
    std::vector<double> matrix(total * total, 0.0);
    std::vector<double> rhs(total, 0.0);
    std::size_t row = 0;
    for (const SourcePanels& target : solved) {
        for (std::size_t i = 0; i < target.size(); ++i, ++row) {
            const Vec2 point = target.controlPoints[i];
            const Complex normal = toComplex(target.normals[i]);
            // v . n = Re(conj(W) conj(n)) = Re(W n).
            rhs[row] = -(external.complexVelocity(point) * normal).real();
            std::size_t column = 0;
            for (const SourcePanels& source : solved) {
                const std::size_t n = source.size();
                for (std::size_t j = 0; j < n; ++j, ++column) {
                    if (&source == &target && j == i) {
                        matrix[row * total + column] = 0.5; // the panel's own sheet, just outside it
                    } else {
                        const Complex w = unitPanelVelocity(source.vertices[j], source.vertices[(j + 1) % n], point);
                        matrix[row * total + column] = (w * normal).real();
                    }
                }
            }
        }
    }

    const std::vector<double> strengths = solveDense(std::move(matrix), std::move(rhs));
    std::size_t offset = 0;
    for (SourcePanels& s : solved) {
        s.strengths.assign(strengths.begin() + static_cast<std::ptrdiff_t>(offset),
                           strengths.begin() + static_cast<std::ptrdiff_t>(offset + s.size()));
        offset += s.size();
        finishBody(s);
    }
    return solved;
}

double normalFlowResidual(const SourcePanels& body, const FlowField& field) {
    double normal = 0.0;
    double scale = 0.0;
    for (std::size_t i = 0; i < body.size(); ++i) {
        const Complex w = field.complexVelocity(body.controlPoints[i]);
        normal = std::max(normal, std::abs((w * toComplex(body.normals[i])).real()));
        scale = std::max(scale, std::abs(w));
    }
    return scale > 0.0 ? normal / scale : normal;
}

Complex exactEllipseVelocity(double a, double b, Vec2 centre, double angle, Vec2 freeStream, double circulation,
                             Vec2 p) {
    const Complex z = toComplex(rotated(p - centre, -angle));
    const Complex stream = std::conj(toComplex(rotated(freeStream, -angle))); // W of the free stream, body frame
    const double k2 = 0.25 * (a * a - b * b);
    const double r2 = 0.25 * (a + b) * (a + b);
    const Complex root = std::sqrt(z * z - 4.0 * k2);
    const Complex zeta1 = 0.5 * (z + root);
    const Complex zeta2 = 0.5 * (z - root);
    const Complex zeta = std::abs(zeta1) >= std::abs(zeta2) ? zeta1 : zeta2; // outside the circle |zeta| = R
    const Complex dFdZeta =
        stream - std::conj(stream) * r2 / (zeta * zeta) + Complex{0.0, circulation / (2.0 * pi)} / zeta;
    const Complex dzdZeta = 1.0 - k2 / (zeta * zeta);
    return dFdZeta / dzdZeta * std::polar(1.0, -angle);
}

double exactEllipseSurfaceSpeed(double a, double b, double angle, Vec2 freeStream, double circulation, double eta) {
    const Vec2 local = rotated(freeStream, -angle);
    const double alpha = std::atan2(local.y, local.x);
    const double radius = 0.5 * (a + b);
    const double s = std::sin(eta);
    const double c = std::cos(eta);
    const double tangential = 2.0 * length(freeStream) * std::sin(eta - alpha) + circulation / (2.0 * pi * radius);
    return std::abs(tangential) * radius / std::sqrt(a * a * s * s + b * b * c * c);
}

double exactEllipseSourceStrength(double a, double b, double angle, Vec2 freeStream, double eta) {
    const Vec2 q = rotated(freeStream, -angle);
    const double radius = 0.5 * (a + b);
    const Vec2 inside{2.0 * radius * q.x / a, 2.0 * radius * q.y / b};
    const Vec2 normal = Vec2{std::cos(eta) / a, std::sin(eta) / b} / length(Vec2{std::cos(eta) / a, std::sin(eta) / b});
    return -dot(inside, normal);
}

} // namespace pfl
