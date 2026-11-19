#include "core/panel_method.hpp"

#include "core/flow_field.hpp"
#include "core/panel_geometry.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace pfl {
namespace {

constexpr double pi = std::numbers::pi;

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

// Solves the dense system a x = b (row-major n x n) for `m` right-hand sides at once (b row-major
// n x m) by Gaussian elimination with partial pivoting. Returns x, n x m.
std::vector<double> solveDense(std::vector<double> a, std::vector<double> b, std::size_t m) {
    const std::size_t n = m > 0 ? b.size() / m : 0;
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t pivot = col;
        for (std::size_t row = col + 1; row < n; ++row) {
            if (std::abs(a[row * n + col]) > std::abs(a[pivot * n + col])) pivot = row;
        }
        if (pivot != col) {
            for (std::size_t k = 0; k < n; ++k) std::swap(a[col * n + k], a[pivot * n + k]);
            for (std::size_t r = 0; r < m; ++r) std::swap(b[col * m + r], b[pivot * m + r]);
        }
        const double diagonal = a[col * n + col];
        if (diagonal == 0.0) continue; // singular (overlapping bodies): leave that unknown at 0
        for (std::size_t row = col + 1; row < n; ++row) {
            const double f = a[row * n + col] / diagonal;
            if (f == 0.0) continue;
            for (std::size_t k = col; k < n; ++k) a[row * n + k] -= f * a[col * n + k];
            for (std::size_t r = 0; r < m; ++r) b[row * m + r] -= f * b[col * m + r];
        }
    }
    std::vector<double> x(n * m, 0.0);
    for (std::size_t i = n; i-- > 0;) {
        const double diagonal = a[i * n + i];
        for (std::size_t r = 0; r < m; ++r) {
            double sum = b[i * m + r];
            for (std::size_t k = i + 1; k < n; ++k) sum -= a[i * n + k] * x[k * m + r];
            x[i * m + r] = diagonal != 0.0 ? sum / diagonal : 0.0;
        }
    }
    return x;
}

bool crossesProperly(Vec2 a, Vec2 b, Vec2 c, Vec2 d) {
    const double o1 = cross(b - a, c - a);
    const double o2 = cross(b - a, d - a);
    const double o3 = cross(d - c, a - c);
    const double o4 = cross(d - c, b - c);
    return ((o1 > 0.0) != (o2 > 0.0)) && ((o3 > 0.0) != (o4 > 0.0)) && o1 != 0.0 && o2 != 0.0;
}

// Every vertex sees `centre` without crossing a panel (and `centre` is inside).
bool isStarShaped(const std::vector<Vec2>& v, Vec2 centre) {
    const std::size_t n = v.size();
    if (!insidePolygon(v, centre)) return false;
    for (std::size_t k = 0; k < n; ++k) {
        for (std::size_t j = 0; j < n; ++j) {
            if (j == k || (j + 1) % n == k) continue; // the two panels that end at vertex k
            if (crossesProperly(v[k], centre, v[j], v[(j + 1) % n])) return false;
        }
    }
    return true;
}

// The vertices in the body's own frame, before it is turned and moved.
std::vector<Vec2> localVertices(const PanelBody& body) {
    const int n = body.panels;
    const double a = body.semiAxisA;
    const double b = body.semiAxisB;
    std::vector<Vec2> local;
    local.reserve(static_cast<std::size_t>(n));
    switch (body.shape) {
    case PanelShape::Ellipse:
        for (int k = 0; k < n; ++k) {
            const double eta = pi + pi / n - 2.0 * pi * k / n;
            local.push_back({a * std::cos(eta), b * std::sin(eta)});
        }
        break;
    case PanelShape::Rectangle: {
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
        break;
    }
    case PanelShape::Naca:
        // Unit chord from the leading edge at (0, 0): scaled to the chord 2a, the middle of it at 0.
        for (const Vec2 p : nacaVertices(body.nacaCamber / 100.0, body.nacaPosition / 10.0, b / a, n)) {
            local.push_back({(p.x - 0.5) * 2.0 * a, p.y * 2.0 * a});
        }
        break;
    case PanelShape::Custom: {
        // Stretched from the outline's bounding box to [-a, a] x [-b, b].
        const std::vector<Vec2> outline = cleanOutline(body.outline);
        Vec2 lo = outline.front();
        Vec2 hi = outline.front();
        for (const Vec2 p : outline) {
            lo = {std::min(lo.x, p.x), std::min(lo.y, p.y)};
            hi = {std::max(hi.x, p.x), std::max(hi.y, p.y)};
        }
        const Vec2 mid = (lo + hi) * 0.5;
        std::vector<Vec2> scaled;
        scaled.reserve(outline.size());
        for (const Vec2 p : outline) scaled.push_back({(p.x - mid.x) / (0.5 * (hi.x - lo.x)) * a, (p.y - mid.y) / (0.5 * (hi.y - lo.y)) * b});
        local = distributeAlongOutline(scaled, body.smooth, n);
        break;
    }
    }
    return local;
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
    body.panelWeights.assign(n, {});
    for (std::size_t j = 0; j < n; ++j) {
        const Complex t = direction(body.vertices[j], body.vertices[(j + 1) % n]);
        const Complex weight = body.strengths[j] * std::conj(t);
        body.vertexWeights[j] += weight;           // panel j starts at vertex j
        body.vertexWeights[(j + 1) % n] -= weight; // and ends at vertex j + 1
        body.panelWeights[j] = weight / (2.0 * pi);
        body.netSource += body.strengths[j] * body.lengths[j];
    }
    body.starShaped = isStarShaped(body.vertices, body.centre);

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

    // The groups of panels: panel j, from a to b, has log((z - a)/(z - b)) = sum_m ((b - z_q)^m - (a - z_q)^m) / (m (z - z_q)^m).
    body.clusters = clusterPanels(n, [&](std::size_t j) { return std::pair{body.vertices[j], body.vertices[(j + 1) % n]}; });
    for (PanelCluster& c : body.clusters) {
        c.laurent.assign(k_laurentTerms, {});
        for (std::size_t j = c.first; j < c.last; ++j) {
            const Complex a = toComplex(body.vertices[j] - c.centre) / c.radius;
            const Complex b = toComplex(body.vertices[(j + 1) % n] - c.centre) / c.radius;
            Complex powerA = a;
            Complex powerB = b;
            for (int m = 1; m <= k_laurentTerms; ++m) {
                c.laurent[static_cast<std::size_t>(m - 1)] += body.panelWeights[j] * (powerB - powerA) / static_cast<double>(m);
                powerA *= a;
                powerB *= b;
            }
        }
    }
}

} // namespace

const char* shapeName(PanelShape shape) {
    switch (shape) {
    case PanelShape::Ellipse: return "ellipse";
    case PanelShape::Rectangle: return "rectangle";
    case PanelShape::Naca: return "naca";
    case PanelShape::Custom: return "outline";
    }
    return "?";
}

std::string nacaName(const PanelBody& body) {
    const int thickness = static_cast<int>(std::lround(100.0 * body.semiAxisB / body.semiAxisA));
    const int camber = std::abs(body.nacaCamber);
    const std::string digits = std::to_string(camber) + std::to_string(camber == 0 ? 0 : body.nacaPosition) +
                               (thickness < 10 ? "0" : "") + std::to_string(thickness);
    return (body.nacaCamber < 0 ? "inverted NACA " : "NACA ") + digits;
}

std::optional<std::string> panelBodyError(const PanelBody& body) {
    if (!std::isfinite(body.position.x) || !std::isfinite(body.position.y) || !std::isfinite(body.semiAxisA) ||
        !std::isfinite(body.semiAxisB) || !std::isfinite(body.angle)) {
        return "All parameters must be finite numbers.";
    }
    if (!(body.semiAxisA > 0.0) || !(body.semiAxisB > 0.0)) return "The sizes a and b must be positive.";
    if (body.panels < k_minPanels || body.panels > k_maxPanels) {
        return "The number of panels N must be between " + std::to_string(k_minPanels) + " and " +
               std::to_string(k_maxPanels) + ".";
    }
    switch (body.shape) {
    case PanelShape::Ellipse: break;
    case PanelShape::Rectangle:
        if (body.panels % 2 != 0) return "A rectangle needs an even number of panels N (the same number on opposite sides).";
        break;
    case PanelShape::Naca:
        if (body.panels % 2 != 0) return "An airfoil needs an even number of panels N (half on each surface).";
        if (body.nacaCamber < -9 || body.nacaCamber > 9) return "The camber m must be 0 to 9 % of the chord.";
        if (body.nacaCamber != 0 && (body.nacaPosition < 1 || body.nacaPosition > 9)) {
            return "The position p of the camber must be 1 to 9 tenths of the chord.";
        }
        if (body.semiAxisB > 0.4 * body.semiAxisA) return "The thickness t/c = b/a must be at most 40%.";
        break;
    case PanelShape::Custom: return outlineError(body.outline, body.smooth, body.panels);
    }
    return std::nullopt;
}

std::vector<Vec2> panelVertices(const PanelBody& body) {
    const std::vector<Vec2> local = localVertices(body);
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
    if (radius < k_farFieldRatio * distance) return laurentVelocity(laurent, radius, dz, distance);
    // Near the body: each group of panels by its own series where p is far enough from it, else panel
    // by panel. Each panel's logarithm is cut only along the panel itself; it is written out (ln |r|
    // from |r|^2, no complex division), since this loop is where particles and streamlines spend
    // their time near a body. On the panel itself (its control point) the ratio is negative real, and
    // the limit from outside, on the right of the clockwise panel, has the argument -pi.
    const std::size_t n = vertices.size();
    double re = 0.0;
    double im = 0.0;
    Complex series{};
    for (const PanelCluster& c : clusters) {
        const Vec2 d = p - c.centre;
        const double d2 = lengthSquared(d);
        if (useSeries(c.radius, d2)) {
            series += laurentVelocity(c.laurent, c.radius, toComplex(d), std::sqrt(d2));
            continue;
        }
        Vec2 u1 = p - vertices[c.first];
        double l1 = lengthSquared(u1);
        for (std::size_t j = c.first; j < c.last; ++j) {
            const Vec2 u2 = p - vertices[(j + 1) % n];
            const double l2 = lengthSquared(u2);
            // log(u1 / u2) = ln(|u1| / |u2|) + i arg(u1 conj(u2)).
            const double x = dot(u1, u2);
            const double y = cross(u2, u1);
            const double arg = x < 0.0 && std::abs(y) <= 1e-12 * -x ? -pi : std::atan2(y, x);
            const double lnr = 0.5 * std::log(l1 / l2);
            const Complex w = panelWeights[j];
            re += w.real() * lnr - w.imag() * arg;
            im += w.imag() * lnr + w.real() * arg;
            u1 = u2;
            l1 = l2;
        }
    }
    return Complex{re, im} + series;
}

Complex SourcePanels::complexPotential(Vec2 p) const {
    if (!starShaped) {
        // F_j = (lambda_j / 2 pi) (conj(t_j) (u1 log u1 - u2 log u2) - L_j), u = z - the panel's ends.
        Complex sum{};
        const std::size_t n = vertices.size();
        for (std::size_t j = 0; j < n; ++j) {
            const Complex u1 = toComplex(p - vertices[j]);
            const Complex u2 = toComplex(p - vertices[(j + 1) % n]);
            sum += panelWeights[j] * (u1 * std::log(u1) - u2 * std::log(u2)) - strengths[j] * lengths[j] / (2.0 * pi);
        }
        return sum;
    }
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
    return solvePanelMethods(bodies, elements, {}, {}, external).bodies;
}

PanelSolution solvePanelMethods(std::span<const PanelBody> bodies, std::span<const std::size_t> bodyElements,
                                std::span<const ThinAirfoil> airfoils, std::span<const std::size_t> airfoilElements,
                                const FlowField& external) {
    PanelSolution solution;
    std::vector<SourcePanels>& solved = solution.bodies;
    std::vector<VortexSheet>& sheets = solution.sheets;
    solved.resize(bodies.size());
    std::size_t total = 0;
    for (std::size_t b = 0; b < bodies.size(); ++b) {
        SourcePanels& s = solved[b];
        s.element = b < bodyElements.size() ? bodyElements[b] : b;
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
    // Unknowns of an airfoil: gamma at its vertices 0 .. N-1; gamma_N = 0 (Kutta).
    for (std::size_t a = 0; a < airfoils.size(); ++a) {
        sheets.push_back(sheetGeometry(airfoils[a], a < airfoilElements.size() ? airfoilElements[a] : a));
        total += sheets.back().size();
    }

    // Rows: the normal velocity at every control point, the bodies' first. Columns: a body panel with
    // unit source strength, or an airfoil vertex with unit gamma. Three right-hand sides: the flow
    // outside the panels, and unit streams along +x and +y (for the airfoils' alpha_L0 and lift slope).
    struct Row {
        Vec2 point;
        Complex normal;
        const SourcePanels* body;  // the body or sheet the control point is on, and its panel
        const VortexSheet* sheet;
        std::size_t panel;
    };
    std::vector<Row> rows;
    rows.reserve(total);
    for (const SourcePanels& body : solved) {
        for (std::size_t i = 0; i < body.size(); ++i) rows.push_back({body.controlPoints[i], toComplex(body.normals[i]), &body, nullptr, i});
    }
    for (const VortexSheet& sheet : sheets) {
        for (std::size_t i = 0; i < sheet.size(); ++i) {
            rows.push_back({sheet.controlPoints[i], toComplex(sheet.normals[i]), nullptr, &sheet, i});
        }
    }

    constexpr std::size_t sides = 3;
    std::vector<double> matrix(total * total, 0.0);
    std::vector<double> rhs(total * sides, 0.0);
    for (std::size_t row = 0; row < rows.size(); ++row) {
        const Row& r = rows[row];
        double* line = &matrix[row * total];
        // v . n = Re(conj(W) conj(n)) = Re(W n). A stream along +x has W = 1, along +y W = -i.
        rhs[row * sides + 0] = -(external.complexVelocity(r.point) * r.normal).real();
        rhs[row * sides + 1] = -r.normal.real();
        rhs[row * sides + 2] = -r.normal.imag();
        std::size_t column = 0;
        for (const SourcePanels& source : solved) {
            const std::size_t n = source.size();
            for (std::size_t j = 0; j < n; ++j, ++column) {
                if (&source == r.body && j == r.panel) {
                    line[column] = 0.5; // the panel's own sheet, just outside it
                } else {
                    const Complex w = unitPanelVelocity(source.vertices[j], source.vertices[(j + 1) % n], r.point);
                    line[column] = (w * r.normal).real();
                }
            }
        }
        for (const VortexSheet& sheet : sheets) {
            const std::size_t n = sheet.size();
            // On its own panel (r.sheet, r.panel) the influence is the principal value, which
            // unitVortexPanelVelocity() takes there.
            for (std::size_t j = 0; j < n; ++j) {
                const PanelInfluence w = unitVortexPanelVelocity(sheet.vertices[j], sheet.vertices[j + 1], r.point);
                line[column + j] += (w.first * r.normal).real();
                if (j + 1 < n) line[column + j + 1] += (w.second * r.normal).real();
            }
            column += n;
        }
    }

    const std::vector<double> x = solveDense(std::move(matrix), std::move(rhs), sides);
    const auto unknown = [&x](std::size_t k, std::size_t side) { return x[k * sides + side]; };
    std::size_t offset = 0;
    for (SourcePanels& s : solved) {
        s.strengths.resize(s.size());
        for (std::size_t j = 0; j < s.size(); ++j) s.strengths[j] = unknown(offset + j, 0);
        offset += s.size();
        finishBody(s);
    }
    // alpha_L0 and the lift slope of an airfoil follow from the unit streams only when the rest of
    // the flow is uniform too.
    const bool uniformOnly = external.rows.empty() && std::ranges::all_of(external.elements, [](const FlowElement& e) {
                                 return std::holds_alternative<UniformFlow>(e);
                             });
    for (VortexSheet& sheet : sheets) {
        const std::size_t n = sheet.size();
        std::vector<double> gamma(n + 1, 0.0);
        Vec2 unit;
        for (std::size_t k = 0; k < n; ++k) {
            gamma[k] = unknown(offset + k, 0);
            // The circulation of the unit-stream solutions, with the trapezoidal sum of linear gamma.
            const double weight = 0.5 * (sheet.lengths[k] + (k > 0 ? sheet.lengths[k - 1] : 0.0));
            unit += Vec2{unknown(offset + k, 1), unknown(offset + k, 2)} * weight;
        }
        offset += n;
        setSheetStrengths(sheet, std::move(gamma));
        if (uniformOnly) sheet.unitCirculation = unit;
    }
    return solution;
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
