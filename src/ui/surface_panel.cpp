#define IMGUI_DEFINE_MATH_OPERATORS
#include "ui/surface_panel.hpp"

#include "ui/style.hpp"

#include <implot.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <format>
#include <limits>
#include <numbers>
#include <string>

namespace pfl {
namespace {

constexpr ImU32 k_exactColor = IM_COL32(205, 205, 200, 230);
constexpr ImU32 k_zeroColor = IM_COL32(150, 150, 150, 110);
constexpr float k_pickRadiusPx = 14.0f;
constexpr double nan = std::numeric_limits<double>::quiet_NaN();

const char* axisLabel(int axis) {
    constexpr const char* labels[] = {"θ: angle around the body from the free stream (°)", "x/c along the free stream",
                                      "s: arc length from the stagnation point"};
    return labels[axis];
}

// One body's curve on the chosen axis; `index` maps each value back to the surface point, and NaN
// entries break the line (where the angle wraps from 360° to 0°).
struct Curve {
    std::vector<double> x;
    std::vector<double> cp;
    std::vector<double> exact; // empty unless the body is a circular cylinder
    std::vector<std::size_t> index;
};

Curve makeCurve(const BodySurface& body, int axis, const std::optional<CylinderReference>& cylinder) {
    Curve curve;
    const bool panelExact = !cylinder && !body.cpExact.empty();
    const std::vector<double>& source = axis == 0 ? body.angle : axis == 1 ? body.streamwise : body.arcLength;
    for (std::size_t k = 0; k < body.points.size(); ++k) {
        if (axis == 0 && k > 0 && std::abs(source[k] - source[k - 1]) > 180.0) {
            curve.x.push_back(nan);
            curve.cp.push_back(nan);
            if (cylinder || panelExact) curve.exact.push_back(nan);
            curve.index.push_back(k);
        }
        curve.x.push_back(source[k]);
        curve.cp.push_back(body.cp[k]);
        if (cylinder) {
            const Vec2 r = body.points[k] - cylinder->centre;
            curve.exact.push_back(cylinder->cp(std::atan2(r.y, r.x)));
        } else if (panelExact) {
            curve.exact.push_back(body.cpExact[k]);
        }
        curve.index.push_back(k);
    }
    return curve;
}

std::string bodyName(std::size_t b, const BodySurface& body) {
    return std::format("Body {}{}", b + 1, body.panelBody ? " (panels)" : body.closed ? "" : " (half-body)");
}

// The label of a body's closed-form curve.
const char* exactLabel(const BodySurface& body, const FlowModel& model) {
    if (!body.panelBody) return "Exact (cylinder)";
    const PanelBody& shape = model.field.panelBodies[*body.panelBody].body;
    return shape.semiAxisA == shape.semiAxisB ? "Exact (circle)" : "Exact (ellipse)";
}

// What the panel method says about a panel body.
void drawPanelInfo(const BodySurface& body, const FlowModel& model) {
    const SourcePanels& panels = model.field.panelBodies[*body.panelBody];
    double worst = 0.0;
    if (!body.cpExact.empty()) {
        for (std::size_t k = 0; k < body.cp.size(); ++k) worst = std::max(worst, std::abs(body.cp[k] - body.cpExact[k]));
    }
    ImGui::Text("Source panels: N = %zu, Σλs / (V∞ · perimeter) = %.2e", panels.size(),
                std::abs(panels.netSource) / (model.freeStreamSpeed * panels.perimeter()));
    ImGui::SetItemTooltip("The net source strength of the panels: zero for the exact solution of a closed body,\n"
                          "so a measure of the discretization error. It shrinks as N grows.");
    if (!body.cpExact.empty()) ImGui::Text("max |Cp − Cp exact| at the control points = %.2e", worst);
}

void forceRow(const char* name, const std::string& computed, const std::string& exact) {
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::TextUnformatted(name);
    ImGui::TableNextColumn();
    ImGui::TextUnformatted(computed.c_str());
    ImGui::TableNextColumn();
    ImGui::TextUnformatted(exact.c_str());
}

void drawForces(std::size_t b, const BodySurface& body, const FlowModel& model) {
    const BodyForces f = bodyForces(body, model);
    const double v = model.freeStreamSpeed;
    const auto minimum = std::ranges::min_element(body.cp);
    const Vec2 at = body.points[static_cast<std::size_t>(minimum - body.cp.begin())];

    ImGui::Spacing();
    ImGui::Text("%s: chord c = %.4g, Γ inside = %.4g", bodyName(b, body).c_str(), f.chord, f.circulation);
    ImGui::Text("Cp min = %.4g at (%.3g, %.3g)", *minimum, at.x, at.y);
    if (body.panelBody) drawPanelInfo(body, model);
    if (ImGui::BeginTable(("##forces" + std::to_string(b)).c_str(), 3,
                          ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_BordersInnerH)) {
        ImGui::TableSetupColumn("per unit span");
        ImGui::TableSetupColumn("from Cp");
        ImGui::TableSetupColumn(f.isolated ? "exact (isolated body)" : "exact");
        ImGui::TableHeadersRow();
        const bool exact = f.isolated;
        forceRow("Drag D'/q∞", std::format("{:.2e}", f.drag), exact ? "0 (d'Alembert)" : "-");
        forceRow("Lift L'/q∞", std::format("{:.5g}", f.lift),
                 exact ? std::format("2Γ/V∞ = {:.5g} (Kutta-Joukowski)", 2.0 * f.circulation / v) : "-");
        forceRow("c_l = L'/(q∞ c)", std::format("{:.5g}", f.lift / f.chord),
                 exact ? std::format("{:.5g}", 2.0 * f.circulation / (v * f.chord)) : "-");
        ImGui::EndTable();
    }
    if (!f.isolated) {
        ImGui::TextWrapped("Other singular points lie outside this body (other bodies, wall images or sector "
                           "flows), so the isolated-body values do not apply.");
    }
    if (body.closureGap > 0.0) {
        ImGui::TextWrapped("Nearly closed: the upper and lower surfaces miss the rear stagnation point by %.2g "
                           "(%.1f%% of the chord) because other bodies make the flow around this one asymmetric. "
                           "A short straight segment closes the surface there.",
                           body.closureGap, 100.0 * body.closureGap / f.chord);
    }
    if (body.stagnationPoints.size() == 1) {
        ImGui::TextWrapped("This surface is the dividing streamline through the single stagnation point. It "
                           "loops around the body together with the fluid circulating with it; the force is "
                           "the same as on the body inside.");
    }
}

// The source strength of every panel of the panel bodies along the surface, on the axis of the Cp
// plot, over the exact sheet where it is known. Returns the surface point hovered, if any.
std::optional<Vec2> drawStrengthPlot(const std::vector<BodySurface>& bodies, int axis, const FlowModel& model) {
    const float s = uiScale();
    const double speed = model.freeStreamSpeed;
    std::optional<Vec2> picked;
    ImGui::SeparatorText("Source strength of the panels");
    ImGui::SetItemTooltip("λ: volume flow per unit length out of each panel (> 0 a source, < 0 a sink), solved so that\n"
                          "no fluid crosses the panel at its control point. The exact curve is the continuous source\n"
                          "sheet that makes the same flow outside the body: λ = −v·n of the uniform flow inside it.");
    const ImVec2 plotSize(-1.0f, std::max(150.0f * s, 0.42f * ImGui::GetContentRegionAvail().y));
    constexpr ImPlotFlags plotFlags = ImPlotFlags_NoTitle | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect;
    if (!ImPlot::BeginPlot("##panel-strength", plotSize, plotFlags)) return picked;
    ImPlot::SetupAxis(ImAxis_X1, axis == 0 ? "θ (°)" : axis == 1 ? "x/c" : "s", ImPlotAxisFlags_AutoFit);
    ImPlot::SetupAxis(ImAxis_Y1, "λ / V∞", ImPlotAxisFlags_AutoFit);
    if (axis == 0) ImPlot::SetupAxisTicks(ImAxis_X1, 0.0, 360.0, 5);
    ImPlot::SetupLegend(ImPlotLocation_South, ImPlotLegendFlags_Outside | ImPlotLegendFlags_Horizontal);
    const double zero = 0.0;
    ImPlot::PlotInfLines("##lambda0", &zero, 1,
                         {ImPlotProp_LineColor, k_zeroColor, ImPlotProp_Flags,
                          static_cast<int>(ImPlotItemFlags_NoLegend) | static_cast<int>(ImPlotInfLinesFlags_Horizontal)});

    const ImVec2 mouse = ImGui::GetMousePos();
    float best = k_pickRadiusPx * s;
    std::optional<std::pair<std::size_t, std::size_t>> hit; // body, point
    for (std::size_t b = 0; b < bodies.size() && b < k_bodyColors.size(); ++b) {
        const BodySurface& body = bodies[b];
        if (!body.panelBody || body.strength.size() != body.points.size()) continue;
        const std::vector<double>& source = axis == 0 ? body.angle : axis == 1 ? body.streamwise : body.arcLength;
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> exact;
        std::vector<std::size_t> index;
        const bool withExact = body.strengthExact.size() == body.points.size();
        for (std::size_t k = 0; k < body.points.size(); ++k) {
            if (axis == 0 && k > 0 && std::abs(source[k] - source[k - 1]) > 180.0) {
                x.push_back(nan);
                y.push_back(nan);
                if (withExact) exact.push_back(nan);
                index.push_back(k);
            }
            x.push_back(source[k]);
            y.push_back(body.strength[k] / speed);
            if (withExact) exact.push_back(body.strengthExact[k] / speed);
            index.push_back(k);
        }
        const int n = static_cast<int>(x.size());
        if (withExact) {
            ImPlot::PlotLine("Exact sheet", x.data(), exact.data(), n,
                             {ImPlotProp_LineColor, k_exactColor, ImPlotProp_LineWeight, 5.0f * s});
        }
        ImPlot::PlotLine(bodyName(b, body).c_str(), x.data(), y.data(), n,
                         {ImPlotProp_LineColor, k_bodyColors[b], ImPlotProp_LineWeight, 2.0f * s, ImPlotProp_Marker,
                          ImPlotMarker_Circle, ImPlotProp_MarkerSize, 2.5f * s});
        if (!ImPlot::IsPlotHovered()) continue;
        for (std::size_t i = 0; i < x.size(); ++i) {
            if (std::isnan(x[i])) continue;
            const ImVec2 d = ImPlot::PlotToPixels(x[i], y[i]) - mouse;
            const float distance = std::sqrt(d.x * d.x + d.y * d.y);
            if (distance < best) {
                best = distance;
                hit = {b, index[i]};
            }
        }
    }
    if (hit) {
        const auto [b, k] = *hit;
        const BodySurface& body = bodies[b];
        const double x = axis == 0 ? body.angle[k] : axis == 1 ? body.streamwise[k] : body.arcLength[k];
        picked = body.points[k];
        ImDrawList* dl = ImPlot::GetPlotDrawList();
        ImPlot::PushPlotClipRect();
        dl->AddCircle(ImPlot::PlotToPixels(x, body.strength[k] / speed), 6.0f * s, k_bodyColors[b], 0, 2.0f * s);
        ImPlot::PopPlotClipRect();
        const std::size_t n = body.points.size() - 1;
        std::string text = std::format("{}, panel {} of {}\n{} = {:.4g}\nλ / V∞ = {:.5f}", bodyName(b, body),
                                       (body.firstPanel + k) % n + 1, n, axis == 0 ? "θ" : axis == 1 ? "x/c" : "s", x,
                                       body.strength[k] / speed);
        if (body.strengthExact.size() == body.points.size()) {
            text += std::format("  (exact sheet {:.5f})", body.strengthExact[k] / speed);
        }
        text += std::format("\nλ / 2πV∞ = {:.4f}", body.strength[k] / (2.0 * std::numbers::pi * speed));
        ImGui::SetTooltip("%s", text.c_str());
    }
    ImPlot::EndPlot();
    return picked;
}

// True if the airfoil is alone in a uniform stream, which is what thin-airfoil theory describes.
bool aloneInStream(const FlowModel& model) {
    const FlowField& field = model.field;
    return field.vortexSheets.size() == 1 && field.panelBodies.empty() && field.rows.empty() &&
           std::ranges::all_of(field.elements, [](const FlowElement& e) { return std::holds_alternative<UniformFlow>(e); });
}

void tableRow(const char* name, const std::string& panels, const std::string& theory, const char* help) {
    ImGui::TableNextRow();
    ImGui::TableNextColumn();
    ImGui::TextUnformatted(name);
    if (help != nullptr) ImGui::SetItemTooltip("%s", help);
    ImGui::TableNextColumn();
    ImGui::TextUnformatted(panels.c_str());
    ImGui::TableNextColumn();
    ImGui::TextUnformatted(theory.c_str());
}

// One airfoil: Delta Cp along the chord against thin-airfoil theory, the theory's Fourier
// coefficients, and lift, moment and alpha_L0 from the panels beside the theory. Returns the vertex
// hovered in the plot, if any.
std::optional<Vec2> drawAirfoil(std::size_t a, const AirfoilLoading& loading, const FlowModel& model) {
    const float s = uiScale();
    const double degrees = 180.0 / std::numbers::pi;
    std::optional<Vec2> picked;
    const VortexSheet& sheet = model.field.vortexSheets[loading.sheet];
    const ThinAirfoilTheory& theory = loading.theory;
    const std::string id = std::to_string(a);

    ImGui::SeparatorText(std::format("Airfoil {}: {} (#{})", a + 1, airfoilName(sheet.airfoil), sheet.element + 1).c_str());
    ImGui::Text("c = %.4g, α = %.4g°, N = %zu vortex panels, V∞ = %.4g", sheet.airfoil.chord, loading.alpha * degrees, sheet.size(),
                loading.speed);

    // The theory curve, finely sampled in theta (dense at both edges), and the plot's range from
    // everything but the leading-edge peak, where both go to infinity.
    std::vector<double> theoryX;
    std::vector<double> theoryY;
    constexpr int samples = 400;
    double lo = 0.0;
    double hi = 0.0;
    for (int k = 1; k <= samples; ++k) {
        const double x = 0.5 * (1.0 - std::cos(std::numbers::pi * k / samples));
        theoryX.push_back(x);
        theoryY.push_back(theory.deltaCp(x));
        if (x >= 0.03) {
            lo = std::min(lo, theoryY.back());
            hi = std::max(hi, theoryY.back());
        }
    }
    for (std::size_t k = 0; k < loading.x.size(); ++k) {
        if (loading.x[k] < 0.03) continue;
        lo = std::min(lo, loading.deltaCp[k]);
        hi = std::max(hi, loading.deltaCp[k]);
    }
    const double margin = 0.12 * std::max(hi - lo, 1e-3);

    const ImVec2 plotSize(-1.0f, std::max(160.0f * s, 0.3f * ImGui::GetContentRegionAvail().y));
    constexpr ImPlotFlags plotFlags = ImPlotFlags_NoTitle | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect;
    if (ImPlot::BeginPlot(("##dcp" + id).c_str(), plotSize, plotFlags)) {
        ImPlot::SetupAxis(ImAxis_X1, "x/c along the chord");
        ImPlot::SetupAxis(ImAxis_Y1, "ΔCp = 2γ/V∞");
        ImPlot::SetupAxisLimits(ImAxis_X1, -0.02, 1.02, ImPlotCond_Always);
        ImPlot::SetupAxisLimits(ImAxis_Y1, lo - margin, hi + margin, ImPlotCond_Always);
        ImPlot::SetupLegend(ImPlotLocation_NorthEast);
        const double zero = 0.0;
        ImPlot::PlotInfLines("##dcp0", &zero, 1,
                             {ImPlotProp_LineColor, k_zeroColor, ImPlotProp_Flags,
                              static_cast<int>(ImPlotItemFlags_NoLegend) | static_cast<int>(ImPlotInfLinesFlags_Horizontal)});
        ImPlot::PlotLine("Thin-airfoil theory", theoryX.data(), theoryY.data(), static_cast<int>(theoryX.size()),
                         {ImPlotProp_LineColor, k_exactColor, ImPlotProp_LineWeight, 5.0f * s});
        const ImU32 colour = k_bodyColors[a % k_bodyColors.size()];
        ImPlot::PlotLine("Vortex panels", loading.x.data(), loading.deltaCp.data(), static_cast<int>(loading.x.size()),
                         {ImPlotProp_LineColor, colour, ImPlotProp_LineWeight, 2.0f * s, ImPlotProp_Marker, ImPlotMarker_Circle,
                          ImPlotProp_MarkerSize, 2.5f * s});
        if (ImPlot::IsPlotHovered()) {
            const ImVec2 mouse = ImGui::GetMousePos();
            float best = k_pickRadiusPx * s;
            std::optional<std::size_t> hit;
            for (std::size_t k = 0; k < loading.x.size(); ++k) {
                const ImVec2 d = ImPlot::PlotToPixels(loading.x[k], loading.deltaCp[k]) - mouse;
                const float distance = std::sqrt(d.x * d.x + d.y * d.y);
                if (distance < best) {
                    best = distance;
                    hit = k;
                }
            }
            if (hit) {
                const std::size_t k = *hit;
                picked = loading.points[k];
                ImDrawList* dl = ImPlot::GetPlotDrawList();
                ImPlot::PushPlotClipRect();
                dl->AddCircle(ImPlot::PlotToPixels(loading.x[k], loading.deltaCp[k]), 6.0f * s, colour, 0, 2.0f * s);
                ImPlot::PopPlotClipRect();
                std::string text = std::format("Airfoil {}, vertex {} of {}\nx/c = {:.4f}\nΔCp = 2γ/V∞ = {:.5f}", a + 1, k + 1,
                                               loading.x.size(), loading.x[k], loading.deltaCp[k]);
                if (k > 0) text += std::format("  (theory {:.5f})", theory.deltaCp(loading.x[k]));
                if (k == 0) text += "\nThe leading edge: infinite in theory, finite on the panels.";
                if (k + 1 == loading.x.size()) text += "\nThe trailing edge: 0 by the Kutta condition.";
                ImGui::SetTooltip("%s", text.c_str());
            }
        }
        ImPlot::EndPlot();
    }

    // The theory's Fourier coefficients.
    if (ImGui::BeginTable(("##fourier" + id).c_str(), 4, ImGuiTableFlags_SizingStretchSame | ImGuiTableFlags_BordersInnerH)) {
        for (const char* name : {"A₀", "A₁", "A₂", "A₃"}) ImGui::TableSetupColumn(name);
        ImGui::TableHeadersRow();
        ImGui::TableNextRow();
        for (std::size_t n = 0; n < 4; ++n) {
            ImGui::TableNextColumn();
            ImGui::Text("%.5f", theory.A[n]);
        }
        ImGui::EndTable();
    }
    ImGui::SetItemTooltip("Thin-airfoil theory: γ(θ) = 2V∞[A₀(1 + cos θ)/sin θ + Σ Aₙ sin nθ], x = (c/2)(1 − cos θ),\n"
                          "A₀ = α − (1/π)∫ dz/dx dθ₀, Aₙ = (2/π)∫ dz/dx cos nθ₀ dθ₀ (0 to π). A₀ holds the angle of attack,\n"
                          "the others the camber only.");

    // Panels beside theory.
    if (ImGui::BeginTable(("##airfoil" + id).c_str(), 3, ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_BordersInnerH)) {
        ImGui::TableSetupColumn("");
        ImGui::TableSetupColumn("vortex panels");
        ImGui::TableSetupColumn("thin-airfoil theory");
        ImGui::TableHeadersRow();
        const double c = sheet.airfoil.chord;
        tableRow("Γ / (V∞ c)", std::format("{:.5f}", loading.circulation / (loading.speed * c)),
                 std::format("π(A₀ + A₁/2) = {:.5f}", std::numbers::pi * (theory.A[0] + 0.5 * theory.A[1])),
                 "Circulation Γ = ∫γ ds along the camber line, > 0 clockwise.");
        tableRow("c_l = 2Γ/(V∞c)", std::format("{:.5f}", loading.clKuttaJoukowski), std::format("2π(A₀ + A₁/2) = {:.5f}", theory.cl),
                 "Kutta–Joukowski: L′ = ρV∞Γ, so c_l = L′/(q∞c) = 2Γ/(V∞c).");
        tableRow("c_l = ∫ΔCp d(x/c)", std::format("{:.5f}", loading.clPressure), std::format("{:.5f}", theory.cl),
                 "The pressure difference integrated along the chord, as in thin-airfoil theory. γ is per unit length of\n"
                 "the camber line and dx = ds cos ε (ε its slope), so with camber this is a little below 2Γ/(V∞c): the two\n"
                 "differ by terms of order ε², which the theory leaves out. On a flat plate they agree exactly.");
        tableRow("c_m,c/4", std::format("{:.5f}", loading.cmQuarter + 0.0), std::format("(π/4)(A₂ − A₁) = {:.5f}", theory.cmQuarter + 0.0),
                 "Moment about the quarter chord, nose up > 0: −(2/V∞c²)∫γ (x − c/4) dx.");
        // + 0.0 turns -0 into 0.
        tableRow("α_L0", loading.alphaZeroLift ? std::format("{:.4f}°", *loading.alphaZeroLift * degrees + 0.0) : std::string("—"),
                 std::format("{:.4f}°", theory.alphaZeroLift * degrees + 0.0),
                 "Zero-lift angle. Panels: where Γ = 0, from the solutions for unit streams along x and y (only when the rest\n"
                 "of the flow is uniform). Theory: −(1/π)∫ dz/dx (cos θ₀ − 1) dθ₀.");
        tableRow("dc_l/dα", loading.liftSlope ? std::format("{:.4f} /rad", *loading.liftSlope) : std::string("—"), "2π = 6.2832 /rad",
                 "Lift slope at this α. For a flat plate the exact potential flow gives 2π cos α (c_l = 2π sin α).");
        tableRow("γ at the trailing edge", std::format("0 (last panel {:.2g} V∞)", loading.trailingEdgeGamma), "0",
                 "The Kutta condition: the flow leaves the trailing edge smoothly, γ(TE) = 0. The panels set the vertex value to 0;\n"
                 "at the last control point γ falls like √(c − x), to 0 as N grows.");
        ImGui::EndTable();
    }
    if (!aloneInStream(model)) {
        ImGui::TextWrapped("Other elements or bodies are in the flow; thin-airfoil theory is for the airfoil alone in a uniform "
                           "stream, so the two columns need not agree.");
    }
    return picked;
}

} // namespace

float surfacePanelWidth(float available) {
    const float s = uiScale();
    return std::clamp(0.34f * available, 320.0f * s, 560.0f * s);
}

std::optional<Vec2> SurfacePanel::draw(ImVec2 pos, ImVec2 size, const std::vector<BodySurface>& bodies,
                                       const std::vector<AirfoilLoading>& airfoils, const FlowModel& model) {
    ImGui::SetNextWindowPos(pos);
    ImGui::SetNextWindowSize(size);
    constexpr ImGuiWindowFlags flags = ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse |
                                       ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoDocking;
    std::optional<Vec2> picked;
    if (!ImGui::Begin("Surface pressure", nullptr, flags)) {
        ImGui::End();
        return picked;
    }
    const float s = uiScale();

    if (!(model.freeStreamSpeed > 0.0)) {
        ImGui::TextWrapped("Surface Cp needs a free stream. Add a uniform flow: Cp = 1 - |V|²/V∞² uses its speed "
                           "as the reference.");
    } else if (bodies.empty() && airfoils.empty()) {
        ImGui::TextWrapped("No body found. Bodies are bounded by the dividing streamlines through stagnation "
                           "points: add a doublet (cylinder), a source (half-body) or a source and a sink "
                           "(Rankine oval) to the uniform flow.\n\nA body is only found while its dividing "
                           "streamlines close (within 10%%). Asymmetric flows can leave them open: a source and "
                           "sink at incidence, for example, let part of the source's fluid escape downstream.");
    } else {
        for (std::size_t a = 0; a < airfoils.size(); ++a) {
            if (const std::optional<Vec2> point = drawAirfoil(a, airfoils[a], model)) picked = point;
        }
        if (!airfoils.empty() && !bodies.empty()) ImGui::Spacing();
    }
    if (model.freeStreamSpeed > 0.0 && !bodies.empty()) {
        const bool anyClosed = std::ranges::any_of(bodies, &BodySurface::closed);
        const bool anyOpen = !std::ranges::all_of(bodies, &BodySurface::closed);
        int axis = anyClosed ? static_cast<int>(m_axis) : 2;
        ImGui::SetNextItemWidth(-FLT_MIN);
        if (ImGui::BeginCombo("##axis", axisLabel(axis))) {
            for (int a = 0; a < 3; ++a) {
                if (ImGui::Selectable(axisLabel(a), a == axis, anyClosed || a == 2 ? 0 : ImGuiSelectableFlags_Disabled)) {
                    m_axis = static_cast<Axis>(a);
                    axis = a;
                }
            }
            ImGui::EndCombo();
        }

        const bool anyPanels = std::ranges::any_of(bodies, [](const BodySurface& b) { return b.panelBody.has_value(); });
        const ImVec2 plotSize(-1.0f, std::max(anyPanels ? 180.0f * s : 220.0f * s,
                                              (anyPanels ? 0.38f : 0.55f) * ImGui::GetContentRegionAvail().y));
        constexpr ImPlotFlags plotFlags = ImPlotFlags_NoTitle | ImPlotFlags_NoMenus | ImPlotFlags_NoBoxSelect;
        if (ImPlot::BeginPlot("##surface-cp", plotSize, plotFlags)) {
            ImPlot::SetupAxis(ImAxis_X1, axis == 0 ? "θ (°)" : axis == 1 ? "x/c" : "s", ImPlotAxisFlags_AutoFit);
            // Inverted: suction (Cp < 0) up, the usual way to plot pressure distributions.
            ImPlot::SetupAxis(ImAxis_Y1, "Cp", ImPlotAxisFlags_Invert | ImPlotAxisFlags_AutoFit);
            if (axis == 0) ImPlot::SetupAxisTicks(ImAxis_X1, 0.0, 360.0, 5);
            ImPlot::SetupLegend(ImPlotLocation_South, ImPlotLegendFlags_Outside | ImPlotLegendFlags_Horizontal);

            const double zero = 0.0;
            ImPlot::PlotInfLines("##cp0", &zero, 1,
                                 {ImPlotProp_LineColor, k_zeroColor, ImPlotProp_Flags,
                                  static_cast<int>(ImPlotItemFlags_NoLegend) | static_cast<int>(ImPlotInfLinesFlags_Horizontal)});

            const ImVec2 mouse = ImGui::GetMousePos();
            float best = k_pickRadiusPx * s;
            std::optional<std::pair<std::size_t, std::size_t>> hit; // body, point
            for (std::size_t b = 0; b < bodies.size() && b < k_bodyColors.size(); ++b) {
                const BodySurface& body = bodies[b];
                if (!body.closed && axis != 2) continue;
                const std::optional<CylinderReference> cylinder = cylinderReference(body, model);
                const Curve curve = makeCurve(body, axis, cylinder);
                const int n = static_cast<int>(curve.x.size());
                if (!curve.exact.empty()) {
                    ImPlot::PlotLine(exactLabel(body, model), curve.x.data(), curve.exact.data(), n,
                                     {ImPlotProp_LineColor, k_exactColor, ImPlotProp_LineWeight, 5.0f * s});
                }
                // Panel bodies: a marker at every control point, where the method gives Cp.
                ImPlot::PlotLine(bodyName(b, body).c_str(), curve.x.data(), curve.cp.data(), n,
                                 {ImPlotProp_LineColor, k_bodyColors[b], ImPlotProp_LineWeight, 2.0f * s,
                                  ImPlotProp_Marker, body.panelBody ? ImPlotMarker_Circle : ImPlotMarker_None,
                                  ImPlotProp_MarkerSize, 2.5f * s});
                if (!ImPlot::IsPlotHovered()) continue;
                for (std::size_t i = 0; i < curve.x.size(); ++i) {
                    if (std::isnan(curve.x[i])) continue;
                    const ImVec2 d = ImPlot::PlotToPixels(curve.x[i], curve.cp[i]) - mouse;
                    const float distance = std::sqrt(d.x * d.x + d.y * d.y);
                    if (distance < best) {
                        best = distance;
                        hit = {b, curve.index[i]};
                    }
                }
            }

            if (hit) {
                const auto [b, k] = *hit;
                const BodySurface& body = bodies[b];
                const double x = axis == 0 ? body.angle[k] : axis == 1 ? body.streamwise[k] : body.arcLength[k];
                picked = body.points[k];
                ImDrawList* dl = ImPlot::GetPlotDrawList();
                ImPlot::PushPlotClipRect();
                dl->AddCircle(ImPlot::PlotToPixels(x, body.cp[k]), 6.0f * s, k_bodyColors[b], 0, 2.0f * s);
                ImPlot::PopPlotClipRect();
                std::string text = std::format("{}\n{} = {:.4g}\nCp = {:.4f}", bodyName(b, body), axis == 0 ? "θ" : axis == 1 ? "x/c" : "s", x, body.cp[k]);
                if (const std::optional<CylinderReference> cylinder = cylinderReference(body, model)) {
                    const Vec2 r = body.points[k] - cylinder->centre;
                    text += std::format("  (exact {:.4f})", cylinder->cp(std::atan2(r.y, r.x)));
                } else if (!body.cpExact.empty()) {
                    text += std::format("  (exact {:.4f})", body.cpExact[k]);
                }
                text += std::format("\nat ({:.4g}, {:.4g})", body.points[k].x, body.points[k].y);
                ImGui::SetTooltip("%s", text.c_str());
            }
            ImPlot::EndPlot();
        }
        if (anyOpen && axis != 2) ImGui::TextDisabled("Half-bodies are open: they are plotted against s only.");
        if (bodies.size() > k_bodyColors.size()) ImGui::TextDisabled("Only the first %zu bodies are plotted.", k_bodyColors.size());
        if (anyPanels) {
            if (const std::optional<Vec2> point = drawStrengthPlot(bodies, axis, model)) picked = point;
        }

        if (anyClosed) {
            ImGui::SeparatorText("Forces from the surface pressure");
            for (std::size_t b = 0; b < bodies.size(); ++b) {
                if (bodies[b].closed) drawForces(b, bodies[b], model);
            }
            ImGui::Spacing();
            ImGui::TextDisabled("F/q∞ = -∫ Cp n ds around the surface (trapezoidal rule), q∞ = ½ρV∞², V∞ = %.4g.", model.freeStreamSpeed);
        }
    }
    ImGui::End();
    return picked;
}

} // namespace pfl
