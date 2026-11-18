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

} // namespace

float surfacePanelWidth(float available) {
    const float s = uiScale();
    return std::clamp(0.34f * available, 320.0f * s, 560.0f * s);
}

std::optional<Vec2> SurfacePanel::draw(ImVec2 pos, ImVec2 size, const std::vector<BodySurface>& bodies,
                                       const FlowModel& model) {
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
    } else if (bodies.empty()) {
        ImGui::TextWrapped("No body found. Bodies are bounded by the dividing streamlines through stagnation "
                           "points: add a doublet (cylinder), a source (half-body) or a source and a sink "
                           "(Rankine oval) to the uniform flow.\n\nA body is only found while its dividing "
                           "streamlines close (within 10%%). Asymmetric flows can leave them open: a source and "
                           "sink at incidence, for example, let part of the source's fluid escape downstream.");
    } else {
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
