#include "ui/status_bar.hpp"

#include "ui/style.hpp"

#include <algorithm>
#include <format>

namespace pfl {
namespace {

constexpr ImVec4 k_warningColor{1.0f, 0.8f, 0.35f, 1.0f};

void separator() {
    ImGui::SameLine();
    ImGui::TextDisabled("|");
    ImGui::SameLine();
}

} // namespace

float statusBarHeight() { return ImGui::GetFrameHeight() + 4.0f * uiScale(); }

void drawStatusBar(ImVec2 pos, ImVec2 size, const StatusInfo& info) {
    ImGui::SetNextWindowPos(pos);
    ImGui::SetNextWindowSize(size);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(8.0f * uiScale(), 2.0f * uiScale()));
    constexpr ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
                                       ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoDocking |
                                       ImGuiWindowFlags_NoBringToFrontOnFocus;
    ImGui::Begin("##status", nullptr, flags);
    ImGui::PopStyleVar();
    ImGui::AlignTextToFramePadding();

    if (info.probe) {
        const ProbeInfo& p = *info.probe;
        std::string text = std::format("x = {:.4f}   y = {:.4f}", p.position.x, p.position.y);
        if (p.insideBody) {
            text += std::format("   (inside source-panel body #{})", *p.insideBody);
        } else if (!p.inDomain) {
            text += "   (outside the flow domain)";
        } else {
            // + 0.0 turns -0 into 0 (e.g. v = -Im W of an empty field).
            text += std::format("   u = {:.4g}   v = {:.4g}   |V| = {:.4g}", p.velocity.x + 0.0, p.velocity.y + 0.0,
                                length(p.velocity));
            if (p.pressureCoefficient) text += std::format("   Cp = {:.4g}", *p.pressureCoefficient);
        }
        ImGui::TextUnformatted(text.c_str());
    } else {
        ImGui::TextDisabled("Move the mouse over the plot to probe the velocity.");
    }

    separator();
    ImGui::Text("%zu streamlines", info.streamlineCount);

    if (!info.walls.empty()) {
        separator();
        const auto worst = std::ranges::max_element(info.walls, {}, &WallStatus::relativeNormalSpeed);
        ImGui::Text("Walls: max |v·n| / |V| = %.1e", worst->relativeNormalSpeed);
        if (ImGui::BeginItemTooltip()) {
            ImGui::TextUnformatted("Normal velocity sampled along each symmetry wall in view:");
            for (const WallStatus& w : info.walls) {
                ImGui::Text("Wall #%zu: max |v·n| / max |V| = %.2e (%d points)", w.number, w.relativeNormalSpeed,
                            w.samples);
            }
            ImGui::EndTooltip();
        }
    }

    if (!info.panelBodies.empty()) {
        separator();
        const auto worst = std::ranges::max_element(info.panelBodies, {}, &PanelStatus::relativeNormalSpeed);
        ImGui::Text("Panels: max |v·n| / |V| = %.1e", worst->relativeNormalSpeed);
        if (ImGui::BeginItemTooltip()) {
            ImGui::TextUnformatted("Flow through each source-panel body and vortex-sheet airfoil at its control points; the "
                                   "net source\nstrength Σλs of a body (0 for the exact solution of a closed body; it "
                                   "shrinks as N grows), the circulation of an airfoil:");
            for (const PanelStatus& b : info.panelBodies) {
                if (b.airfoil) {
                    ImGui::Text("Airfoil #%zu (N = %d): max |v·n| / max |V| = %.2e,  Γ = %.5g,  γ(TE) = 0 (Kutta)", b.number,
                                b.panels, b.relativeNormalSpeed, b.circulation);
                } else {
                    ImGui::Text("Body #%zu (N = %d): max |v·n| / max |V| = %.2e,  Σλs / (V∞ · perimeter) = %.2e", b.number,
                                b.panels, b.relativeNormalSpeed, b.relativeNetSource);
                }
            }
            ImGui::EndTooltip();
        }
    }

    if (!info.warnings.empty()) {
        separator();
        ImGui::TextColored(k_warningColor, "(!) %zu warning%s", info.warnings.size(),
                           info.warnings.size() == 1 ? "" : "s");
        if (ImGui::BeginItemTooltip()) {
            ImGui::PushTextWrapPos(ImGui::GetFontSize() * 36.0f);
            for (const std::string& w : info.warnings) ImGui::BulletText("%s", w.c_str());
            ImGui::PopTextWrapPos();
            ImGui::EndTooltip();
        }
    }
    ImGui::End();
}

} // namespace pfl
