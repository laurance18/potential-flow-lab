#include "ui/parameter_dialog.hpp"

#include "ui/style.hpp"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <numbers>
#include <string>
#include <vector>

namespace pfl {
namespace {

constexpr const char* k_popupId = "###element-parameters";
constexpr double k_degrees = 180.0 / std::numbers::pi;
constexpr ImVec4 k_errorColor{1.0f, 0.45f, 0.4f, 1.0f};
constexpr ImVec4 k_warningColor{1.0f, 0.8f, 0.35f, 1.0f};

bool inputValue(const char* label, double& value, const char* help) {
    ImGui::SetNextItemWidth(170.0f * uiScale());
    const bool changed = ImGui::InputDouble(label, &value, 0.0, 0.0, "%.6g");
    ImGui::SetItemTooltip("%s", help);
    return changed;
}

void inputAngle(const char* label, double& radians, const char* help) {
    double degrees = radians * k_degrees;
    if (inputValue(label, degrees, help)) radians = degrees / k_degrees;
}

void inputPoint(Vec2& p, const char* what) {
    inputValue("x₀", p.x, std::format("x coordinate of the {}", what).c_str());
    inputValue("y₀", p.y, std::format("y coordinate of the {}", what).c_str());
}

void note(const std::string& text) {
    ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + 380.0f * uiScale());
    ImGui::TextDisabled("%s", text.c_str());
    ImGui::PopTextWrapPos();
}

// Sum of the scene's uniform flows (for the doublet's cylinder-radius hint).
double freeStreamSpeed(const Scene& scene) {
    Complex w{};
    for (const SceneElement& e : scene.elements) {
        if (const auto* u = std::get_if<UniformFlow>(&e)) w += complexVelocity(*u, {});
    }
    return std::abs(w);
}

void drawPanelBodyFields(PanelBody& body) {
    const float width = 170.0f * uiScale();
    int shape = static_cast<int>(body.shape);
    ImGui::SetNextItemWidth(width);
    if (ImGui::Combo("Shape", &shape, "Ellipse (circle when a = b)\0Rectangle\0NACA 4-digit airfoil\0Outline (drawn or imported)\0")) {
        const PanelShape before = body.shape;
        body.shape = static_cast<PanelShape>(shape);
        if (body.shape == PanelShape::Naca && before != PanelShape::Naca) body.semiAxisB = 0.12 * body.semiAxisA;
        const bool even = body.shape == PanelShape::Rectangle || body.shape == PanelShape::Naca;
        if (even && body.panels % 2 != 0) ++body.panels;
    }
    const bool even = body.shape == PanelShape::Rectangle || body.shape == PanelShape::Naca;

    switch (body.shape) {
    case PanelShape::Ellipse:
        inputValue("Semi-axis a", body.semiAxisA, "Half the length along the body's own x axis (a circle's radius when a = b)");
        inputValue("Semi-axis b", body.semiAxisB, "Half the length across the body's own x axis");
        break;
    case PanelShape::Rectangle:
        inputValue("Half-width a", body.semiAxisA, "Half the length along the body's own x axis");
        inputValue("Half-height b", body.semiAxisB, "Half the length across it");
        break;
    case PanelShape::Naca: {
        ImGui::SetNextItemWidth(width);
        ImGui::InputInt("Camber m [%]", &body.nacaCamber);
        ImGui::SetItemTooltip("Maximum camber, % of the chord: the first digit of NACA m p tt (0 for a symmetric airfoil)");
        ImGui::BeginDisabled(body.nacaCamber == 0);
        ImGui::SetNextItemWidth(width);
        ImGui::InputInt("Its position p [c/10]", &body.nacaPosition);
        ImGui::SetItemTooltip("Where the camber is largest, in tenths of the chord from the leading edge: the second digit");
        ImGui::EndDisabled();
        double chord = 2.0 * body.semiAxisA;
        double thickness = 100.0 * body.semiAxisB / body.semiAxisA;
        if (inputValue("Chord c", chord, "Length from the leading to the trailing edge (c = 2a); the thickness ratio stays")) {
            body.semiAxisA = 0.5 * chord;
            body.semiAxisB = thickness / 100.0 * body.semiAxisA;
        }
        if (inputValue("Thickness t/c [%]", thickness, "Largest thickness, % of the chord: the last two digits (b = (t/c) a)")) {
            body.semiAxisB = thickness / 100.0 * body.semiAxisA;
        }
        break;
    }
    case PanelShape::Custom: {
        ImGui::Text("Outline: %s, %zu points", body.name.empty() ? "custom" : body.name.c_str(), body.outline.size());
        ImGui::Checkbox("Smooth outline", &body.smooth);
        ImGui::SetItemTooltip("On: a smooth curve through the points, keeping only sharp corners (over 100°, such as a "
                              "trailing edge).\nOff: straight lines from point to point, a corner at each.");
        inputValue("Half-width a", body.semiAxisA, "The outline is stretched to fit 2a along the body's own x axis");
        inputValue("Half-height b", body.semiAxisB, "... and 2b across it");
        if (body.outline.size() >= 3 && ImGui::Button("Keep the outline's proportions")) {
            double xLo = 1e300, xHi = -1e300, yLo = 1e300, yHi = -1e300;
            for (const Vec2 p : body.outline) {
                xLo = std::min(xLo, p.x);
                xHi = std::max(xHi, p.x);
                yLo = std::min(yLo, p.y);
                yHi = std::max(yHi, p.y);
            }
            if (xHi > xLo) body.semiAxisB = body.semiAxisA * (yHi - yLo) / (xHi - xLo);
        }
        ImGui::SetItemTooltip("Sets b so the outline is not stretched more one way than the other.");
        break;
    }
    }
    inputPoint(body.position, body.shape == PanelShape::Naca ? "middle of the chord" : "body's centre");
    inputAngle("Rotation θ [°]", body.angle, "Turns the body counter-clockwise about its centre");
    ImGui::SetNextItemWidth(width);
    ImGui::InputInt("Panels N", &body.panels, even ? 2 : 1, 10);
    ImGui::SetItemTooltip("Number of straight source panels around the surface (%d to %d%s)", k_minPanels, k_maxPanels,
                          even ? ", even" : "");

    std::string text = "Source panel method: the surface is cut into N straight panels, each a source sheet of its own "
                       "constant strength λ. The strengths are solved so that no fluid crosses any panel at its midpoint "
                       "(the control point), in the flow of all the other elements. The body is non-lifting: add a vortex "
                       "inside it to give it circulation. Inside the body the flow is not shown. More panels: closer to "
                       "the exact solution.";
    if (body.shape == PanelShape::Naca) {
        text = std::format("This is the {}. With a stream along +x the angle of attack is α = −θ. ", nacaName(body)) +
               "Panels are cosine-spaced along the chord, short at both edges. Source panels give an airfoil no lift "
               "(no circulation, no Kutta condition): that needs a vortex sheet.";
    } else if (body.shape == PanelShape::Custom) {
        text = "The points only give the shape: the N panels are laid along it, a vertex at every corner, shorter where "
               "it curves. Draw an outline with the pencil tool below the element icons, or import airfoil coordinates "
               "(File → Import airfoil).";
    }
    note(text);
}

void drawFields(SceneElement& element, const Scene& scene) {
    if (auto* uniform = std::get_if<UniformFlow>(&element)) {
        inputValue("Speed U∞", uniform->speed, "Free-stream speed (not negative)");
        inputAngle("Direction α [°]", uniform->angle, "Measured counter-clockwise from +x");
        note(std::format("u = U∞ cos α = {:.4g},  v = U∞ sin α = {:.4g}", uniform->speed * std::cos(uniform->angle),
                         uniform->speed * std::sin(uniform->angle)));
    } else if (auto* source = std::get_if<SourceSink>(&element)) {
        inputValue("Strength m = Λ", source->strength, "Volume flow rate per unit depth: > 0 source, < 0 sink");
        inputPoint(source->position, "source / sink");
        note("Λ > 0 is a source, Λ < 0 a sink. Radial velocity V_r = Λ / (2π r).");
    } else if (auto* vortex = std::get_if<Vortex>(&element)) {
        inputValue("Circulation Γ", vortex->circulation, "> 0 clockwise, < 0 counter-clockwise");
        inputPoint(vortex->position, "vortex");
        note("Γ > 0 turns clockwise, Γ < 0 counter-clockwise (Anderson's convention). V_θ = −Γ / (2π r).");
    } else if (auto* doublet = std::get_if<Doublet>(&element)) {
        inputValue("Strength μ", doublet->strength, "μ = κ / 2π, as in F = μ e^(iβ) / (z − z₀)");
        inputPoint(doublet->position, "doublet");
        inputAngle("Axis angle β [°]", doublet->angle, "Rotates the doublet axis counter-clockwise");
        std::string text = "F = μ e^(iβ) / (z − z₀). With a uniform flow at angle α, β = α makes a circular cylinder.";
        if (const double u = freeStreamSpeed(scene); u > 0.0 && doublet->strength > 0.0) {
            text += std::format(" With U∞ = {:.4g}: radius a = √(μ/U∞) = {:.4g}.", u, std::sqrt(doublet->strength / u));
        }
        note(text);
    } else if (auto* sector = std::get_if<SectorFlow>(&element)) {
        inputValue("Strength V", sector->strength, "Scaling constant in F = V (ζ e^(−iθ₁))ⁿ");
        inputPoint(sector->vertex, "wedge vertex");
        inputAngle("First wall θ₁ [°]", sector->wall1, "Angle of the first wall, counter-clockwise from +x");
        inputAngle("Second wall θ₂ [°]", sector->wall2, "Angle of the second wall; 0 < θ₂ − θ₁ ≤ 360°");
        if (sector->opening() > 0.0) {
            note(std::format("Opening β = θ₂ − θ₁ = {:.4g}°, n = π/β = {:.4g}. The flow exists only inside the "
                             "wedge; the plot greys out the rest. Other elements are added without "
                             "restriction, so they generally cross the wedge walls.",
                             sector->opening() * k_degrees, sector->exponent()));
        }
    } else if (auto* body = std::get_if<PanelBody>(&element)) {
        drawPanelBodyFields(*body);
    } else if (auto* wall = std::get_if<SymmetryWall>(&element)) {
        inputPoint(wall->point, "point on the wall");
        inputAngle("Angle φ [°]", wall->angle, "Direction of the wall line, counter-clockwise from +x");
        note("An infinite straight wall. Zero normal velocity is imposed with images: sources, vortices and "
             "doublets are mirrored across it (vortices with −Γ), and uniform flow keeps only its component "
             "along the wall. Sector flows are not mirrored.");
    }
}

} // namespace

Scene applyEdit(Scene scene, const ElementEdit& edit) {
    if (edit.index && *edit.index < scene.elements.size()) {
        scene.elements[*edit.index] = edit.element;
    } else {
        scene.elements.push_back(edit.element);
    }
    return scene;
}

void ParameterDialog::openNew(SceneElement element) { open(std::nullopt, std::move(element)); }

void ParameterDialog::openEdit(std::size_t index, SceneElement element) { open(index, std::move(element)); }

void ParameterDialog::open(std::optional<std::size_t> index, SceneElement element) {
    m_index = index;
    m_pending = std::move(element);
    m_openRequested = true;
}

std::optional<ElementEdit> ParameterDialog::draw(const Scene& scene) {
    if (m_openRequested) {
        ImGui::OpenPopup(k_popupId);
        m_openRequested = false;
    }

    const char* kind = kindName(kindOf(m_pending));
    const std::string title = (m_index ? std::format("Edit #{}  {}", *m_index + 1, kind) : std::format("Add {}", kind)) +
                              k_popupId;
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    if (!ImGui::BeginPopupModal(title.c_str(), nullptr, ImGuiWindowFlags_AlwaysAutoResize)) return std::nullopt;

    std::optional<ElementEdit> result;
    drawFields(m_pending, scene);

    const std::optional<std::string> error = sceneElementError(m_pending);
    const ElementEdit edit{m_index, m_pending};
    if (error) {
        ImGui::Spacing();
        ImGui::PushStyleColor(ImGuiCol_Text, k_errorColor);
        ImGui::TextWrapped("%s", error->c_str());
        ImGui::PopStyleColor();
    } else {
        // Warnings this change would introduce.
        // The warnings do not depend on the panel solution, so the panels are not solved here.
        const std::vector<std::string> before = sceneWarnings(scene, buildFlowModel(scene, false));
        const Scene after = applyEdit(scene, edit);
        ImGui::PushStyleColor(ImGuiCol_Text, k_warningColor);
        ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + 380.0f * uiScale());
        for (const std::string& warning : sceneWarnings(after, buildFlowModel(after, false))) {
            if (std::ranges::find(before, warning) == before.end()) ImGui::TextWrapped("%s", warning.c_str());
        }
        ImGui::PopTextWrapPos();
        ImGui::PopStyleColor();
    }

    ImGui::Spacing();
    const ImVec2 buttonSize(120.0f * uiScale(), 0.0f);
    ImGui::BeginDisabled(error.has_value());
    if (ImGui::Button("Apply", buttonSize)) {
        result = edit;
        ImGui::CloseCurrentPopup();
    }
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (ImGui::Button("Cancel", buttonSize) || ImGui::IsKeyPressed(ImGuiKey_Escape)) ImGui::CloseCurrentPopup();
    ImGui::EndPopup();
    return result;
}

} // namespace pfl
