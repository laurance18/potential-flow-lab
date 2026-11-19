#define IMGUI_DEFINE_MATH_OPERATORS
#include "ui/toolbar.hpp"

#include "ui/style.hpp"

namespace pfl {
namespace {

constexpr float k_buttonSize = 44.0f;
constexpr float k_padding = 6.0f;

} // namespace

float toolbarWidth() { return (k_buttonSize + 2.0f * k_padding) * uiScale(); }

ToolbarAction drawToolbar(ImVec2 pos, ImVec2 size, const ToolbarState& state) {
    ToolbarAction action;
    const float s = uiScale();
    const float button = k_buttonSize * s;

    ImGui::SetNextWindowPos(pos);
    ImGui::SetNextWindowSize(size);
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(k_padding * s, k_padding * s));
    ImGui::PushStyleVar(ImGuiStyleVar_ItemSpacing, ImVec2(0.0f, k_padding * s));
    constexpr ImGuiWindowFlags flags = ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_NoMove |
                                       ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoDocking |
                                       ImGuiWindowFlags_NoBringToFrontOnFocus;
    ImGui::Begin("##toolbar", nullptr, flags);

    ImDrawList* dl = ImGui::GetWindowDrawList();
    // The button just submitted, with its background; returns its centre for the icon.
    const auto buttonBackground = [&](bool active) {
        const ImVec2 min = ImGui::GetItemRectMin();
        const ImVec2 max = ImGui::GetItemRectMax();
        const ImGuiCol background = ImGui::IsItemActive() || active ? ImGuiCol_ButtonActive
                                    : ImGui::IsItemHovered()        ? ImGuiCol_ButtonHovered
                                                                    : ImGuiCol_Button;
        dl->AddRectFilled(min, max, ImGui::GetColorU32(background), 6.0f * s);
        return (min + max) * 0.5f;
    };
    if (ImGui::InvisibleButton("##click-to-place", ImVec2(button, button))) action.toggleClickToPlace = true;
    drawPlaceToolIcon(dl, buttonBackground(state.clickToPlace), 0.78f * button, state.clickToPlace);
    ImGui::SetItemTooltip("Click to place: %s\n"
                          "On: after choosing an element below, click on the plot where it goes; its other\n"
                          "parameters are set in the dialog that opens (the position can still be typed there).\n"
                          "Shift snaps to a grid, Esc cancels. Off: the dialog opens at once, at the view's centre.\n"
                          "A uniform flow has no position: its dialog always opens at once.",
                          state.clickToPlace ? "on" : "off");
    ImGui::Separator();

    for (const ElementKind kind : k_elementKinds) {
        if (isPanelKind(kind) && !isPanelKind(k_elementKinds[static_cast<std::size_t>(kind) - 1])) ImGui::Separator(); // panel methods
        ImGui::PushID(static_cast<int>(kind));
        if (ImGui::InvisibleButton("##tool", ImVec2(button, button))) action.element = kind;
        drawToolIcon(dl, buttonBackground(state.placing == kind), 0.78f * button, kind);
        if (kind == ElementKind::ThinAirfoil) {
            ImGui::SetItemTooltip("%s\nThin-airfoil theory: a vortex sheet on the camber line, with the Kutta condition "
                                  "at the trailing edge.", kindName(kind));
        } else {
            ImGui::SetItemTooltip("%s", kindName(kind));
        }
        ImGui::PopID();
        if (kind == ElementKind::PanelBody) {
            if (ImGui::InvisibleButton("##draw-body", ImVec2(button, button))) action.drawBody = true;
            drawOutlineToolIcon(dl, buttonBackground(state.drawingBody), 0.78f * button);
            ImGui::SetItemTooltip("Draw a source-panel body\nClick points around the shape (or drag to sketch it), then "
                                  "click the first point to close it.");
        }
    }

    ImGui::Separator();
    if (ImGui::InvisibleButton("##rake", ImVec2(button, button))) action.placeRake = true;
    drawDyeRakeIcon(dl, buttonBackground(state.placingRake), 0.78f * button);
    ImGui::SetItemTooltip("Dye rake\nDrag across the plot to place a line of dye emitters.");

    ImGui::End();
    ImGui::PopStyleVar(2);
    return action;
}

} // namespace pfl
