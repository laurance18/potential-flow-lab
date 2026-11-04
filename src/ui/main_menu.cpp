#include "ui/main_menu.hpp"

#include "io/examples.hpp"
#include "ui/style.hpp"

#include <imgui.h>

#include <format>

namespace pfl {
namespace {

constexpr const char* k_discardPopup = "Unsaved changes";

void drawEditMenu(const MainMenuState& state, MainMenuActions& actions) {
    if (!ImGui::BeginMenu("Edit")) return;
    if (ImGui::MenuItem("Undo", "Ctrl+Z", false, state.canUndo)) actions.undo = true;
    ImGui::SetItemTooltip("Take back the last change to the elements or dye rakes (a whole drag, slider move or\n"
                          "sweep is one change). The view and the View settings are not affected.");
    if (ImGui::MenuItem("Redo", "Ctrl+Y", false, state.canRedo)) actions.redo = true;
    ImGui::EndMenu();
}

void drawExamplesMenu(MainMenuActions& actions) {
    if (!ImGui::BeginMenu("Examples")) return;
    const std::span<const ExampleScene> examples = exampleScenes();
    for (std::size_t i = 0; i < examples.size(); ++i) {
        if (ImGui::MenuItem(examples[i].title)) actions.example = i;
        ImGui::SetItemTooltip("%s", examples[i].description);
    }
    ImGui::EndMenu();
}

void drawViewMenu(PlotSettings& plot, bool parameterPanel, MainMenuActions& actions) {
    if (!ImGui::BeginMenu("View")) return;
    if (ImGui::MenuItem("Parameters", "Ctrl+P", parameterPanel)) actions.toggleParameterPanel = true;
    ImGui::SetItemTooltip("A panel with a live slider for every element parameter, and the parameter sweep\n"
                          "(animate one parameter, or export the stagnation points and forces along it).");
    if (ImGui::MenuItem("Reset view")) actions.resetView = true;
    ImGui::SetItemTooltip("Show the default region -4 < x < 4, -2.5 < y < 2.5.\n"
                          "Scroll to zoom and drag to pan the plot.");
    ImGui::Separator();

    ImGui::SetNextItemWidth(180.0f * uiScale());
    actions.plotChanged |= ImGui::SliderInt("Streamline density", &plot.streamlineDensity, 5, 80);
    ImGui::SetItemTooltip("Number of streamlines across the height of the view.");
    actions.plotChanged |= ImGui::Checkbox("Direction arrows", &plot.showArrows);
    ImGui::BeginDisabled(!plot.showArrows);
    float spacing = static_cast<float>(plot.arrowSpacingPx);
    ImGui::SetNextItemWidth(180.0f * uiScale());
    if (ImGui::SliderFloat("Arrow spacing", &spacing, 40.0f, 400.0f, "%.0f px")) {
        plot.arrowSpacingPx = spacing;
        actions.plotChanged = true;
    }
    ImGui::EndDisabled();
    actions.plotChanged |= ImGui::Checkbox("Show wall images", &plot.showImages);
    ImGui::SetItemTooltip("Draw the mirrored sources, vortices and doublets that enforce the symmetry walls.");
    actions.plotChanged |= ImGui::Checkbox("Equipotential lines", &plot.showEquipotentials);
    ImGui::SetItemTooltip("Lines of constant velocity potential φ, crossing every streamline at right angles.");
    int background = static_cast<int>(plot.background);
    ImGui::SetNextItemWidth(180.0f * uiScale());
    if (ImGui::Combo("Background", &background, "None\0Speed contours |V|\0Pressure coefficient Cp\0")) {
        plot.background = static_cast<Background>(background);
        actions.plotChanged = true;
    }
    ImGui::SetItemTooltip("Filled bands under the streamlines:\n"
                          "|V| from 0 to twice the typical speed in the view, or\n"
                          "Cp = 1 - |V|²/V∞² (Bernoulli, free stream as reference): red above and blue below\n"
                          "the free-stream pressure. Cp needs a uniform flow.");
    actions.plotChanged |= ImGui::Checkbox("Surface pressure", &plot.showSurfacePressure);
    ImGui::SetItemTooltip("Find the bodies (closed or half-bodies bounded by dividing streamlines), colour their\n"
                          "surfaces by Cp with pressure arrows, and open the panel with Cp along the surface\n"
                          "and the integrated lift and drag.");
    actions.plotChanged |= ImGui::Checkbox("Stagnation points", &plot.showStagnationPoints);
    ImGui::SetItemTooltip("Mark the points where the velocity is zero (Cp = 1) with diamonds.");
    actions.plotChanged |= ImGui::Checkbox("Tracer particles", &plot.showParticles);
    ImGui::BeginDisabled(!plot.showParticles);
    actions.plotChanged |= ImGui::Checkbox("Particle trails", &plot.particleTrails);
    ImGui::EndDisabled();
    ImGui::EndMenu();
}

void drawVerifyMenu(MainMenuActions& actions) {
    if (!ImGui::BeginMenu("Verify")) return;
    if (ImGui::MenuItem("Physics checks...", "F5")) actions.verify = true;
    ImGui::SetItemTooltip("Run the program against the exact results of potential-flow theory (elementary flows,\n"
                          "superposition, walls, surface pressure, lift and drag, streamlines) and see each one pass.");
    ImGui::EndMenu();
}

// Right end of the menu bar: the recording in progress, or the notice.
void drawMenuBarStatus(const MainMenuState& state, MainMenuActions& actions) {
    const ImGuiStyle& style = ImGui::GetStyle();
    if (state.recording) {
        const std::string text = std::format("Recording  {:.1f} / {:.1f} s", state.recording->seconds, state.recording->total);
        const float width = ImGui::GetFontSize() + ImGui::CalcTextSize(text.c_str()).x + ImGui::CalcTextSize("Stop").x +
                            2.0f * style.FramePadding.x + 3.0f * style.ItemSpacing.x;
        ImGui::SetCursorPosX(ImGui::GetWindowWidth() - width);
        // A red dot, like a camera's.
        const ImVec2 p = ImGui::GetCursorScreenPos();
        const float r = 0.3f * ImGui::GetFontSize();
        ImGui::GetWindowDrawList()->AddCircleFilled(ImVec2(p.x + r, p.y + 0.5f * ImGui::GetFrameHeight()), r,
                                                    IM_COL32(235, 60, 60, 255));
        ImGui::Dummy(ImVec2(2.0f * r, ImGui::GetFrameHeight()));
        ImGui::TextUnformatted(text.c_str());
        if (ImGui::SmallButton("Stop")) actions.stopRecording = true;
        ImGui::SetItemTooltip("Stop now and keep what has been recorded.");
    } else if (!state.notice.empty()) {
        const float width = ImGui::CalcTextSize(state.notice.c_str()).x + 2.0f * style.ItemSpacing.x;
        ImGui::SetCursorPosX(ImGui::GetWindowWidth() - width);
        ImGui::TextDisabled("%s", state.notice.c_str());
    }
}

} // namespace

MainMenuActions drawMainMenu(const Scene& scene, PlotSettings& plot, const MainMenuState& state) {
    MainMenuActions actions;
    const bool canReloadShaders = state.canReloadShaders;
    constexpr ImGuiInputFlags route = ImGuiInputFlags_RouteGlobal;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_N, route)) actions.newScene = true;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_O, route)) actions.open = true;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_S, route)) actions.save = true;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiMod_Shift | ImGuiKey_S, route)) actions.saveAs = true;
    if (canReloadShaders && ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_R, route)) actions.reloadShaders = true;
    if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_P, route)) actions.toggleParameterPanel = true;
    if (ImGui::Shortcut(ImGuiKey_F5, route)) actions.verify = true;
    // Not while a dialog or menu is open; a text field being edited takes Ctrl+Z for itself.
    if (!ImGui::IsPopupOpen("", ImGuiPopupFlags_AnyPopupId)) {
        if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_Z, route)) actions.undo = true;
        if (ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiKey_Y, route) ||
            ImGui::Shortcut(ImGuiMod_Ctrl | ImGuiMod_Shift | ImGuiKey_Z, route)) {
            actions.redo = true;
        }
    }

    if (ImGui::BeginMainMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("New", "Ctrl+N")) actions.newScene = true;
            ImGui::SetItemTooltip("Start an empty scene with default plot settings.");
            if (ImGui::MenuItem("Open...", "Ctrl+O")) actions.open = true;
            if (ImGui::MenuItem("Save", "Ctrl+S")) actions.save = true;
            if (ImGui::MenuItem("Save As...", "Ctrl+Shift+S")) actions.saveAs = true;
            ImGui::Separator();
            const bool hasFlow = !scene.elements.empty();
            if (ImGui::MenuItem("Export image (PNG)...", nullptr, false, hasFlow)) actions.exportImage = true;
            ImGui::SetItemTooltip("Save the plot, axes included, as it is shown now.");
            if (ImGui::MenuItem("Export field data (CSV)...", nullptr, false, hasFlow)) actions.exportField = true;
            ImGui::SetItemTooltip("x, y, u, v, |V|, φ, ψ and Cp on a 201-column grid over the current view.");
            if (ImGui::MenuItem("Export lines (CSV)...", nullptr, false, hasFlow)) actions.exportLines = true;
            ImGui::SetItemTooltip("The points of every streamline and equipotential line shown now.");
            if (ImGui::MenuItem("Export surface Cp (CSV)...", nullptr, false, hasFlow)) actions.exportSurface = true;
            ImGui::SetItemTooltip("x, y, arc length, x/c, angle and Cp (and the exact Cp for a cylinder) at every\n"
                                  "point of every body surface found in the flow.");
            if (state.recording) {
                if (ImGui::MenuItem("Stop recording")) actions.stopRecording = true;
            } else if (ImGui::MenuItem("Record animation (MP4 / GIF)...", nullptr, false, hasFlow)) {
                actions.record = true;
            }
            ImGui::SetItemTooltip("Record the animated plot (particles, dye, live edits, a parameter sweep) as an MP4\n"
                                  "video or an animated GIF.");
            ImGui::Separator();
            if (ImGui::MenuItem("Clear All", nullptr, false, !scene.elements.empty() || !scene.rakes.empty())) {
                actions.clearAll = true;
            }
            ImGui::SetItemTooltip("Remove every element and dye rake but keep the view and plot settings.");
            ImGui::Separator();
            if (canReloadShaders && ImGui::MenuItem("Reload shaders", "Ctrl+R")) actions.reloadShaders = true;
            if (ImGui::MenuItem("Quit", "Alt+F4")) actions.quit = true;
            ImGui::EndMenu();
        }
        drawEditMenu(state, actions);
        drawExamplesMenu(actions);
        actions.elements = drawElementsMenu(scene);
        drawViewMenu(plot, state.parameterPanel, actions);
        drawVerifyMenu(actions);
        drawMenuBarStatus(state, actions);
        ImGui::EndMainMenuBar();
    }
    return actions;
}

void openDiscardChangesDialog() { ImGui::OpenPopup(k_discardPopup); }

std::optional<bool> drawDiscardChangesDialog() {
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    if (!ImGui::BeginPopupModal(k_discardPopup, nullptr, ImGuiWindowFlags_AlwaysAutoResize)) return std::nullopt;

    std::optional<bool> choice;
    ImGui::TextUnformatted("The scene has unsaved changes. Discard them?");
    ImGui::Spacing();
    const ImVec2 buttonSize(120.0f * uiScale(), 0.0f);
    if (ImGui::Button("Discard", buttonSize)) choice = true;
    ImGui::SameLine();
    if (ImGui::Button("Cancel", buttonSize) || ImGui::IsKeyPressed(ImGuiKey_Escape)) choice = false;
    if (choice) ImGui::CloseCurrentPopup();
    ImGui::EndPopup();
    return choice;
}

} // namespace pfl
