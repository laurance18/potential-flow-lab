// style.hpp: colours and small vector glyphs shared by the toolbar icons and the plot markers.
// Glyphs are drawn with ImDrawList in screen pixels, so no image assets are needed.
#pragma once

#include "core/scene.hpp"

#include <imgui.h>

#include <array>
#include <cstddef>

namespace pfl {

inline constexpr ImU32 k_sourceColor = IM_COL32(95, 205, 120, 255);
inline constexpr ImU32 k_sinkColor = IM_COL32(235, 95, 85, 255);
inline constexpr ImU32 k_vortexColor = IM_COL32(245, 185, 65, 255);
inline constexpr ImU32 k_doubletColor = IM_COL32(195, 140, 245, 255);
inline constexpr ImU32 k_sectorColor = IM_COL32(90, 205, 225, 255);
inline constexpr ImU32 k_wallColor = IM_COL32(235, 235, 235, 255);
inline constexpr ImU32 k_uniformColor = IM_COL32(235, 235, 235, 255);
inline constexpr ImU32 k_panelColor = IM_COL32(150, 195, 255, 255);     // panel bodies: outline and markers
inline constexpr ImU32 k_panelFillColor = IM_COL32(44, 50, 64, 255);   // ... and the solid inside
inline constexpr ImU32 k_sheetColor = IM_COL32(255, 205, 140, 255);     // thin airfoils: the camber line

// Bodies in the surface-pressure panel, in this fixed order (identity, not rank).
inline constexpr std::array<ImU32, 8> k_bodyColors = {
    IM_COL32(57, 135, 229, 255), IM_COL32(217, 89, 38, 255),  IM_COL32(25, 158, 112, 255), IM_COL32(201, 133, 0, 255),
    IM_COL32(213, 81, 129, 255), IM_COL32(0, 131, 0, 255),    IM_COL32(144, 133, 233, 255), IM_COL32(230, 103, 103, 255)};

// Dye rakes, by rake index (cycling): bright hues that stand out from the blue streamlines.
inline constexpr std::array<ImU32, 6> k_dyeColors = {IM_COL32(255, 90, 170, 255), IM_COL32(110, 235, 120, 255),
                                                     IM_COL32(255, 175, 50, 255), IM_COL32(80, 215, 255, 255),
                                                     IM_COL32(205, 135, 255, 255), IM_COL32(250, 240, 90, 255)};
inline ImU32 dyeColor(std::size_t rake) { return k_dyeColors[rake % k_dyeColors.size()]; }

// Diverging colour for a pressure coefficient: a neutral grey at Cp = 0 (free-stream pressure),
// red towards Cp = 1 (stagnation) and blue towards `lowest` (< 0, the strongest suction shown).
// Values beyond either end take the end colour.
ImU32 pressureColor(double cp, double lowest, float alpha = 1.0f);

// DPI scale of the UI (1 at 96 dpi).
inline float uiScale() { return ImGui::GetStyle().FontScaleDpi; }

// `colour` with its alpha multiplied by `alpha`.
ImU32 fadeColor(ImU32 colour, float alpha);

// Filled arrow head with its tip at `tip`, pointing along the unit vector `dir`.
void drawArrowHead(ImDrawList* dl, ImVec2 tip, ImVec2 dir, float size, ImU32 colour);

// Circular arrow around `centre`; clockwise on screen when `clockwise` (a vortex with Gamma > 0).
void drawVortexGlyph(ImDrawList* dl, ImVec2 centre, float radius, bool clockwise, ImU32 colour);

// Two circles tangent to the doublet axis `axis` (unit screen vector, beta = 0 points right) and an
// arrow along -axis, the direction of Anderson's doublet arrow (sink to source).
void drawDoubletGlyph(ImDrawList* dl, ImVec2 centre, float radius, ImVec2 axis, ImU32 colour);

// Toolbar icon for drawing a source-panel body: a pencil finishing a closed outline.
void drawOutlineToolIcon(ImDrawList* dl, ImVec2 centre, float size);

// Toolbar icon for the dye rake tool: a rake with streaks trailing from it.
void drawDyeRakeIcon(ImDrawList* dl, ImVec2 centre, float size);

// Toolbar icon for an element type, fitted in a square of side `size` centred on `centre`.
void drawToolIcon(ImDrawList* dl, ImVec2 centre, float size, ElementKind kind);

} // namespace pfl
