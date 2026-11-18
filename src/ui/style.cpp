#define IMGUI_DEFINE_MATH_OPERATORS
#include "ui/style.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace pfl {
namespace {

constexpr float k_pi = std::numbers::pi_v<float>;

ImVec2 polar(float r, float angle) { return {r * std::cos(angle), r * std::sin(angle)}; }

} // namespace

ImU32 fadeColor(ImU32 colour, float alpha) {
    const auto a = static_cast<ImU32>(static_cast<float>((colour >> IM_COL32_A_SHIFT) & 0xFF) * alpha);
    return (colour & ~IM_COL32_A_MASK) | (a << IM_COL32_A_SHIFT);
}

ImU32 pressureColor(double cp, double lowest, float alpha) {
    constexpr ImVec4 neutral{56.0f / 255.0f, 56.0f / 255.0f, 53.0f / 255.0f, 1.0f};
    constexpr ImVec4 pressure{230.0f / 255.0f, 103.0f / 255.0f, 103.0f / 255.0f, 1.0f};
    constexpr ImVec4 suction{85.0f / 255.0f, 152.0f / 255.0f, 231.0f / 255.0f, 1.0f};
    const bool positive = cp >= 0.0;
    const float t = static_cast<float>(std::clamp(positive ? cp : cp / lowest, 0.0, 1.0));
    const ImVec4& pole = positive ? pressure : suction;
    return ImGui::ColorConvertFloat4ToU32(
        {neutral.x + (pole.x - neutral.x) * t, neutral.y + (pole.y - neutral.y) * t, neutral.z + (pole.z - neutral.z) * t, alpha});
}

void drawArrowHead(ImDrawList* dl, ImVec2 tip, ImVec2 dir, float size, ImU32 colour) {
    const ImVec2 normal{-dir.y, dir.x};
    const ImVec2 base = tip - dir * size;
    dl->AddTriangleFilled(tip, base + normal * (0.45f * size), base - normal * (0.45f * size), colour);
}

void drawVortexGlyph(ImDrawList* dl, ImVec2 centre, float radius, bool clockwise, ImU32 colour) {
    const float s = uiScale();
    // Screen angles grow clockwise because screen y points down.
    const float a0 = 0.6f;
    const float a1 = a0 + 1.6f * k_pi;
    dl->PathArcTo(centre, radius, a0, a1, 24);
    dl->PathStroke(colour, ImDrawFlags_None, 2.0f * s);
    const float end = clockwise ? a1 : a0;
    const ImVec2 tangent = polar(1.0f, end + 0.5f * k_pi) * (clockwise ? 1.0f : -1.0f);
    drawArrowHead(dl, centre + polar(radius, end) + tangent * (3.0f * s), tangent, 6.0f * s, colour);
    dl->AddCircleFilled(centre, 2.0f * s, colour);
}

void drawDoubletGlyph(ImDrawList* dl, ImVec2 centre, float radius, ImVec2 axis, ImU32 colour) {
    const float s = uiScale();
    const ImVec2 normal{-axis.y, axis.x};
    dl->AddCircle(centre + normal * radius, radius, colour, 0, 1.5f * s);
    dl->AddCircle(centre - normal * radius, radius, colour, 0, 1.5f * s);
    dl->AddLine(centre + axis * (1.4f * radius), centre - axis * (1.4f * radius), colour, 1.8f * s);
    drawArrowHead(dl, centre - axis * (2.2f * radius), -axis, 0.9f * radius + 2.0f * s, colour);
}

void drawToolIcon(ImDrawList* dl, ImVec2 c, float size, ElementKind kind) {
    const float s = uiScale();
    const float h = 0.5f * size;
    switch (kind) {
    case ElementKind::Uniform:
        for (const float dy : {-0.45f, 0.0f, 0.45f}) {
            const ImVec2 a = c + ImVec2(-0.7f * h, dy * h);
            const ImVec2 b = c + ImVec2(0.55f * h, dy * h);
            dl->AddLine(a, b, k_uniformColor, 1.8f * s);
            drawArrowHead(dl, b + ImVec2(0.25f * h, 0.0f), ImVec2(1.0f, 0.0f), 0.3f * h, k_uniformColor);
        }
        break;
    case ElementKind::SourceSink:
        for (int k = 0; k < 8; ++k) {
            const float a = k * 0.25f * k_pi;
            const ImVec2 dir = polar(1.0f, a);
            dl->AddLine(c + dir * (0.3f * h), c + dir * (0.7f * h), k_sourceColor, 1.5f * s);
            drawArrowHead(dl, c + dir * (0.9f * h), dir, 0.22f * h, k_sourceColor);
        }
        dl->AddCircleFilled(c, 0.2f * h, k_sourceColor);
        break;
    case ElementKind::Vortex: drawVortexGlyph(dl, c, 0.62f * h, true, k_vortexColor); break;
    case ElementKind::Doublet: drawDoubletGlyph(dl, c, 0.36f * h, ImVec2(1.0f, 0.0f), k_doubletColor); break;
    case ElementKind::Sector: {
        // 90-degree corner with two hyperbolic streamlines x y = const.
        const ImVec2 o = c + ImVec2(-0.75f * h, 0.75f * h);
        dl->AddLine(o, o + ImVec2(1.6f * h, 0.0f), k_sectorColor, 2.2f * s);
        dl->AddLine(o, o - ImVec2(0.0f, 1.6f * h), k_sectorColor, 2.2f * s);
        for (const float k : {0.18f, 0.5f}) {
            for (int i = 0; i < 12; ++i) {
                const float x = (k / 1.5f) + (1.5f - k / 1.5f) * i / 11.0f;
                dl->PathLineTo(o + ImVec2(x * h, -(k / x) * h));
            }
            dl->PathStroke(k_sectorColor, ImDrawFlags_None, 1.4f * s);
        }
        break;
    }
    case ElementKind::Wall: {
        const ImVec2 a = c + ImVec2(-0.8f * h, 0.1f * h);
        const ImVec2 b = c + ImVec2(0.8f * h, 0.1f * h);
        dl->AddLine(a, b, k_wallColor, 2.5f * s);
        for (int i = 0; i < 6; ++i) {
            const ImVec2 p = a + (b - a) * ((i + 0.5f) / 6.0f);
            dl->AddLine(p, p + ImVec2(-0.25f * h, 0.35f * h), fadeColor(k_wallColor, 0.7f), 1.2f * s);
        }
        // A streamline above the wall, mirrored by a faint one below.
        dl->AddLine(c + ImVec2(-0.8f * h, -0.45f * h), c + ImVec2(0.8f * h, -0.45f * h), k_uniformColor, 1.2f * s);
        drawArrowHead(dl, c + ImVec2(0.85f * h, -0.45f * h), ImVec2(1.0f, 0.0f), 0.25f * h, k_uniformColor);
        break;
    }
    case ElementKind::PanelBody: {
        // An ellipse cut into straight panels, each with a small source arrow out of its midpoint.
        constexpr int panels = 10;
        std::array<ImVec2, panels> v;
        for (int k = 0; k < panels; ++k) v[static_cast<std::size_t>(k)] = c + ImVec2(0.62f * h * std::cos(k * 2.0f * k_pi / panels), 0.4f * h * std::sin(k * 2.0f * k_pi / panels));
        dl->AddConvexPolyFilled(v.data(), panels, k_panelFillColor);
        dl->AddPolyline(v.data(), panels, k_panelColor, ImDrawFlags_Closed, 1.6f * s);
        for (int k = 0; k < panels; k += 2) {
            const ImVec2 a = v[static_cast<std::size_t>(k)];
            const ImVec2 b = v[static_cast<std::size_t>((k + 1) % panels)];
            const ImVec2 m = (a + b) * 0.5f;
            ImVec2 n(b.y - a.y, a.x - b.x); // outward: the points run clockwise on screen (y down)
            n = n / std::sqrt(n.x * n.x + n.y * n.y);
            dl->AddLine(m, m + n * (0.22f * h), k_sourceColor, 1.3f * s);
            drawArrowHead(dl, m + n * (0.38f * h), n, 0.18f * h, k_sourceColor);
        }
        break;
    }
    }
}

void drawDyeRakeIcon(ImDrawList* dl, ImVec2 c, float size) {
    const float s = uiScale();
    const float h = 0.5f * size;
    const ImVec2 top = c + ImVec2(-0.7f * h, -0.75f * h);
    const ImVec2 bottom = c + ImVec2(-0.7f * h, 0.75f * h);
    dl->AddLine(top, bottom, IM_COL32(200, 200, 210, 255), 1.6f * s);
    for (int i = 0; i < 3; ++i) {
        const ImVec2 o = top + (bottom - top) * ((i + 0.5f) / 3.0f);
        const ImU32 colour = k_dyeColors[static_cast<std::size_t>(i)];
        for (int k = 0; k <= 12; ++k) {
            const float t = k / 12.0f;
            dl->PathLineTo(o + ImVec2(1.5f * h * t, 0.12f * h * std::sin(t * 2.0f * k_pi)));
        }
        dl->PathStroke(colour, ImDrawFlags_None, 2.0f * s);
        dl->AddCircleFilled(o, 2.2f * s, colour);
    }
}

} // namespace pfl
