#include "ui/record_dialog.hpp"

#include "ui/style.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <span>

namespace pfl {
namespace {

constexpr const char* k_popupId = "Record animation";
constexpr std::array k_mp4Rates = {24, 25, 30, 50, 60};
constexpr std::array k_gifRates = {10, 15, 20, 25, 30, 50}; // GIF delays are whole 1/100 s, at least 2
constexpr std::array k_scales = {100, 75, 50};
constexpr double k_minSeconds = 0.5;
constexpr double k_maxSeconds = 600.0;

std::span<const int> frameRates(VideoFormat format) {
    if (format == VideoFormat::Gif) return k_gifRates;
    return k_mp4Rates;
}

int defaultFrameRate(VideoFormat format) { return format == VideoFormat::Gif ? 20 : 30; }

} // namespace

VideoSize videoSize(int width, int height, const RecordSettings& settings) {
    const auto scaled = [&](int pixels) {
        const int n = std::max(static_cast<int>(std::lround(pixels * settings.scalePercent / 100.0)), 2);
        return settings.format == VideoFormat::Mp4 ? n & ~1 : n;
    };
    return {scaled(width), scaled(height)};
}

int frameCount(const RecordSettings& settings) {
    return std::max(static_cast<int>(std::lround(settings.seconds * settings.fps)), 1);
}

void RecordDialog::open(std::string sweep) {
    m_sweep = std::move(sweep);
    if (m_sweep.empty()) m_settings.sweep = false;
    m_openRequested = true;
}

std::optional<RecordSettings> RecordDialog::draw(ImVec2 plotPixels) {
    if (m_openRequested) {
        ImGui::OpenPopup(k_popupId);
        m_openRequested = false;
    }
    ImGui::SetNextWindowPos(ImGui::GetMainViewport()->GetCenter(), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    if (!ImGui::BeginPopupModal(k_popupId, nullptr, ImGuiWindowFlags_AlwaysAutoResize)) return std::nullopt;

    const float s = uiScale();
    const float fieldWidth = 170.0f * s;
    RecordSettings& r = m_settings;

    int format = static_cast<int>(r.format);
    bool formatChanged = ImGui::RadioButton("MP4 video (H.264)", &format, static_cast<int>(VideoFormat::Mp4));
    ImGui::SetItemTooltip("Small files with smooth colours; plays everywhere (PowerPoint, browsers, phones).");
    ImGui::SameLine();
    formatChanged |= ImGui::RadioButton("Animated GIF", &format, static_cast<int>(VideoFormat::Gif));
    ImGui::SetItemTooltip("Loops forever and plays wherever an image can go (web pages, chat, Markdown).\n"
                          "Each frame has at most 255 colours; larger files than MP4.");
    if (formatChanged) {
        r.format = static_cast<VideoFormat>(format);
        const std::span<const int> rates = frameRates(r.format);
        if (std::ranges::find(rates, r.fps) == rates.end()) r.fps = defaultFrameRate(r.format);
    }

    float seconds = static_cast<float>(r.seconds);
    ImGui::SetNextItemWidth(fieldWidth);
    if (ImGui::SliderFloat("Length", &seconds, 1.0f, 30.0f, "%.1f s")) {
        r.seconds = std::clamp(static_cast<double>(seconds), k_minSeconds, k_maxSeconds);
    }
    ImGui::SetItemTooltip("Seconds of video. Ctrl+click to type a longer one (up to 10 minutes).");

    ImGui::SetNextItemWidth(fieldWidth);
    if (ImGui::BeginCombo("Frame rate", std::format("{} fps", r.fps).c_str())) {
        for (const int fps : frameRates(r.format)) {
            if (ImGui::Selectable(std::format("{} fps", fps).c_str(), fps == r.fps)) r.fps = fps;
        }
        ImGui::EndCombo();
    }

    ImGui::SetNextItemWidth(fieldWidth);
    if (ImGui::BeginCombo("Size", std::format("{}%", r.scalePercent).c_str())) {
        for (const int scale : k_scales) {
            if (ImGui::Selectable(std::format("{}%", scale).c_str(), scale == r.scalePercent)) r.scalePercent = scale;
        }
        ImGui::EndCombo();
    }
    ImGui::SetItemTooltip("Of the plot as shown on screen. Smaller GIFs load faster.");
    const int width = static_cast<int>(std::lround(plotPixels.x));
    const int height = static_cast<int>(std::lround(plotPixels.y));
    const VideoSize size = videoSize(width, height, r);
    ImGui::TextDisabled("%d × %d pixels, %d frames", size.width, size.height, frameCount(r));

    ImGui::Spacing();
    ImGui::Checkbox("Start from fresh particles and dye", &r.restart);
    ImGui::SetItemTooltip("Respawn the tracer particles and clear the dye streaks when the recording starts,\n"
                          "so that the video shows the dye leaving the rakes.");
    ImGui::BeginDisabled(m_sweep.empty());
    const std::string sweepLabel = m_sweep.empty() ? std::string("Play the parameter sweep") : "Sweep " + m_sweep;
    ImGui::Checkbox(sweepLabel.c_str(), &r.sweep);
    ImGui::EndDisabled();
    ImGui::SetItemTooltip(m_sweep.empty() ? "Set up a sweep in the Parameters panel (View → Parameters) first."
                                          : "Step the parameter once from From to To over the length of the video.");

    ImGui::Spacing();
    ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + 380.0f * s);
    ImGui::TextDisabled("The plot is recorded as drawn, axes included. You can keep working while it records: "
                        "drag elements or move sliders and the video shows it. The animation runs at the "
                        "video's frame rate, so the video plays at the right speed even if the screen lags.");
    ImGui::PopTextWrapPos();

    ImGui::Spacing();
    std::optional<RecordSettings> result;
    const ImVec2 buttonSize(120.0f * s, 0.0f);
    if (ImGui::Button("Record...", buttonSize)) {
        result = r;
        ImGui::CloseCurrentPopup();
    }
    ImGui::SameLine();
    if (ImGui::Button("Cancel", buttonSize) || ImGui::IsKeyPressed(ImGuiKey_Escape)) ImGui::CloseCurrentPopup();
    ImGui::EndPopup();
    return result;
}

} // namespace pfl
