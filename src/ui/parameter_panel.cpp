#include "ui/parameter_panel.hpp"

#include "ui/style.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <format>
#include <numbers>

namespace pfl {
namespace {

using std::numbers::pi;
constexpr double k_degrees = 180.0 / pi;
constexpr ImVec4 k_sweepColor{0.45f, 0.8f, 1.0f, 1.0f};
constexpr int k_maxSteps = 1001;

bool isAngle(const ElementParameter& p) { return p.kind == ParameterKind::Angle; }
double toDisplay(const ElementParameter& p, double value) { return isAngle(p) ? value * k_degrees : value; }
double fromDisplay(const ElementParameter& p, double value) { return isAngle(p) ? value / k_degrees : value; }

std::string formatValue(const ElementParameter& p, double value) {
    return isAngle(p) ? std::format("{:.4g}°", value * k_degrees) : std::format("{:.4g}", value);
}

// The smallest 1, 2 or 5 × 10ⁿ that is at least x (> 0).
double niceCeil(double x) {
    const double scale = std::pow(10.0, std::floor(std::log10(x)));
    for (const double m : {1.0, 2.0, 5.0}) {
        if (m * scale >= x * (1.0 - 1e-12)) return m * scale;
    }
    return 10.0 * scale;
}

// Slider limits for a parameter showing `value` (see parameter_panel.hpp).
std::pair<double, double> defaultRange(ParameterKind kind, double value, const Rect& view) {
    std::pair<double, double> r;
    switch (kind) {
    case ParameterKind::Speed: r = {0.0, niceCeil(2.0 * std::max(value, 1.0))}; break;
    case ParameterKind::Strength: {
        const double limit = niceCeil(2.0 * std::max(std::abs(value), 1.0));
        r = {-limit, limit};
        break;
    }
    case ParameterKind::Angle: {
        const double limit = std::abs(value) <= pi * (1.0 + 1e-12) ? pi : 2.0 * pi;
        r = {-limit, limit};
        break;
    }
    case ParameterKind::X: r = {view.min.x, view.max.x}; break;
    case ParameterKind::Y: r = {view.min.y, view.max.y}; break;
    }
    return {std::min(r.first, value), std::max(r.second, value)};
}

std::vector<ElementKind> kindsOf(const Scene& scene) {
    std::vector<ElementKind> kinds;
    for (const SceneElement& e : scene.elements) kinds.push_back(kindOf(e));
    return kinds;
}

} // namespace

float parameterPanelWidth(float available) {
    const float s = uiScale();
    return std::clamp(0.24f * available, 280.0f * s, 420.0f * s);
}

void ParameterPanel::reset() {
    m_kinds.clear();
    m_ranges.clear();
    m_sweep.reset();
    m_playing = false;
}

void ParameterPanel::sync(const Scene& scene) {
    std::vector<ElementKind> kinds = kindsOf(scene);
    // Elements added at the end leave the existing ones (and their ranges and sweep) as they were.
    const bool appended = kinds.size() >= m_kinds.size() && std::equal(m_kinds.begin(), m_kinds.end(), kinds.begin());
    if (!appended) reset();
    m_kinds = std::move(kinds);
}

std::optional<ElementEdit> ParameterPanel::advance(const Scene& scene, double dt) {
    sync(scene);
    if (!m_playing) return std::nullopt;
    const std::optional<SweepSetup> setup = sweep(scene);
    if (!setup) {
        m_playing = false;
        return std::nullopt;
    }
    m_time += dt;
    double f = 1.0; // fraction of the way from From to To
    if (m_once) {
        if (m_onceSeconds > 0.0) f = std::min(m_time / m_onceSeconds, 1.0);
        if (f >= 1.0) m_playing = false;
    } else if (setup->seconds > 0.0) {
        const double passes = m_time / setup->seconds;
        f = setup->backAndForth ? 1.0 - std::abs(std::fmod(passes, 2.0) - 1.0) : std::fmod(passes, 1.0);
    }
    const ParameterRef ref = setup->parameter;
    const double value = setup->from + f * (setup->to - setup->from);
    if (findParameter(scene, ref)->get(scene.elements[ref.element]) == value) return std::nullopt;
    std::optional<SceneElement> element = withParameter(scene, ref, value);
    if (!element) return std::nullopt; // e.g. a sector wall swept past the other one
    return ElementEdit{ref.element, std::move(*element)};
}

void ParameterPanel::playOnce(double seconds) {
    if (!m_sweep) return;
    m_playing = true;
    m_once = true;
    m_onceSeconds = seconds;
    m_time = 0.0;
}

std::optional<SweepSetup> ParameterPanel::sweep(const Scene& scene) const {
    if (!m_sweep) return std::nullopt;
    const ParameterRef ref = m_sweep->parameter;
    if (findParameter(scene, ref) == nullptr || ref.element >= m_kinds.size() ||
        kindOf(scene.elements[ref.element]) != m_kinds[ref.element]) {
        return std::nullopt;
    }
    return m_sweep;
}

std::string ParameterPanel::sweepDescription(const Scene& scene) const {
    const std::optional<SweepSetup> setup = sweep(scene);
    if (!setup) return {};
    const ElementParameter& p = *findParameter(scene, setup->parameter);
    return std::format("{} from {} to {}", parameterName(scene, setup->parameter), formatValue(p, setup->from),
                       formatValue(p, setup->to));
}

ParameterPanel::Range& ParameterPanel::range(ParameterRef ref, const ElementParameter& parameter, double value,
                                             const Rect& view) {
    const auto [it, added] = m_ranges.try_emplace({ref.element, ref.parameter});
    Range& r = it->second;
    if (added || value < r.low || value > r.high) {
        const auto [low, high] = defaultRange(parameter.kind, value, view);
        r = {added ? low : std::min(low, r.low), added ? high : std::max(high, r.high)};
    }
    return r;
}

void ParameterPanel::chooseSweep(ParameterRef ref, double value, const Range& range) {
    if (m_sweep && m_sweep->parameter == ref) return;
    SweepSetup setup = m_sweep.value_or(SweepSetup{}); // keeps steps, duration and back-and-forth
    setup.parameter = ref;
    setup.from = value;
    setup.to = range.high - value >= value - range.low ? range.high : range.low;
    if (setup.to == setup.from) setup.to = setup.from + 1.0;
    m_sweep = setup;
    m_playing = false;
}

bool ParameterPanel::drawSlider(ParameterRef ref, const ElementParameter& parameter, SceneElement& element,
                                const Rect& view) {
    const bool swept = m_sweep && m_sweep->parameter == ref;
    const double value = parameter.get(element);
    const Range& r = range(ref, parameter, value, view);

    ImGui::PushID(static_cast<int>(ref.parameter));
    ImGui::AlignTextToFramePadding();
    ImGui::TextColored(swept ? k_sweepColor : ImGui::GetStyleColorVec4(ImGuiCol_Text), "%s", parameter.label);
    ImGui::SetItemTooltip("%s", parameter.help);
    ImGui::SameLine(ImGui::GetFontSize() * 2.2f);
    ImGui::SetNextItemWidth(-FLT_MIN);
    double shown = toDisplay(parameter, value);
    const double low = toDisplay(parameter, r.low);
    const double high = toDisplay(parameter, r.high);
    ImGui::BeginDisabled(swept && m_playing);
    const bool changed = ImGui::SliderScalar("##value", ImGuiDataType_Double, &shown, &low, &high,
                                             isAngle(parameter) ? "%.1f°" : "%.4g");
    ImGui::EndDisabled();
    ImGui::SetItemTooltip("%s\n%s", parameter.help,
                          swept && m_playing ? "The sweep is moving it; stop the sweep to set it here."
                                             : "Ctrl+click to type a value. Right-click to sweep it.");
    if (ImGui::BeginPopupContextItem("##sweep")) {
        if (ImGui::MenuItem("Sweep this parameter", nullptr, swept)) chooseSweep(ref, value, r);
        ImGui::EndPopup();
    }
    ImGui::PopID();
    if (changed) parameter.set(element, fromDisplay(parameter, shown));
    return changed;
}

ParameterPanelResult ParameterPanel::draw(ImVec2 pos, ImVec2 size, const Scene& scene, const Rect& view,
                                          std::optional<float> exportProgress) {
    sync(scene);
    ParameterPanelResult result;
    ImGui::SetNextWindowPos(pos);
    ImGui::SetNextWindowSize(size);
    constexpr ImGuiWindowFlags flags = ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoCollapse |
                                       ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoDocking;
    if (!ImGui::Begin("Parameters", nullptr, flags)) {
        ImGui::End();
        return result;
    }

    if (scene.elements.empty()) ImGui::TextWrapped("No elements yet. Add one with the toolbar on the left.");
    for (std::size_t i = 0; i < scene.elements.size(); ++i) {
        const ElementKind kind = kindOf(scene.elements[i]);
        ImGui::SetNextItemOpen(true, ImGuiCond_Once);
        const bool open = ImGui::CollapsingHeader(std::format("#{}  {}###element{}", i + 1, kindName(kind), i).c_str());
        const ImVec2 min = ImGui::GetItemRectMin();
        const ImVec2 max = ImGui::GetItemRectMax();
        const float h = max.y - min.y;
        drawToolIcon(ImGui::GetWindowDrawList(), ImVec2(max.x - 0.8f * h, 0.5f * (min.y + max.y)), 0.75f * h, kind);
        if (!open) continue;

        ImGui::PushID(static_cast<int>(i));
        SceneElement element = scene.elements[i];
        bool changed = false;
        const std::span<const ElementParameter> parameters = elementParameters(kind);
        for (std::size_t p = 0; p < parameters.size(); ++p) changed |= drawSlider({i, p}, parameters[p], element, view);
        ImGui::PopID();
        // Values that make the element invalid (a sector's θ₂ dragged below θ₁, ...) are not applied.
        if (changed && !sceneElementError(element)) result.edits.push_back({i, std::move(element)});
    }

    ImGui::Spacing();
    ImGui::SeparatorText("Sweep");
    drawSweep(scene, view, exportProgress, result);
    ImGui::End();
    return result;
}

void ParameterPanel::drawSweep(const Scene& scene, const Rect& view, std::optional<float> exportProgress,
                               ParameterPanelResult& result) {
    const std::optional<SweepSetup> setup = sweep(scene);
    ImGui::SetNextItemWidth(-FLT_MIN);
    if (ImGui::BeginCombo("##sweep-parameter", setup ? parameterName(scene, setup->parameter).c_str() : "Choose a parameter...")) {
        for (std::size_t i = 0; i < scene.elements.size(); ++i) {
            const std::span<const ElementParameter> parameters = elementParameters(kindOf(scene.elements[i]));
            for (std::size_t p = 0; p < parameters.size(); ++p) {
                const ParameterRef ref{i, p};
                const bool selected = setup && setup->parameter == ref;
                if (ImGui::Selectable(parameterName(scene, ref).c_str(), selected)) {
                    const double value = parameters[p].get(scene.elements[i]);
                    chooseSweep(ref, value, range(ref, parameters[p], value, view));
                }
            }
        }
        ImGui::EndCombo();
    }
    if (!m_sweep || !setup) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
        ImGui::TextWrapped("Pick the parameter to sweep, or right-click its slider.");
        ImGui::PopStyleColor();
        return;
    }

    SweepSetup& s = *m_sweep;
    const ElementParameter& parameter = *findParameter(scene, s.parameter);
    const float labelWidth = ImGui::GetFontSize() * 4.2f;
    const auto label = [labelWidth](const char* text) {
        ImGui::AlignTextToFramePadding();
        ImGui::TextUnformatted(text);
        ImGui::SameLine(labelWidth);
        ImGui::SetNextItemWidth(-FLT_MIN);
    };
    const char* format = isAngle(parameter) ? "%.4g°" : "%.6g";
    double from = toDisplay(parameter, s.from);
    label("From");
    if (ImGui::InputDouble("##from", &from, 0.0, 0.0, format)) s.from = fromDisplay(parameter, from);
    double to = toDisplay(parameter, s.to);
    label("To");
    if (ImGui::InputDouble("##to", &to, 0.0, 0.0, format)) s.to = fromDisplay(parameter, to);

    label("Duration");
    float seconds = static_cast<float>(s.seconds);
    if (ImGui::SliderFloat("##duration", &seconds, 1.0f, 30.0f, "%.1f s")) s.seconds = std::max(seconds, 0.1f);
    ImGui::SetItemTooltip("Seconds of animation for one pass from From to To while playing.");
    ImGui::Checkbox("Back and forth", &s.backAndForth);
    ImGui::SetItemTooltip("Play From → To → From ...; otherwise From → To, then start again at From.");

    if (ImGui::Button(m_playing ? "Stop" : "Play", ImVec2(-FLT_MIN, 0.0f))) {
        m_playing = !m_playing;
        m_once = false;
        m_time = 0.0;
    }
    ImGui::SetItemTooltip("Animate the parameter: the flow follows live. Stopping leaves it where it is.");

    ImGui::Spacing();
    label("Steps");
    if (ImGui::InputInt("##steps", &s.steps)) s.steps = std::clamp(s.steps, 2, k_maxSteps);
    ImGui::SetItemTooltip("Values in the CSV export, From and To included.");
    if (exportProgress) {
        const float buttonWidth = 0.65f * ImGui::GetContentRegionAvail().x;
        ImGui::ProgressBar(*exportProgress, ImVec2(buttonWidth, 0.0f));
        ImGui::SameLine();
        if (ImGui::Button("Cancel", ImVec2(-FLT_MIN, 0.0f))) result.cancelExport = true;
    } else {
        if (ImGui::Button("Export CSV...", ImVec2(-FLT_MIN, 0.0f))) result.exportSweep = true;
        ImGui::SetItemTooltip("At each step: the stagnation points in the view and every body found, with its lift,\n"
                              "drag, c_l, c_d and lowest Cp (and the Kutta-Joukowski lift for an isolated body).\n"
                              "Computed one step per frame over the current view.");
    }
}

} // namespace pfl
