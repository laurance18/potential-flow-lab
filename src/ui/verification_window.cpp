#define IMGUI_DEFINE_MATH_OPERATORS
#include "ui/verification_window.hpp"

#include "core/flow_field.hpp"
#include "ui/style.hpp"

#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <exception>
#include <format>
#include <limits>
#include <numbers>
#include <span>

namespace pfl {
namespace {

constexpr const char* k_windowTitle = "Physics verification";
constexpr double k_minCheckSeconds = 0.14; // each check shows as running at least this long
constexpr double k_popSeconds = 0.25;      // a finished check's icon pops in
constexpr double k_meterSeconds = 0.7;     // and its bar fills
constexpr double k_maxDigits = 16.0;       // a full bar: agreement to double precision
constexpr double k_machinePrecision = 1e-13;
constexpr double k_tracerSpeed = 1.4;      // badge radii per second in the free stream
constexpr double k_trailInterval = 1.0 / 30.0;
constexpr double k_celebrationCirculation = 1.8 * std::numbers::pi; // stagnation points 26° below the axis
constexpr float k_pi = std::numbers::pi_v<float>;

constexpr ImU32 k_passColor = IM_COL32(95, 205, 120, 255);
constexpr ImU32 k_failColor = IM_COL32(235, 95, 85, 255);
constexpr ImU32 k_runColor = IM_COL32(110, 170, 245, 255);
constexpr ImU32 k_pendingColor = IM_COL32(255, 255, 255, 70);
constexpr ImU32 k_trackColor = IM_COL32(255, 255, 255, 24);
constexpr ImU32 k_tickColor = IM_COL32(255, 255, 255, 190);
constexpr ImU32 k_bandColor = IM_COL32(11, 13, 19, 255);
constexpr ImU32 k_badgeColor = IM_COL32(24, 27, 37, 255);
constexpr ImU32 k_iconInkColor = IM_COL32(14, 26, 18, 255);
constexpr ImVec4 k_idleTracer{0.62f, 0.72f, 0.88f, 1.0f};
constexpr ImVec4 k_passTracer{0.45f, 0.86f, 0.56f, 1.0f};
constexpr ImVec4 k_failTracer{0.93f, 0.47f, 0.42f, 1.0f};

float clamp01(double t) { return static_cast<float>(std::clamp(t, 0.0, 1.0)); }
float easeOut(double t) {
    const float u = 1.0f - clamp01(t);
    return 1.0f - u * u * u;
}
// Overshoots a little before settling: a pop.
float easeOutBack(double t) {
    constexpr float c1 = 1.70158f;
    const float u = clamp01(t) - 1.0f;
    return 1.0f + (c1 + 1.0f) * u * u * u + c1 * u * u;
}

ImVec4 mix(ImVec4 a, ImVec4 b, float t) {
    return {a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, a.z + (b.z - a.z) * t, a.w + (b.w - a.w) * t};
}

double digitsOf(double error) { return error <= 0.0 ? k_maxDigits : std::clamp(-std::log10(error), 0.0, k_maxDigits); }

std::string errorText(const VerificationResult& r) {
    if (!std::isfinite(r.error)) return "failed";
    if (r.error == 0.0) return "exact";
    return std::format("{:.1e}", r.error);
}

std::string toleranceText(double tolerance) { return std::format("{:.0e}", tolerance); }

float screenLength(ImVec2 v) { return std::sqrt(v.x * v.x + v.y * v.y); }

// A tick in a box of side `size` round `c`, drawn up to `progress` (0 to 1) of its length.
void drawTick(ImDrawList* dl, ImVec2 c, float size, float progress, ImU32 colour, float thickness) {
    if (progress <= 0.0f) return;
    const ImVec2 a = c + ImVec2(-0.34f, 0.02f) * size;
    const ImVec2 b = c + ImVec2(-0.10f, 0.26f) * size;
    const ImVec2 d = c + ImVec2(0.36f, -0.24f) * size;
    const float first = screenLength(b - a);
    const float along = progress * (first + screenLength(d - b));
    dl->PathLineTo(a);
    if (along <= first) {
        dl->PathLineTo(a + (b - a) * (along / first));
    } else {
        dl->PathLineTo(b);
        dl->PathLineTo(b + (d - b) * ((along - first) / screenLength(d - b)));
    }
    dl->PathStroke(colour, ImDrawFlags_None, thickness);
}

void drawCross(ImDrawList* dl, ImVec2 c, float size, ImU32 colour, float thickness) {
    const float h = 0.27f * size;
    dl->AddLine(c + ImVec2(-h, -h), c + ImVec2(h, h), colour, thickness);
    dl->AddLine(c + ImVec2(-h, h), c + ImVec2(h, -h), colour, thickness);
}

// The status of one check: a faint ring while it waits, a spinner while it runs, then a green
// disc with a tick (or a red one with a cross) that pops in with a ripple.
void drawStatus(ImDrawList* dl, ImVec2 c, float r, const std::optional<VerificationResult>& result, bool active,
                double since, double now) {
    const float s = uiScale();
    if (!result) {
        if (!active) {
            dl->AddCircle(c, 0.62f * r, k_pendingColor, 0, 1.5f * s);
            return;
        }
        const float a0 = static_cast<float>(std::fmod(now * 5.0, 2.0 * std::numbers::pi));
        dl->AddCircle(c, 0.72f * r, k_trackColor, 0, 2.2f * s);
        dl->PathArcTo(c, 0.72f * r, a0, a0 + 1.35f * k_pi, 20);
        dl->PathStroke(k_runColor, ImDrawFlags_None, 2.2f * s);
        return;
    }
    const double t = now - since;
    const ImU32 colour = result->passed ? k_passColor : k_failColor;
    if (t < 0.6) dl->AddCircle(c, r * (0.9f + 0.9f * easeOut(t / 0.6)), fadeColor(colour, 1.0f - clamp01(t / 0.6)), 0, 2.0f * s);
    dl->AddCircleFilled(c, 0.82f * r * easeOutBack(t / k_popSeconds), colour);
    if (result->passed) {
        drawTick(dl, c, 1.5f * r, easeOut((t - 0.08) / k_popSeconds), k_iconInkColor, 2.3f * s);
    } else if (t > 0.08) {
        drawCross(dl, c, 1.5f * r, IM_COL32(40, 12, 12, 255), 2.3f * s);
    }
}

// The error on a log scale: the bar fills to -log10(error) of 16 digits, and a tick marks the
// tolerance. While the check runs a highlight sweeps along it.
void drawMeter(ImDrawList* dl, ImVec2 min, ImVec2 size, const std::optional<VerificationResult>& result, bool active,
               double since, double now) {
    const float s = uiScale();
    const float radius = 0.5f * size.y;
    const ImVec2 max = min + size;
    dl->AddRectFilled(min, max, k_trackColor, radius);
    if (!result) {
        if (active) {
            const float x = static_cast<float>(std::fmod(now * 0.9, 1.0));
            const float left = min.x + size.x * std::max(0.0f, x - 0.25f);
            const float right = min.x + size.x * std::min(1.0f, x + 0.05f);
            dl->AddRectFilled(ImVec2(left, min.y), ImVec2(right, max.y), fadeColor(k_runColor, 0.45f), radius);
        }
        return;
    }
    const ImU32 colour = result->passed ? k_passColor : k_failColor;
    const float fill = static_cast<float>(digitsOf(result->error) / k_maxDigits) * easeOut((now - since) / k_meterSeconds);
    if (fill > 0.0f) dl->AddRectFilled(min, ImVec2(min.x + std::max(size.x * fill, size.y), max.y), colour, radius);
    if (result->tolerance > 0.0) {
        const float x = min.x + size.x * static_cast<float>(digitsOf(result->tolerance) / k_maxDigits);
        dl->AddLine(ImVec2(x, min.y - 3.0f * s), ImVec2(x, max.y + 3.0f * s), k_tickColor, 1.5f * s);
    }
}

void centredText(const char* text, ImU32 colour) {
    const float width = ImGui::CalcTextSize(text).x;
    ImGui::SetCursorPosX(std::max(ImGui::GetCursorPosX(), 0.5f * (ImGui::GetWindowWidth() - width)));
    ImGui::PushStyleColor(ImGuiCol_Text, colour);
    ImGui::TextUnformatted(text);
    ImGui::PopStyleColor();
}

void wrappedDim(const std::string& text) {
    ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
    ImGui::TextWrapped("%s", text.c_str());
    ImGui::PopStyleColor();
}

} // namespace

VerificationWindow::VerificationWindow() : m_entries(verificationChecks().size()), m_current(m_entries.size()) {}

// A check still running finishes first (the future's destructor waits); none takes long.
VerificationWindow::~VerificationWindow() = default;

bool VerificationWindow::complete() const {
    return std::ranges::all_of(m_entries, [](const Entry& e) { return e.result.has_value(); });
}

int VerificationWindow::passedCount() const {
    return static_cast<int>(std::ranges::count_if(m_entries, [](const Entry& e) { return e.result && e.result->passed; }));
}

void VerificationWindow::open() {
    m_open = true;
    ImGui::SetWindowFocus(k_windowTitle);
    if (!running() && !complete()) start();
}

void VerificationWindow::start() {
    for (Entry& e : m_entries) e = Entry{};
    m_shownProgress = 0.0;
    const double now = ImGui::GetTime();
    m_runStartedAt = now;
    launch(0, now);
}

void VerificationWindow::launch(std::size_t index, double now) {
    m_current = index;
    m_checkStartedAt = now;
    m_followRunning = true;
    m_pending = std::async(std::launch::async, verificationChecks()[index].run);
}

void VerificationWindow::poll(double now) {
    if (!running() || now - m_checkStartedAt < k_minCheckSeconds) return;
    if (m_pending.wait_for(std::chrono::seconds(0)) != std::future_status::ready) return;
    Entry& entry = m_entries[m_current];
    try {
        entry.result = m_pending.get();
    } catch (const std::exception& e) {
        entry.result = VerificationResult{.error = std::numeric_limits<double>::infinity(), .note = e.what()};
    }
    entry.finishedAt = now;
    if (!entry.result->passed) entry.expanded = true; // show straight away what went wrong
    if (m_current + 1 < m_entries.size()) {
        launch(m_current + 1, now);
    } else {
        m_current = m_entries.size();
        m_runFinishedAt = now;
    }
}

double VerificationWindow::random() {
    m_seed = m_seed * 1'664'525u + 1'013'904'223u;
    return static_cast<double>(m_seed >> 8) / static_cast<double>(1u << 24);
}

std::optional<VerificationShow> VerificationWindow::draw() {
    const double now = ImGui::GetTime();
    const double dt = m_lastTime < 0.0 ? 0.0 : std::min(now - m_lastTime, 1.0 / 20.0);
    m_lastTime = now;
    if (!m_open) {
        // Closing abandons the run; the check in progress finishes on its own.
        m_current = m_entries.size();
        return std::nullopt;
    }
    poll(now);

    const float s = uiScale();
    // Inside the area below the menu bar, so the window never hides the menus.
    const ImGuiViewport* viewport = ImGui::GetMainViewport();
    const ImVec2 work = viewport->WorkSize;
    ImGui::SetNextWindowSize(ImVec2(std::min(700.0f * s, 0.9f * work.x), std::min(780.0f * s, work.y - 24.0f * s)),
                             ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowPos(viewport->WorkPos + work * 0.5f, ImGuiCond_FirstUseEver, ImVec2(0.5f, 0.5f));    ImGui::SetNextWindowSizeConstraints(ImVec2(520.0f * s, 460.0f * s), ImVec2(FLT_MAX, FLT_MAX));
    std::optional<VerificationShow> show;
    if (ImGui::Begin(k_windowTitle, &m_open, ImGuiWindowFlags_NoCollapse)) {
        drawHeader(now, dt);
        drawSummary(now);
        ImGui::Separator();
        show = drawChecks(now);
    }
    ImGui::End();
    return show;
}

void VerificationWindow::drawHeader(double now, double dt) {
    const float s = uiScale();
    const ImVec2 origin = ImGui::GetCursorScreenPos();
    const ImVec2 size(ImGui::GetContentRegionAvail().x, 136.0f * s);
    ImGui::Dummy(size);
    if (size.x <= 0.0f) return;

    const float unit = 34.0f * s; // the badge's radius, and the flow's length scale
    const double halfWidth = 0.5 * size.x / unit;
    const double halfHeight = 0.5 * size.y / unit;
    const bool done = complete() && !running();
    const bool allPassed = done && passedCount() == static_cast<int>(m_entries.size());

    const double progress = static_cast<double>(std::ranges::count_if(m_entries, [](const Entry& e) { return e.result.has_value(); })) /
                            static_cast<double>(m_entries.size());
    m_shownProgress += (progress - m_shownProgress) * std::min(1.0, dt * 8.0);
    const double circulation = allPassed ? k_celebrationCirculation : 0.0;
    m_circulation += (circulation - m_circulation) * std::min(1.0, dt * 1.2);
    updateTracers(dt, halfWidth, halfHeight);

    const float celebrate = done ? easeOut((now - m_runFinishedAt) / 1.2) : 0.0f;
    const ImVec4 tint = mix(k_idleTracer, allPassed ? k_passTracer : k_failTracer, celebrate);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 max = origin + size;
    dl->PushClipRect(origin, max, true);
    dl->AddRectFilled(origin, max, k_bandColor, 6.0f * s);
    const ImVec2 centre = origin + size * 0.5f;
    drawTracers(dl, centre, unit, halfWidth, ImGui::ColorConvertFloat4ToU32(tint));
    drawBadge(dl, centre, unit, now);
    dl->PopClipRect();
    if (ImGui::IsItemHovered()) {
        ImGui::SetTooltip("Flow past a cylinder (uniform flow + doublet), computed by the program's own solver.\n"
                          "Once every check has passed it picks up circulation: the stagnation points move down\n"
                          "and the cylinder carries lift.");
    }
}

void VerificationWindow::updateTracers(double dt, double halfWidth, double halfHeight) {
    const FlowField field{.elements = {UniformFlow{}, Doublet{.strength = 1.0}, Vortex{.circulation = m_circulation}}};
    const auto spawn = [&](Tracer& t, bool anywhere) {
        t = Tracer{};
        for (int attempt = 0; attempt < 20; ++attempt) {
            t.position = {anywhere ? (2.0 * random() - 1.0) * halfWidth : -halfWidth - 0.2 - 0.6 * random(),
                          (2.0 * random() - 1.0) * halfHeight};
            if (length(t.position) > 1.08) break;
        }
        t.age = anywhere ? 1.0 : 0.0; // fresh ones fade in; the first fill is already there
    };
    const std::size_t target = static_cast<std::size_t>(std::clamp(5.0 * halfWidth * halfHeight, 40.0, 280.0));
    while (m_tracers.size() < target) spawn(m_tracers.emplace_back(), true);
    m_tracers.resize(target);

    m_trailClock += dt;
    const bool sampleTrail = m_trailClock >= k_trailInterval;
    if (sampleTrail) m_trailClock = std::fmod(m_trailClock, k_trailInterval);
    for (Tracer& t : m_tracers) {
        if (sampleTrail) {
            if (t.trailCount < k_trailLength) {
                t.trail[t.trailCount++] = t.position;
            } else {
                std::shift_left(t.trail.begin(), t.trail.end(), 1);
                t.trail.back() = t.position;
            }
        }
        // Midpoint rule; the flow is smooth outside the badge.
        const double h = dt * k_tracerSpeed;
        const Vec2 mid = t.position + field.velocity(t.position) * (0.5 * h);
        t.position += field.velocity(mid) * h;
        if (const double r = length(t.position); r < 1.0) t.position = t.position * (1.0001 / r);
        t.age += dt;
        if (t.position.x > halfWidth + 0.3 || std::abs(t.position.y) > halfHeight + 0.6 || t.age > 30.0) spawn(t, false);
    }
}

void VerificationWindow::drawTracers(ImDrawList* dl, ImVec2 centre, float unit, double halfWidth, ImU32 colour) const {
    const float s = uiScale();
    const auto screen = [&](Vec2 p) { return centre + ImVec2(static_cast<float>(p.x) * unit, -static_cast<float>(p.y) * unit); };
    for (const Tracer& t : m_tracers) {
        // Fade in after spawning and towards both ends of the band.
        const float alpha = clamp01(t.age / 0.8) * clamp01((halfWidth + 0.2 - std::abs(t.position.x)) / 1.2);
        if (alpha <= 0.0f) continue;
        ImVec2 previous = screen(t.position);
        for (std::size_t k = t.trailCount; k-- > 0;) {
            const ImVec2 p = screen(t.trail[k]);
            const float fade = 1.0f - static_cast<float>(t.trailCount - k) / static_cast<float>(k_trailLength + 1);
            dl->AddLine(previous, p, fadeColor(colour, 0.5f * alpha * fade), 1.4f * s);
            previous = p;
        }
        dl->AddCircleFilled(screen(t.position), 1.6f * s, fadeColor(colour, 0.95f * alpha));
    }
}

void VerificationWindow::drawBadge(ImDrawList* dl, ImVec2 centre, float radius, double now) const {
    const float s = uiScale();
    const bool done = complete() && !running();
    const int total = static_cast<int>(m_entries.size());
    const int passed = passedCount();
    const bool allPassed = done && passed == total;
    const double since = now - m_runFinishedAt;

    dl->AddCircleFilled(centre, radius, k_badgeColor, 64);
    dl->AddCircle(centre, radius, IM_COL32(255, 255, 255, 40), 64, 1.0f * s);
    const float ring = radius - 6.0f * s;
    dl->AddCircle(centre, ring, k_trackColor, 64, 4.0f * s);
    const ImU32 ringColour = !done ? k_runColor : allPassed ? k_passColor : k_failColor;
    const float sweep = 2.0f * k_pi * static_cast<float>(m_shownProgress);
    if (sweep > 0.01f) {
        dl->PathArcTo(centre, ring, -0.5f * k_pi, -0.5f * k_pi + sweep, 64);
        dl->PathStroke(ringColour, ImDrawFlags_None, 4.0f * s);
    }
    if (done && since < 1.0) {
        // A ripple leaves the badge as the run completes.
        dl->AddCircle(centre, radius * (1.0f + 1.4f * easeOut(since)), fadeColor(ringColour, 0.8f * (1.0f - clamp01(since))), 64,
                      2.5f * s);
    }

    if (allPassed) {
        drawTick(dl, centre, 1.25f * radius, easeOut((since - 0.1) / 0.45), k_passColor, 4.5f * s);
        return;
    }
    // While it runs, and after a failure: how many passed.
    const std::string text = std::format("{}/{}", passed, total);
    ImFont* font = ImGui::GetFont();
    const float size = 1.25f * ImGui::GetFontSize();
    const ImVec2 extent = font->CalcTextSizeA(size, FLT_MAX, 0.0f, text.c_str());
    dl->AddText(font, size, centre - extent * 0.5f, done ? k_failColor : IM_COL32(230, 235, 245, 255), text.c_str());
}

void VerificationWindow::drawSummary(double now) {
    const float s = uiScale();
    const std::span<const VerificationCheck> checks = verificationChecks();
    const int total = static_cast<int>(m_entries.size());
    const int passed = passedCount();
    const bool done = complete() && !running();

    ImGui::Spacing();
    std::string title;
    ImU32 colour = IM_COL32(235, 238, 245, 255);
    if (running()) {
        title = "Checking the physics…";
    } else if (!done) {
        title = "Stopped before the end";
    } else if (passed == total) {
        title = std::format("All {} checks passed", total);
        colour = ImGui::ColorConvertFloat4ToU32(mix(ImVec4(0.92f, 0.93f, 0.96f, 1.0f), ImVec4(0.55f, 0.88f, 0.62f, 1.0f),
                                                    easeOut((now - m_runFinishedAt) / 0.8)));
    } else {
        title = std::format("{} of {} checks failed", total - passed, total);
        colour = k_failColor;
    }
    ImGui::PushFont(nullptr, ImGui::GetStyle().FontSizeBase * 1.45f);
    centredText(title.c_str(), colour);
    ImGui::PopFont();

    std::string line;
    if (running()) {
        line = std::format("Check {} of {}: {}", m_current + 1, total, checks[m_current].title);
    } else if (done) {
        int exact = 0;
        long long samples = 0;
        for (const Entry& e : m_entries) {
            exact += e.result->error < k_machinePrecision;
            samples += e.result->samples;
        }
        line = passed == total
                   ? std::format("{} agree with theory to machine precision; the rest to the accuracy of their integration.",
                                 exact == total ? "All" : std::format("{} of {}", exact, total))
                   : "Open the failed checks below to see where the program and theory part ways.";
        centredText(line.c_str(), ImGui::GetColorU32(ImGuiCol_TextDisabled));
        std::string grouped = std::to_string(samples);
        for (int i = static_cast<int>(grouped.size()) - 3; i > 0; i -= 3) grouped.insert(static_cast<std::size_t>(i), ",");
        line = std::format("{} values compared in {:.1f} s", grouped, m_runFinishedAt - m_runStartedAt);
    } else {
        line = "Closed before every check had run.";
    }
    centredText(line.c_str(), ImGui::GetColorU32(ImGuiCol_TextDisabled));

    ImGui::Spacing();
    const ImVec2 button(150.0f * s, 0.0f);
    const float buttonsWidth = 2.0f * button.x + ImGui::GetStyle().ItemSpacing.x;
    ImGui::SetCursorPosX(std::max(ImGui::GetCursorPosX(), 0.5f * (ImGui::GetWindowWidth() - buttonsWidth)));
    ImGui::BeginDisabled(running());
    if (ImGui::Button(running() ? "Running…" : done ? "Run again" : "Run", button)) start();
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::BeginDisabled(!done);
    if (ImGui::Button(now - m_copiedAt < 1.5 ? "Copied" : "Copy report", button)) {
        ImGui::SetClipboardText(report().c_str());
        m_copiedAt = now;
    }
    ImGui::EndDisabled();
    ImGui::SetItemTooltip("Copy every check, with its numbers, as plain text.");
    ImGui::Spacing();
    ImGui::PushTextWrapPos(0.0f);
    wrappedDim("Each check builds a flow with the program and compares it with the exact result from potential-flow "
               "theory. The bars show −log₁₀ of the error, about the number of digits that agree; the white tick marks "
               "the tolerance. Click a check for its numbers.");
    ImGui::PopTextWrapPos();
    ImGui::Spacing();
}

std::optional<VerificationShow> VerificationWindow::drawChecks(double now) {
    std::optional<VerificationShow> show;
    if (!ImGui::BeginChild("##checks", ImVec2(0.0f, 0.0f), ImGuiChildFlags_None)) {
        ImGui::EndChild();
        return show;
    }
    const float s = uiScale();
    const std::span<const VerificationCheck> checks = verificationChecks();
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const float rowHeight = 1.3f * ImGui::GetFrameHeight();
    const float meterWidth = 130.0f * s;
    const float errorWidth = ImGui::CalcTextSize("8.8e-88").x + 10.0f * s;
    // Open checks and the hovered one are lit faintly, so the status colours stay the loudest thing.
    ImGui::PushStyleColor(ImGuiCol_Header, IM_COL32(255, 255, 255, 14));
    ImGui::PushStyleColor(ImGuiCol_HeaderHovered, IM_COL32(255, 255, 255, 24));
    ImGui::PushStyleColor(ImGuiCol_HeaderActive, IM_COL32(255, 255, 255, 34));

    for (std::size_t i = 0; i < checks.size(); ++i) {
        const VerificationCheck& check = checks[i];
        Entry& entry = m_entries[i];
        if (i == 0 || std::string_view(check.group) != checks[i - 1].group) {
            if (i > 0) ImGui::Spacing();
            int inGroup = 0;
            int passedInGroup = 0;
            for (std::size_t k = 0; k < checks.size(); ++k) {
                if (std::string_view(checks[k].group) != check.group) continue;
                ++inGroup;
                passedInGroup += m_entries[k].result && m_entries[k].result->passed;
            }
            const float right = ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x;
            ImGui::TextDisabled("%s", check.group);
            const std::string tally = std::format("{} / {}", passedInGroup, inGroup);
            ImGui::SameLine(right - ImGui::CalcTextSize(tally.c_str()).x);
            ImGui::PushStyleColor(ImGuiCol_Text, passedInGroup == inGroup ? k_passColor : ImGui::GetColorU32(ImGuiCol_TextDisabled));
            ImGui::TextUnformatted(tally.c_str());
            ImGui::PopStyleColor();
        }

        ImGui::PushID(static_cast<int>(i));
        const ImVec2 min = ImGui::GetCursorScreenPos();
        const bool active = m_current == i;
        if (ImGui::Selectable("##row", entry.expanded, ImGuiSelectableFlags_None, ImVec2(0.0f, rowHeight)) && entry.result) {
            entry.expanded = !entry.expanded;
        }
        const bool hovered = ImGui::IsItemHovered();
        const ImVec2 max = ImGui::GetItemRectMax();
        if (active && m_followRunning) {
            if (max.y > ImGui::GetWindowPos().y + ImGui::GetWindowHeight()) ImGui::SetScrollHereY(1.0f);
            m_followRunning = false;
        }

        const float midY = 0.5f * (min.y + max.y);
        const float iconRadius = 0.36f * rowHeight;
        drawStatus(dl, ImVec2(min.x + 6.0f * s + iconRadius, midY), iconRadius, entry.result, active, entry.finishedAt, now);

        // Right to left: the chevron, the error and the bar.
        const float chevronX = max.x - 10.0f * s;
        if (entry.result) {
            const float h = 3.5f * s;
            const ImU32 c = ImGui::GetColorU32(ImGuiCol_TextDisabled);
            if (entry.expanded) {
                dl->AddTriangleFilled(ImVec2(chevronX - h, midY - 0.5f * h), ImVec2(chevronX + h, midY - 0.5f * h),
                                      ImVec2(chevronX, midY + 0.7f * h), c);
            } else {
                dl->AddTriangleFilled(ImVec2(chevronX - 0.5f * h, midY - h), ImVec2(chevronX - 0.5f * h, midY + h),
                                      ImVec2(chevronX + 0.7f * h, midY), c);
            }
        }
        const float errorX = chevronX - 12.0f * s - errorWidth;
        if (entry.result) {
            const std::string text = errorText(*entry.result);
            const ImU32 c = entry.result->passed ? ImGui::GetColorU32(ImGuiCol_Text) : k_failColor;
            dl->AddText(ImVec2(errorX + errorWidth - ImGui::CalcTextSize(text.c_str()).x, midY - 0.5f * ImGui::GetFontSize()),
                        fadeColor(c, clamp01((now - entry.finishedAt) / 0.3)), text.c_str());
        }
        const ImVec2 meterMin(errorX - 10.0f * s - meterWidth, midY - 4.0f * s);
        const ImVec2 meterSize(meterWidth, 8.0f * s);
        drawMeter(dl, meterMin, meterSize, entry.result, active, entry.finishedAt, now);

        const float textX = min.x + 14.0f * s + 2.0f * iconRadius;
        const ImU32 textColour = entry.result || active ? ImGui::GetColorU32(ImGuiCol_Text) : ImGui::GetColorU32(ImGuiCol_TextDisabled);
        dl->PushClipRect(ImVec2(textX, min.y), ImVec2(meterMin.x - 8.0f * s, max.y), true);
        dl->AddText(ImVec2(textX, midY - 0.5f * ImGui::GetFontSize()), textColour, check.title);
        dl->PopClipRect();

        if (hovered && entry.result) {
            const VerificationResult& r = *entry.result;
            ImGui::SetTooltip("%s\n\nError %s, tolerance %s: %s.\nClick for the numbers.", check.theory, errorText(r).c_str(),
                              toleranceText(r.tolerance).c_str(), r.passed ? "passed" : "failed");
        } else if (hovered) {
            ImGui::SetTooltip("%s", check.theory);
        }

        if (entry.expanded && entry.result) {
            if (std::optional<VerificationShow> pressed = drawDetails(i)) show = std::move(pressed);
        }
        ImGui::PopID();
    }
    ImGui::PopStyleColor(3);
    ImGui::Spacing();
    ImGui::EndChild();
    return show;
}

std::optional<VerificationShow> VerificationWindow::drawDetails(std::size_t index) {
    const float s = uiScale();
    const VerificationCheck& check = verificationChecks()[index];
    const VerificationResult& r = *m_entries[index].result;
    std::optional<VerificationShow> show;

    const float indent = 14.0f * s + 0.72f * ImGui::GetFrameHeight() * 1.3f;
    ImGui::Indent(indent);
    ImGui::Spacing();
    ImGui::TextWrapped("%s", check.theory);
    wrappedDim(check.method);

    if (!r.comparisons.empty() &&
        ImGui::BeginTable("##values", 3, ImGuiTableFlags_SizingStretchProp | ImGuiTableFlags_BordersInnerH | ImGuiTableFlags_RowBg,
                          ImVec2(ImGui::GetContentRegionAvail().x - 10.0f * s, 0.0f))) {
        ImGui::TableSetupColumn("", ImGuiTableColumnFlags_WidthStretch, 1.5f);
        ImGui::TableSetupColumn("Program", ImGuiTableColumnFlags_WidthStretch, 1.0f);
        ImGui::TableSetupColumn("Theory", ImGuiTableColumnFlags_WidthStretch, 1.1f);
        ImGui::TableHeadersRow();
        for (const Comparison& c : r.comparisons) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::TextWrapped("%s", c.quantity.c_str());
            ImGui::TableNextColumn();
            ImGui::TextWrapped("%s", c.program.c_str());
            ImGui::TableNextColumn();
            ImGui::TextWrapped("%s", c.exact.c_str());
        }
        ImGui::EndTable();
    }

    ImGui::PushStyleColor(ImGuiCol_Text, r.passed ? k_passColor : k_failColor);
    if (std::isfinite(r.error)) {
        ImGui::Text("%s: error %s %s tolerance %s", r.passed ? "Passed" : "Failed", errorText(r).c_str(), r.passed ? "≤" : ">",
                    toleranceText(r.tolerance).c_str());
    } else {
        ImGui::TextUnformatted("Failed: the check could not be completed");
    }
    ImGui::PopStyleColor();
    if (!r.errorLabel.empty()) wrappedDim("Error: " + r.errorLabel + ".");
    if (!r.note.empty()) wrappedDim(r.note);

    if (ImGui::SmallButton("Show in plot")) show = VerificationShow{check.scene(), check.title};
    ImGui::SetItemTooltip("Load this flow into the plot, to look at it and play with it.\n"
                          "(You are asked first if the scene has unsaved changes.)");
    ImGui::Spacing();
    ImGui::Unindent(indent);
    return show;
}

std::string VerificationWindow::report() const {
    const std::span<const VerificationCheck> checks = verificationChecks();
    std::string out = std::format("Potential Flow Lab: physics verification\n{} of {} checks passed.\n", passedCount(),
                                  m_entries.size());
    for (std::size_t i = 0; i < checks.size(); ++i) {
        if (!m_entries[i].result) continue;
        const VerificationResult& r = *m_entries[i].result;
        out += std::format("\n[{}] {} / {}\n", r.passed ? "PASS" : "FAIL", checks[i].group, checks[i].title);
        out += std::format("  Theory: {}\n", checks[i].theory);
        out += std::format("  Error {} (tolerance {}): {}\n", errorText(r), toleranceText(r.tolerance), r.errorLabel);
        for (const Comparison& c : r.comparisons) out += std::format("  {}: program {}, theory {}\n", c.quantity, c.program, c.exact);
        if (!r.note.empty()) out += std::format("  {}\n", r.note);
    }
    return out;
}

} // namespace pfl
