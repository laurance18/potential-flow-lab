// verification_window.hpp: Verify → Physics checks.
//
// Runs every check of core/verification.hpp in turn, on a worker thread so that the window stays
// responsive, and shows each one as it finishes: a tick or a cross, and a bar for the error on a
// log scale (about the number of digits that agree with theory) with the tolerance marked. A
// check opens up to say what theory predicts and how the program was tested, to set the
// program's numbers beside the exact ones, and to load the flow into the plot. Above the list, a
// badge fills as the checks pass, with tracer particles flowing round it: the flow past a cylinder,
// computed by the program's own solver, that gains circulation once every check has passed.
#pragma once

#include "core/scene.hpp"
#include "core/vec2.hpp"
#include "core/verification.hpp"

#include <imgui.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <future>
#include <optional>
#include <string>
#include <vector>

namespace pfl {

// A check's flow, to load into the plot.
struct VerificationShow {
    Scene scene;
    std::string title;
};

class VerificationWindow {
public:
    VerificationWindow();
    ~VerificationWindow();
    VerificationWindow(const VerificationWindow&) = delete;
    VerificationWindow& operator=(const VerificationWindow&) = delete;

    // Shows the window, and runs the checks unless a complete run is already on show.
    void open();

    // Draws the window while it is open. Returns the flow of a check whose "Show in plot" was pressed.
    std::optional<VerificationShow> draw();

private:
    struct Entry {
        std::optional<VerificationResult> result;
        double finishedAt = 0.0; // ImGui::GetTime() when the result came in
        bool expanded = false;
    };

    static constexpr std::size_t k_trailLength = 14;
    struct Tracer {
        Vec2 position;
        std::array<Vec2, k_trailLength> trail{}; // oldest first
        std::size_t trailCount = 0;
        double age = 0.0;
    };

    void start();
    void launch(std::size_t index, double now);
    void poll(double now);
    bool running() const { return m_current < m_entries.size(); }
    bool complete() const;
    int passedCount() const;

    void drawHeader(double now, double dt);
    void updateTracers(double dt, double halfWidth, double halfHeight);
    void drawTracers(ImDrawList* dl, ImVec2 centre, float unit, double halfWidth, ImU32 colour) const;
    void drawBadge(ImDrawList* dl, ImVec2 centre, float radius, double now) const;
    void drawSummary(double now);
    std::optional<VerificationShow> drawChecks(double now);
    std::optional<VerificationShow> drawDetails(std::size_t index);
    std::string report() const;
    double random();

    bool m_open = false;
    std::vector<Entry> m_entries;
    std::future<VerificationResult> m_pending;
    std::size_t m_current = 0;  // the check running; m_entries.size() when none is
    double m_checkStartedAt = 0.0;
    double m_runStartedAt = 0.0;
    double m_runFinishedAt = 0.0;
    bool m_followRunning = false; // scroll the running check into view
    double m_copiedAt = -100.0;

    // Header animation.
    double m_lastTime = -1.0;
    double m_shownProgress = 0.0; // eases towards the fraction of checks done
    double m_circulation = 0.0;   // of the badge's flow; eases up once every check has passed
    double m_trailClock = 0.0;
    std::vector<Tracer> m_tracers;
    std::uint32_t m_seed = 0x9e3779b9u;
};

} // namespace pfl
