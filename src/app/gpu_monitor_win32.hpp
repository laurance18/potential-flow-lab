// gpu_monitor_win32.hpp: GPU load and video memory, read from the Windows performance counters
// ("GPU Engine" and "GPU Process Memory", the ones Task Manager shows) on a background thread
// about once a second, so the main loop never waits for them. Works with any vendor's driver
// (WDDM 2.0+, Windows 10 1709 and newer); on older systems the readings stay empty.
#pragma once

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <optional>
#include <thread>

namespace pfl {

struct GpuReading {
    double totalPercent = 0.0;   // the busiest engine (3D, copy, video, ...) summed over all processes
    double processPercent = 0.0; // the same engine type, this process only
    double processMemoryMB = 0.0; // dedicated video memory this process uses
};

class GpuMonitor {
public:
    GpuMonitor();
    ~GpuMonitor();
    GpuMonitor(const GpuMonitor&) = delete;
    GpuMonitor& operator=(const GpuMonitor&) = delete;

    // The latest reading; empty until the first one, or when the counters are unavailable.
    std::optional<GpuReading> latest() const;

private:
    void run();

    mutable std::mutex m_mutex;
    std::condition_variable m_wake;
    bool m_stop = false;
    std::optional<GpuReading> m_latest;
    std::thread m_thread;
};

} // namespace pfl
