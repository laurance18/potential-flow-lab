#include "app/gpu_monitor_win32.hpp"

#include "app/platform_win32.hpp"

#include <windows.h>
#include <pdh.h>
#include <pdhmsg.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace pfl {
namespace {

constexpr auto k_sampleInterval = std::chrono::milliseconds(1000);

// Every instance of a wildcard counter: name and value. Empty when the counter could not be read.
std::vector<std::pair<std::wstring, double>> readArray(PDH_HCOUNTER counter) {
    std::vector<std::pair<std::wstring, double>> values;
    DWORD bytes = 0;
    DWORD count = 0;
    if (PdhGetFormattedCounterArrayW(counter, PDH_FMT_DOUBLE | PDH_FMT_NOCAP100, &bytes, &count, nullptr) !=
        static_cast<PDH_STATUS>(PDH_MORE_DATA)) {
        return values;
    }
    std::vector<std::uint8_t> buffer(bytes);
    auto* items = reinterpret_cast<PDH_FMT_COUNTERVALUE_ITEM_W*>(buffer.data());
    if (PdhGetFormattedCounterArrayW(counter, PDH_FMT_DOUBLE | PDH_FMT_NOCAP100, &bytes, &count, items) !=
        ERROR_SUCCESS) {
        return values;
    }
    values.reserve(count);
    for (DWORD i = 0; i < count; ++i) {
        if (items[i].FmtValue.CStatus == PDH_CSTATUS_VALID_DATA || items[i].FmtValue.CStatus == PDH_CSTATUS_NEW_DATA) {
            values.emplace_back(items[i].szName, items[i].FmtValue.doubleValue);
        }
    }
    return values;
}

// "pid_1234_luid_0x..._phys_0_eng_3_engtype_3D" → "3D".
std::wstring engineType(const std::wstring& instance) {
    const std::size_t at = instance.find(L"engtype_");
    return at == std::wstring::npos ? std::wstring() : instance.substr(at + 8);
}

} // namespace

GpuMonitor::GpuMonitor() : m_thread([this] { run(); }) {}

GpuMonitor::~GpuMonitor() {
    {
        const std::lock_guard lock(m_mutex);
        m_stop = true;
    }
    m_wake.notify_all();
    m_thread.join();
}

std::optional<GpuReading> GpuMonitor::latest() const {
    const std::lock_guard lock(m_mutex);
    return m_latest;
}

void GpuMonitor::run() {
    PDH_HQUERY query = nullptr;
    if (PdhOpenQueryW(nullptr, 0, &query) != ERROR_SUCCESS) return;
    PDH_HCOUNTER engines = nullptr;
    PDH_HCOUNTER memory = nullptr;
    if (PdhAddEnglishCounterW(query, L"\\GPU Engine(*)\\Utilization Percentage", 0, &engines) != ERROR_SUCCESS) {
        logMessage("GPU performance counters not available; the GPU readout stays empty.");
        PdhCloseQuery(query);
        return;
    }
    if (PdhAddEnglishCounterW(query, L"\\GPU Process Memory(*)\\Dedicated Usage", 0, &memory) != ERROR_SUCCESS) {
        memory = nullptr;
    }
    const std::wstring self = L"pid_" + std::to_wstring(GetCurrentProcessId()) + L"_";
    PdhCollectQueryData(query); // utilisation is a rate: it needs two samples

    std::unique_lock lock(m_mutex);
    while (!m_wake.wait_for(lock, k_sampleInterval, [this] { return m_stop; })) {
        lock.unlock();
        GpuReading reading;
        bool valid = false;
        if (PdhCollectQueryData(query) == ERROR_SUCCESS) {
            // Task Manager's figure: per engine type, the sum over every process; then the busiest type.
            std::map<std::wstring, double> total;
            std::map<std::wstring, double> mine;
            for (const auto& [name, value] : readArray(engines)) {
                const std::wstring type = engineType(name);
                total[type] += value;
                if (name.starts_with(self)) mine[type] += value;
                valid = true;
            }
            for (const auto& [type, value] : total) reading.totalPercent = std::max(reading.totalPercent, value);
            for (const auto& [type, value] : mine) reading.processPercent = std::max(reading.processPercent, value);
            reading.totalPercent = std::min(reading.totalPercent, 100.0);
            reading.processPercent = std::min(reading.processPercent, 100.0);
            if (memory != nullptr) {
                for (const auto& [name, value] : readArray(memory)) {
                    if (name.starts_with(self)) reading.processMemoryMB += value / (1024.0 * 1024.0);
                }
            }
        }
        lock.lock();
        if (valid) m_latest = reading;
    }
    lock.unlock();
    PdhCloseQuery(query);
}

} // namespace pfl
