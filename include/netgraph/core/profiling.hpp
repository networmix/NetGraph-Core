/* Simple profiling utilities for performance analysis.

Enable via environment variable: NGRAPH_CORE_PROFILE=1

When disabled (default), overhead is minimal: ~2-4 CPU cycles per scope
due to a single static bool check with branch prediction.

Usage:
    #include "netgraph/core/profiling.hpp"

    void my_function() {
        NGRAPH_PROFILE_SCOPE("my_function");
        // ... function body ...
    }

From Python:
    import netgraph_core as ngc
    # Run with NGRAPH_CORE_PROFILE=1 python script.py
    if ngc.profiling_enabled():
        ngc.profiling_dump()   # Print stats to stderr
        ngc.profiling_reset()  # Clear stats
*/
#pragma once

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <string>
#include <unordered_map>

namespace netgraph::core {

// Check once at startup, cache result. Returns true if NGRAPH_CORE_PROFILE=1.
inline bool profiling_enabled() noexcept {
    static const bool enabled = [] {
        const char* env = std::getenv("NGRAPH_CORE_PROFILE");
        return env && env[0] == '1';
    }();
    return enabled;
}

// Singleton collecting profiling statistics.
class ProfilingStats {
public:
    static ProfilingStats& instance() {
        static ProfilingStats inst;
        return inst;
    }

    void record(const char* name, double micros) {
        std::lock_guard<std::mutex> lock(mutex_);
        auto& e = stats_[name];
        e.total_us += micros;
        e.count++;
        if (micros < e.min_us) e.min_us = micros;
        if (micros > e.max_us) e.max_us = micros;
    }

    void dump() {
        std::lock_guard<std::mutex> lock(mutex_);
        std::cerr << "\n=== NetGraph-Core Profiling Stats ===\n";
        for (const auto& [name, e] : stats_) {
            double avg = e.count > 0 ? e.total_us / e.count : 0.0;
            std::cerr << name << ":\n"
                      << "  calls: " << e.count << "\n"
                      << "  total: " << e.total_us / 1000.0 << " ms\n"
                      << "  avg:   " << avg << " us\n"
                      << "  min:   " << e.min_us << " us\n"
                      << "  max:   " << e.max_us << " us\n";
        }
        std::cerr << "=====================================\n";
    }

    void reset() {
        std::lock_guard<std::mutex> lock(mutex_);
        stats_.clear();
    }

private:
    ProfilingStats() = default;

    struct Entry {
        double total_us = 0;
        int64_t count = 0;
        double min_us = 1e18;
        double max_us = 0;
    };

    std::mutex mutex_;
    std::unordered_map<std::string, Entry> stats_;
};

// RAII timer. Does nothing if name is nullptr.
class ScopedTimer {
public:
    explicit ScopedTimer(const char* name) noexcept : name_(name) {
        if (name_) {
            start_ = std::chrono::high_resolution_clock::now();
        }
    }

    ~ScopedTimer() {
        if (name_) {
            auto end = std::chrono::high_resolution_clock::now();
            double us = std::chrono::duration<double, std::micro>(end - start_).count();
            ProfilingStats::instance().record(name_, us);
        }
    }

    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;

private:
    const char* name_;
    std::chrono::high_resolution_clock::time_point start_;
};

// Concatenate tokens for unique variable names
#define NGRAPH_CONCAT_IMPL(a, b) a##b
#define NGRAPH_CONCAT(a, b) NGRAPH_CONCAT_IMPL(a, b)

// Main profiling macro. Expands to a ScopedTimer only if profiling is enabled.
// When disabled, the check is a single static bool read (~1-2 cycles).
#define NGRAPH_PROFILE_SCOPE(name) \
    ::netgraph::core::ScopedTimer NGRAPH_CONCAT(_ngraph_timer_, __LINE__)( \
        ::netgraph::core::profiling_enabled() ? (name) : nullptr)

} // namespace netgraph::core
