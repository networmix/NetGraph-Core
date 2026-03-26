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

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace netgraph::core {

// Check once at startup, cache result. Returns true if NGRAPH_CORE_PROFILE=1.
// Defined in profiling.cpp to avoid ODR violations with static library linking.
bool profiling_enabled() noexcept;

// Singleton collecting profiling statistics.
class ProfilingStats {
public:
    // Defined in profiling.cpp to avoid ODR violations with static library linking.
    static ProfilingStats& instance();

    void record(const char* name, double micros);
    void dump() const { dump(std::cerr); }
    void dump(std::ostream& out) const;
    void reset();

private:
    ProfilingStats() = default;
    struct Entry {
        double total_us = 0;
        int64_t count = 0;
        double min_us = 1e18;
        double max_us = 0;
    };

    class LocalStats {
    public:
        explicit LocalStats(ProfilingStats& owner);
        ~LocalStats();

        void record(const char* name, double micros);
        void append_to(std::unordered_map<std::string, Entry>& out) const;
        void clear();

    private:
        ProfilingStats& owner_;
        mutable std::mutex mutex_;
        std::unordered_map<std::string, Entry> stats_;
    };

    LocalStats& local_stats();
    void register_local(LocalStats* local);
    void finalize_local(LocalStats* local);
    static void merge_entry(Entry& dst, const Entry& src);
    std::unordered_map<std::string, Entry> snapshot_unlocked() const;

    mutable std::mutex mutex_;
    std::vector<LocalStats*> active_locals_;
    std::unordered_map<std::string, Entry> finished_stats_;
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
