/*
  Profiling implementation - singleton definitions.

  These are defined in a .cpp file (not inline in the header) to avoid ODR
  violations when linking a static library into a shared library (Python module).
  With inline definitions, separate static variable instances can exist in each
  translation unit, causing profiling data to be recorded in one instance while
  being read from another.
*/
#include "netgraph/core/profiling.hpp"

#include <algorithm>  // std::sort in dump()
#include <cstdlib>

namespace netgraph::core {

bool profiling_enabled() noexcept {
    static const bool enabled = [] {
        const char* env = std::getenv("NGRAPH_CORE_PROFILE");
        return env && env[0] == '1';
    }();
    return enabled;
}

ProfilingStats& ProfilingStats::instance() {
    static auto* inst = new ProfilingStats();
    return *inst;
}

ProfilingStats::LocalStats::LocalStats(ProfilingStats& owner) : owner_(owner) {
    owner_.register_local(this);
}

ProfilingStats::LocalStats::~LocalStats() {
    owner_.finalize_local(this);
}

void ProfilingStats::LocalStats::record(const char* name, double micros) {
    std::lock_guard<std::mutex> lock(mutex_);
    auto& e = stats_[name];
    e.total_us += micros;
    e.count++;
    if (micros < e.min_us) e.min_us = micros;
    if (micros > e.max_us) e.max_us = micros;
}

void ProfilingStats::LocalStats::append_to(std::unordered_map<std::string, Entry>& out) const {
    std::lock_guard<std::mutex> lock(mutex_);
    for (const auto& [name, entry] : stats_) {
        merge_entry(out[name], entry);
    }
}

void ProfilingStats::LocalStats::clear() {
    std::lock_guard<std::mutex> lock(mutex_);
    stats_.clear();
}

ProfilingStats::LocalStats& ProfilingStats::local_stats() {
    static thread_local LocalStats local(*this);
    return local;
}

void ProfilingStats::register_local(LocalStats* local) {
    std::lock_guard<std::mutex> lock(mutex_);
    active_locals_.push_back(local);
}

void ProfilingStats::finalize_local(LocalStats* local) {
    std::lock_guard<std::mutex> lock(mutex_);
    local->append_to(finished_stats_);
    auto it = std::find(active_locals_.begin(), active_locals_.end(), local);
    if (it != active_locals_.end()) {
        active_locals_.erase(it);
    }
}

void ProfilingStats::merge_entry(Entry& dst, const Entry& src) {
    if (src.count == 0) {
        return;
    }
    dst.total_us += src.total_us;
    dst.count += src.count;
    dst.min_us = std::min(dst.min_us, src.min_us);
    dst.max_us = std::max(dst.max_us, src.max_us);
}

std::unordered_map<std::string, ProfilingStats::Entry> ProfilingStats::snapshot_unlocked() const {
    auto snapshot = finished_stats_;
    for (const LocalStats* local : active_locals_) {
        local->append_to(snapshot);
    }
    return snapshot;
}

void ProfilingStats::record(const char* name, double micros) {
    local_stats().record(name, micros);
}

void ProfilingStats::dump(std::ostream& out) const {
    std::lock_guard<std::mutex> lock(mutex_);
    const auto snapshot = snapshot_unlocked();
    std::vector<std::pair<std::string, Entry>> ordered(snapshot.begin(), snapshot.end());
    std::sort(
        ordered.begin(),
        ordered.end(),
        [](const auto& left, const auto& right) { return left.first < right.first; }
    );

    out << "\n=== NetGraph-Core Profiling Stats ===\n";
    for (const auto& [name, entry] : ordered) {
        const double avg = entry.count > 0 ? entry.total_us / static_cast<double>(entry.count) : 0.0;
        out << name << ":\n"
            << "  calls: " << entry.count << "\n"
            << "  total: " << entry.total_us / 1000.0 << " ms\n"
            << "  avg:   " << avg << " us\n"
            << "  min:   " << entry.min_us << " us\n"
            << "  max:   " << entry.max_us << " us\n";
    }
    out << "=====================================\n";
}

void ProfilingStats::reset() {
    std::lock_guard<std::mutex> lock(mutex_);
    finished_stats_.clear();
    for (LocalStats* local : active_locals_) {
        local->clear();
    }
}

} // namespace netgraph::core
