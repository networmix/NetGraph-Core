/* FIFO admission for one cache. A call owns one shared (frozen) or exclusive (building)
   lease; lifecycle operations always request exclusive admission. Workers borrow the call's
   lease instead of enrolling again. No Python code belongs in this machinery. */
#pragma once

#include <condition_variable>
#include <cstddef>
#include <mutex>

#include "netgraph/core/detail/concurrency.hpp"

namespace netgraph::core::detail {

class CacheGate {
  struct Waiter {
    Waiter* next {nullptr};
  };

public:
  CacheGate() = default;
  CacheGate(const CacheGate&) = delete;
  CacheGate& operator=(const CacheGate&) = delete;

  class Lease {
  public:
    explicit Lease(CacheGate& gate, bool lifecycle = false) : gate_(gate) {
      Waiter waiter;
      std::unique_lock lock(gate_.mutex_);
      if (gate_.tail_ != nullptr) gate_.tail_->next = &waiter;
      else gate_.head_ = &waiter;
      gate_.tail_ = &waiter;
      try {
        concurrent_event(ConcurrentStage::CacheQueued, &gate_);
        gate_.ready_.wait(lock, [&] {
          if (gate_.head_ != &waiter || gate_.writer_) return false;
          // The mode cannot change while a reader is active. A writer's mode change
          // happens before it releases this same admission mutex in ~Lease().
          return (!lifecycle && gate_.frozen_) || gate_.readers_ == 0;
        });
      } catch (...) {
        // The waiting node belongs to this stack frame; it must never remain queued
        // when wait throws. Removing it admits the next eligible prefix.
        Waiter** link = &gate_.head_;
        Waiter* previous = nullptr;
        while (*link != &waiter) {
          previous = *link;
          link = &(*link)->next;
        }
        *link = waiter.next;
        if (gate_.tail_ == &waiter) gate_.tail_ = previous;
        gate_.ready_.notify_all();
        throw;
      }
      gate_.head_ = waiter.next;
      if (gate_.head_ == nullptr) gate_.tail_ = nullptr;
      exclusive_ = lifecycle || !gate_.frozen_;
      if (exclusive_) gate_.writer_ = true;
      else ++gate_.readers_;
      // Consecutive frozen readers may enter together; a queued lifecycle request
      // stops that prefix and cannot be overtaken by any later reader or builder.
      gate_.ready_.notify_all();
      lock.unlock();
      try {
        concurrent_event(ConcurrentStage::CacheAcquired, &gate_);
      } catch (...) {
        release();
        throw;
      }
    }
    Lease(const Lease&) = delete;
    Lease& operator=(const Lease&) = delete;
    ~Lease() noexcept { release(); }

  private:
    void release() noexcept {
      const std::lock_guard lock(gate_.mutex_);
      if (exclusive_) gate_.writer_ = false;
      else --gate_.readers_;
      gate_.ready_.notify_all();
    }

    CacheGate& gate_;
    bool exclusive_ {false};
  };

  // These mode accesses require the caller's lease. A shared lease permits only
  // frozen(); setting the mode requires exclusive lifecycle admission.
  [[nodiscard]] bool frozen() const noexcept { return frozen_; }
  void freeze(bool frozen) noexcept { frozen_ = frozen; }

#if NETGRAPH_CORE_PRIVATE_HOOKS
  [[nodiscard]] bool has_waiters() {
    const std::lock_guard lock(mutex_);
    return head_ != nullptr || tail_ != nullptr;
  }
#endif

private:
  std::mutex mutex_;
  std::condition_variable ready_;
  Waiter* head_ {nullptr};
  Waiter* tail_ {nullptr};
  std::size_t readers_ {0};
  bool writer_ {false};
  bool frozen_ {false};
};

} // namespace netgraph::core::detail
