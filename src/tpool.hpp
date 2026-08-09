#ifndef _TPOOL_HPP
#define _TPOOL_HPP

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

// Minimal persistent thread pool for coarse-grained parallel-for sections.
//
// Implementation of parallel_for(n, chunk, fn) is synchronous:
// It invokes fn(i) exactly once for each i in [0, n) using dynamic block scheduling.
// The return happens only after all invocations have completed.
// Hence, the invocation is finished when no worker still references the task.
// The calling (coordinator) thread participates in the work.
// Not re-entrant: fn must not call parallel_for and must be driven froma coordinator thread.
class ThreadPool
{
public:
  explicit ThreadPool(uint32_t nthreads)
  {
    const uint32_t hw = std::max(1u, std::thread::hardware_concurrency());
    nthreads = std::max(1u, std::min(nthreads, hw));
    nworkers = nthreads - 1;
    workers.reserve(nworkers);
    for (uint32_t t = 0; t < nworkers; ++t)
      workers.emplace_back([this]() { worker_loop(); });
  }

  ~ThreadPool()
  {
    {
      std::lock_guard<std::mutex> lock(mtx);
      stop = true;
      ++generation;
    }
    cv.notify_all();
    for (auto& w : workers)
      w.join();
  }

  ThreadPool(const ThreadPool&) = delete;
  ThreadPool& operator=(const ThreadPool&) = delete;

  uint32_t size() const { return nworkers + 1; }

  void parallel_for(uint64_t n, uint64_t chunk, const std::function<void(uint64_t)>& fn)
  {
    if (n == 0) return;
    if (nworkers == 0) {
      for (uint64_t i = 0; i < n; ++i)
        fn(i);
      return;
    }
    chunk = std::max<uint64_t>(chunk, 1);
    next.store(0, std::memory_order_relaxed);
    remaining.store(n, std::memory_order_relaxed);
    active.store(0, std::memory_order_relaxed);
    {
      std::lock_guard<std::mutex> lock(mtx);
      task_n = n;
      task_chunk = chunk;
      task_fn = &fn;
      ++generation;
    }
    cv.notify_all();
    drain(n, chunk, fn);
    // remaining == 0 once every claimed index has been processed; active == 0
    // once no worker can still touch this task's state.
    while (remaining.load(std::memory_order_acquire) != 0 || active.load(std::memory_order_acquire) != 0)
      std::this_thread::yield();
  }

private:
  struct task_snapshot_t
  {
    uint64_t n;
    uint64_t chunk;
    const std::function<void(uint64_t)>* fn;
  };

  static void drain(uint64_t n,
                    uint64_t chunk,
                    const std::function<void(uint64_t)>& fn,
                    std::atomic<uint64_t>& next,
                    std::atomic<uint64_t>& remaining)
  {
    while (true) {
      const uint64_t base = next.fetch_add(chunk, std::memory_order_relaxed);
      if (base >= n) break;
      const uint64_t end = std::min(base + chunk, n);
      for (uint64_t i = base; i < end; ++i)
        fn(i);
      remaining.fetch_sub(end - base, std::memory_order_acq_rel);
    }
  }

  void drain(uint64_t n, uint64_t chunk, const std::function<void(uint64_t)>& fn) { drain(n, chunk, fn, next, remaining); }

  void worker_loop()
  {
    uint64_t seen = 0;
    while (true) {
      task_snapshot_t task;
      {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&]() { return stop || generation != seen; });
        if (stop) return;
        seen = generation;
        task = {task_n, task_chunk, task_fn};
        active.fetch_add(1, std::memory_order_acq_rel);
      }
      drain(task.n, task.chunk, *task.fn);
      active.fetch_sub(1, std::memory_order_acq_rel);
    }
  }

  std::vector<std::thread> workers;
  uint32_t nworkers = 0;
  std::mutex mtx;
  std::condition_variable cv;
  uint64_t generation = 0;
  bool stop = false;
  uint64_t task_n = 0;
  uint64_t task_chunk = 1;
  const std::function<void(uint64_t)>* task_fn = nullptr;
  std::atomic<uint64_t> next{0};
  std::atomic<uint64_t> remaining{0};
  std::atomic<uint64_t> active{0};
};

#endif
