#ifndef _TPOOL_HPP
#define _TPOOL_HPP

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

// Synchronous parallel_for: coordinator participates; fn must not re-enter.
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

    // Hold fn so late workers cannot see a destroyed stack function.
    auto held = std::make_shared<std::function<void(uint64_t)>>(fn);
    next.store(0, std::memory_order_relaxed);
    remaining.store(n, std::memory_order_relaxed);
    acks.store(0, std::memory_order_relaxed);
    {
      std::lock_guard<std::mutex> lock(mtx);
      task_n = n;
      task_chunk = chunk;
      task_fn = held;
      ++generation;
    }
    cv.notify_all();
    drain(n, chunk, *held);
    // Wait until all indices are done and every worker has drained this generation.
    while (remaining.load(std::memory_order_acquire) != 0 || acks.load(std::memory_order_acquire) < nworkers)
      std::this_thread::yield();
  }

private:
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
      uint64_t n = 0;
      uint64_t chunk = 1;
      std::shared_ptr<std::function<void(uint64_t)>> fn;
      {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [&]() { return stop || generation != seen; });
        if (stop) return;
        seen = generation;
        n = task_n;
        chunk = task_chunk;
        fn = task_fn;
      }
      if (fn) drain(n, chunk, *fn);
      acks.fetch_add(1, std::memory_order_acq_rel);
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
  std::shared_ptr<std::function<void(uint64_t)>> task_fn;
  std::atomic<uint64_t> next{0};
  std::atomic<uint64_t> remaining{0};
  std::atomic<uint64_t> acks{0};
};

#endif
