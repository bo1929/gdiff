// ThreadPool: coverage of every index, and exception propagation from workers.
#include "doctest/doctest.h"
#include "msg.hpp"
#include "tpool.hpp"
#include <atomic>
#include <vector>

TEST_SUITE("ThreadPool")
{

  TEST_CASE("parallel_for visits every index once")
  {
    for (uint32_t nthreads : {1u, 4u}) {
      ThreadPool pool(nthreads);
      std::vector<std::atomic<int>> hits(1000);
      pool.parallel_for(hits.size(), 7, [&](uint64_t i) { hits[i].fetch_add(1); });
      for (const auto& h : hits)
        CHECK(h.load() == 1);
    }
  }

  TEST_CASE("an exception from a worker is rethrown on the caller and the pool stays usable")
  {
    for (uint32_t nthreads : {1u, 4u}) {
      ThreadPool pool(nthreads);
      CHECK_THROWS_AS(pool.parallel_for(1000, 1, [](uint64_t i) { if (i == 500) error_exit("boom"); }), fatal_error);

      std::atomic<uint64_t> sum{0};
      pool.parallel_for(100, 3, [&](uint64_t i) { sum.fetch_add(i); });
      CHECK(sum.load() == 4950);
    }
  }

}
