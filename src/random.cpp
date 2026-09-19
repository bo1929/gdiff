#include "random.hpp"

#include <algorithm>
#include <numeric>

uint32_t seed = 0;
uint32_t num_threads = 1;
bool verbose = false;
thread_local std::mt19937 gen;
thread_local std::random_device rd;

void init_thread_rng(uint32_t i)
{
  uint32_t s = seed;
  s ^= i * 0x9E3779B9u;
  s += 0x85EBCA6Bu;
  gen.seed(s);
}

vec<uint64_t> sample_coords(uint64_t npos, uint64_t nsamples, std::mt19937& rng)
{
  assert(npos >= 1);
  const uint64_t n = std::min(npos, nsamples);
  if (n == npos) {
    vec<uint64_t> out_v(npos);
    std::iota(out_v.begin(), out_v.end(), uint64_t(0));
    return out_v;
  }
  vec<uint64_t> out_v;
  out_v.reserve(n);
  std::uniform_int_distribution<uint64_t> pick(0, npos - 1);
  while (out_v.size() < n) {
    const uint64_t need = n - out_v.size();
    for (uint64_t i = 0; i < need; ++i)
      out_v.push_back(pick(rng));
    std::sort(out_v.begin(), out_v.end());
    out_v.erase(std::unique(out_v.begin(), out_v.end()), out_v.end());
  }
  return out_v;
}
