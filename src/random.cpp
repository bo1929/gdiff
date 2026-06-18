#include "random.hpp"

uint32_t seed = 0;
uint32_t num_threads = 1;
str invocation = "";
thread_local std::mt19937 gen;
thread_local std::random_device rd;

void init_thread_rng(const uint32_t tseed)
{
  uint32_t s = seed;
  s ^= tseed * 0x9E3779B9u;
  s += 0x85EBCA6Bu;
  gen.seed(s);
}
