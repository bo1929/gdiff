#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#include "types.hpp"
#include <random>

extern uint32_t seed;
extern thread_local std::mt19937 gen;
extern thread_local std::random_device rd;

// Derive a reproducible per-thread stream from --seed and tseed.
void init_thread_rng(uint32_t tseed = 0);

#endif
