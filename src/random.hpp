#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#include <cassert>
#include <cstdint>
#include <random>
#include "types.hpp"

extern uint32_t seed;
extern uint32_t num_threads;
extern thread_local std::mt19937 gen;
extern thread_local std::random_device rd;

// Derive a reproducible per-thread stream from --seed and thread index.
void init_thread_rng(uint32_t i = 0);

// Draw min(npos, nsamples) distinct coordinates from [0, npos), returned ascending.
// The draw is part of the reproducible output and must not change casually.
vec<uint64_t> sample_coords(uint64_t npos, uint64_t nsamples, std::mt19937& rng);

#endif
