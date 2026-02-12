#ifndef _RANDOM_HPP
#define _RANDOM_HPP

#include <random>

extern uint32_t seed;
extern thread_local std::mt19937 gen;
extern thread_local std::random_device rd;

#endif
