#include "random.hpp"

uint32_t seed = 0;
uint32_t num_threads = 1;
std::string invocation = "";
thread_local std::mt19937 gen;
thread_local std::random_device rd;
