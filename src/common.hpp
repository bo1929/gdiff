#ifndef _COMMON_HPP
#define _COMMON_HPP

#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <locale>
#include <math.h>
#include <memory>
#include <numeric>
#include <ostream>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <type_traits>
#include <utility>
#include <vector>
#include <stdio.h>

#define STRSTREAM_PRECISION 4

#define VERSION "v0.6.0"
#define PRINT_VERSION std::cerr << "krepp version: " << VERSION << std::endl;

extern uint32_t num_threads;
extern std::string invocation;

static std::string vec_to_str(const std::vector<uint8_t>& v)
{
  std::ostringstream oss;
  oss << "[";
  for (size_t i = 0; i < v.size(); ++i) {
    if (i > 0) oss << ", ";
    oss << static_cast<int>(v[i]);
  }
  oss << "]";
  return oss.str();
}

#endif
