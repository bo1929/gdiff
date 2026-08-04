#ifndef _GDIFF_HPP
#define _GDIFF_HPP

#include <cmath>
#include <limits>
#include <atomic>
#include <chrono>
#include <ctime>
#include <mutex>
#include <thread>
#include "msg.hpp"
#include "common.hpp"
#include "types.hpp"
#include "random.hpp"
#include "lshf.hpp"
#include "rqseq.hpp"
#include "map.hpp"
#include "dist.hpp"
#include "detect.hpp"
#include "sketch.hpp"
#include "hm.hpp"
#include "CLI11.hpp"

#define VERSION "v0.0.0"
#define PRINT_VERSION std::cerr << "??? version: " << VERSION << std::endl;
#define STRSTREAM_PRECISION 4

extern uint32_t num_threads;
extern str invocation;

class MergeSC
{
public:
  MergeSC(CLI::App& sc);
  void merge();

private:
  std::filesystem::path output_path;
  std::vector<str> sketch_paths;
};

class InfoSC
{
public:
  InfoSC(CLI::App& sc);
  void info();

private:
  std::filesystem::path sketch_path;
};

#endif
