#ifndef _GDIFF_HPP
#define _GDIFF_HPP

#include <chrono>
#include <ctime>
#include <filesystem>
#include "msg.hpp"
#include "types.hpp"
#include "random.hpp"
#include "map.hpp"
#include "dist.hpp"
#include "detect.hpp"
#include "sketch.hpp"
#include "CLI11.hpp"

#define VERSION "v0.1.0"
#define PRINT_VERSION std::cerr << "gdiff version: " << VERSION << std::endl;

extern uint32_t num_threads;

class MergeSC
{
public:
  MergeSC(CLI::App& sc);
  void merge();

private:
  std::filesystem::path sketch_path;
  vec<str> paths_v;
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
