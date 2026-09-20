#ifndef _GDIFF_HPP
#define _GDIFF_HPP

#include <filesystem>
#include "types.hpp"
#include "CLI11.hpp"

inline constexpr const char* gdiff_version = "v0.2.0-rc";

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
