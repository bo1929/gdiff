#ifndef _MERGE_HPP
#define _MERGE_HPP

#include "CLI11.hpp"
#include "types.hpp"
#include <filesystem>

// Concatenate containers with compatible configurations into one.
class MergeSC
{
public:
  explicit MergeSC(CLI::App& sc);
  void merge();

private:
  vec<str> in_path_v;
  std::filesystem::path out_path;
};

#endif
