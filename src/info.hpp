#ifndef _INFO_HPP
#define _INFO_HPP

#include "CLI11.hpp"
#include <filesystem>

// Print the configuration and per-sketch metadata of a container.
class InfoSC
{
public:
  explicit InfoSC(CLI::App& sc);
  void info();

private:
  std::filesystem::path sketch_path;
};

#endif
