#ifndef _ROLL_HPP
#define _ROLL_HPP

#include "CLI11.hpp"
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>

// Roll a fixed-length window over the query sequences.
// Report the MLE distance at every step.
// No background model, no derivatives and no significance testing.
class RollSC
{
public:
  struct params
  {
    uint64_t tau = 0;      // window length in k-mers (-l, required)
    uint64_t step = 1;     // stride between window starts in k-mers (-s, defaults to tau)
    uint32_t hdist_th = 3; // max Hamming distance for a k-mer to match (--hdist-th)
  };

  explicit RollSC(CLI::App& sc);
  bool validate_configuration();
  void roll();

private:
  std::filesystem::path query_path;
  std::filesystem::path sketch_path;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  params params;
};

#endif
