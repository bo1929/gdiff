#ifndef _DIST_HPP
#define _DIST_HPP

#include "CLI11.hpp"
#include "container.hpp"
#include "estimator.hpp"
#include "sym.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>

// Compare sets of sketches: every pair, both directions, one reconciled row each.
class DistSC
{
public:
  struct params
  {
    bool output_samples = false;
    uint32_t hdist_th = 3;
    double lr_th = lr_th_default;             // TODO:?
    double min_portion = min_portion_default; // TODO:?
  };

  explicit DistSC(CLI::App& sc);
  bool validate_configuration();
  void dist();

private:
  // One sketch of one container, resolved from a positional or a list file.
  struct entry_t
  {
    const Container* file = nullptr;
    uint32_t rix = 0;
    str rname;
  };

  // An unordered pair scheduled as two directional jobs.
  struct pair_t
  {
    uint32_t a = 0;
    uint32_t b = 0;
    direction_t ab;
    direction_t ba;
  };

  // Every unordered pair within set A, or the cross product of A and B.
  void estimate_distances();
  const Container* container_for(const std::filesystem::path& path);
  void resolve_container(const std::filesystem::path& path, vec<uint32_t>& out);
  void resolve_list(const std::filesystem::path& path, vec<uint32_t>& out);
  void write_pair_row(std::ostream& os, const pair_t& pr, const str& name_a, const str& name_b);
  void emit_header(std::ostream& os) const;

  vec<std::unique_ptr<Container>> containers_v; // stable addresses
  vec<entry_t> entries_v;
  vec<uint32_t> set_a_v;
  vec<uint32_t> set_b_v;
  std::filesystem::path container_a_path;
  std::filesystem::path container_b_path;
  std::filesystem::path list_a_path;
  std::filesystem::path list_b_path;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  params params;
};

#endif
