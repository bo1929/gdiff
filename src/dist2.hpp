#ifndef _DIST2_HPP
#define _DIST2_HPP

// Temporary sketch2/dist2 experiment. Delete with sketch2.*, gdiff2.cpp, and makefile gdiff2 blocks.
// Set/pairwise distances over GSK5 bundles (single set = all-vs-all within;
// two sets = cross). Every unordered pair runs both directions; the summary
// row carries both medians plus their average.
#include "CLI11.hpp"
#include "llh.hpp"
#include "sketch2.hpp"
#include "types.hpp"
#include <deque>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>

class Dist2SC
{
public:
  explicit Dist2SC(CLI::App& sc);
  bool validate_configuration();
  void dist();

private:
  // One record of a container, with lazily-loaded windows/buckets.
  struct bundle_t
  {
    std::filesystem::path path;
    sketch2_index_t idx;
  };
  struct member_t
  {
    const bundle_t* bundle = nullptr;
    uint32_t rec = 0; // record index within the bundle
    str rname;
    std::unique_ptr<Sketch2> windows; // resident (both roles act as a query)
    std::unique_ptr<Sketch2> buckets; // transient: one reference at a time
  };

  // Per-direction result shared by both output modes.
  struct dir_out_t
  {
    double d_median = nanx();
    uint64_t n_valid = 0; // number of windows with a valid distance
    uint64_t n_unmapped = 0;
    std::string samples; // rows text when --output-samples (dir column first)
  };

  // One pair (+ both directions) scheduled as a pair of directional jobs.
  struct pair_t
  {
    uint32_t a = 0; // member index in `members` (windows of ab / buckets of ba)
    uint32_t b = 0; // member index (buckets of ab / windows of ba)
    dir_out_t ab;
    dir_out_t ba;
  };

  static dir_out_t
  run_direction(const Sketch2& query, const Sketch2& reference, const str& dir_label, uint32_t hdist_th, bool output_samples);

  void resolve_source(const std::filesystem::path& path, std::vector<uint32_t>& out);
  const bundle_t* bundle_for(const std::filesystem::path& path);

  std::vector<std::unique_ptr<bundle_t>> bundles; // stable addresses
  std::vector<member_t> members;                  // deduplicated over all sources
  std::vector<uint32_t> set_a;                    // member indices (query set)
  std::vector<uint32_t> set_b;                    // member indices (reference set; empty => all-vs-all)
  std::filesystem::path set_a_path;               // positional 1: bundle or list; also the whole set in within mode
  std::filesystem::path set_b_path;               // positional 2; when absent and no --query-list/--reference-list => within
  std::filesystem::path query_list_path;
  std::filesystem::path reference_list_path;
  std::filesystem::path output_path;
  std::ofstream output_file;
  std::ostream* output_stream = &std::cout;
  bool output_samples = false;
  uint32_t hdist_th = 4;
};

#endif