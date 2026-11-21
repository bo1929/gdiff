#ifndef _DIST_HPP
#define _DIST_HPP

#include "CLI11.hpp"
#include "container.hpp"
#include "distance.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <optional>

class ThreadPool;

// One window's estimate.
struct mle_t
{
  double d = nanx();
  double lr_ub = nanx(); // likelihood-ratio statistic of d against the distance ceiling
  char strand = '.';     // strand that supplied d; '.' for canonical sketches
};

// One direction of a pair: the distance of every window of the query to the reference.
struct profile_t
{
  vec<mle_t> windows_v;
  double d_median = nanx();
  double d_mean = nanx();
  double d_highest = nanx();
  uint64_t n_valid = 0; // windows with a distance
};

profile_t obtain_profile(const Sketch& query, const Sketch& reference, uint32_t hdist_th);

// Reconciled summary of one pair: both directions merged rank-wise.
struct summary_t
{
  double d = nanx();       // reported estimate
  vec<double> d_v;         // distances behind `d`
  uint64_t n_na = 0;       // ranks without a distance estimate
  uint64_t n_ub = 0;       // windows whose likelihood-ratio statistic is zero
  uint64_t n_filtered = 0; // windows rejected by the likelihood-ratio filter
};

summary_t summarize_symmetric(vec<mle_t> ab_v, vec<mle_t> ba_v, double lr_th, double min_portion);

// An explicit --hdist-th overrides; otherwise the threshold depends on the reference size.
uint32_t hdist_th_for(uint64_t nvalid_bp, std::optional<uint32_t> hdist_th);

// Compare sets of sketches: every pair, both directions, reconcile data points after sorting.
class DistSC
{
public:
  struct params
  {
    bool output_samples = false;
    std::optional<uint32_t> hdist_th; // unset: chosen per reference by hdist_th_for
    double lr_th = 6.635;
    double min_portion = 0.66;
  };

  explicit DistSC(CLI::App& sc);
  bool validate_configuration();
  void dist();

private:
  // A container opened once, however many times its path is given.
  struct loaded_t
  {
    std::unique_ptr<Container> file; // stable address: entries point at it
    uint32_t first = 0;              // index of its first entry
  };

  // One sketch of one container, resolved from a positional or a list file.
  struct entry_t
  {
    const Container* file = nullptr;
    uint32_t rix = 0;
    str rname;
  };

  // An unordered pair scheduled as two one-way jobs: `ab` has a as the query, `ba` has b.
  struct pair_t
  {
    uint32_t a = 0;
    uint32_t b = 0;
    profile_t ab;
    profile_t ba;
  };

  void resolve_container(const std::filesystem::path& path, vec<uint32_t>& out);
  void resolve_list(const std::filesystem::path& path, vec<uint32_t>& out);
  // `set` is 'a' or 'b'; the set comes from a positional container or from a list file.
  void resolve_set(const std::filesystem::path& container_path,
                   const std::filesystem::path& list_path,
                   char set,
                   vec<uint32_t>& out);

  // Every unordered pair within set A, or the cross product of A and B.
  void estimate_distances();
  void require_comparable_entries() const;
  vec<Sketch> load_query_sketches(ThreadPool& pool);
  vec<pair_t> plan_pairs() const;
  void score_pairs(ThreadPool& pool, const vec<Sketch>& targets_v, vec<pair_t>& pairs_v);

  void write_header(std::ostream& os, bool canonical) const;
  void write_pair_line(std::ostream& os, const pair_t& pr) const;
  void write_samples(std::ostream& os, const pair_t& pr, const vec<Sketch>& targets_v) const;

  vec<loaded_t> loaded_v;
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
