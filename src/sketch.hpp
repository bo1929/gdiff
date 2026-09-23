#ifndef _SKETCH_HPP
#define _SKETCH_HPP

#include "CLI11.hpp"
#include "buckets.hpp"
#include "container.hpp"
#include "lshf.hpp"
#include "msg.hpp"
#include "types.hpp"
#include "windows.hpp"
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>

// A read-only view of one sketch, keeping the container's mapping alive.
class Sketch
{
public:
  Sketch() = default;
  // Built from a FASTA rather than read from a container.
  Sketch(const sketch_config_t& cfg, str rname, Buckets&& buckets, double card, double rho);

  [[nodiscard]] const str& get_rname() const noexcept { return rname; }
  [[nodiscard]] double get_rho() const noexcept { return rho; }
  [[nodiscard]] uint64_t get_nkmers() const noexcept { return nkmers; }
  [[nodiscard]] double get_card() const noexcept { return card; }
  [[nodiscard]] uint64_t get_timestamp() const noexcept { return timestamp; }
  [[nodiscard]] uint64_t get_ntotal_bp() const noexcept { return ntotal_bp; }
  [[nodiscard]] uint64_t get_nvalid_bp() const noexcept { return nvalid_bp; }
  [[nodiscard]] const sketch_config_t& get_config() const noexcept { return cfg; }
  [[nodiscard]] bool is_canonical() const noexcept { return cfg.canonical; }
  [[nodiscard]] uint8_t get_k() const noexcept { return cfg.k; }
  [[nodiscard]] uint8_t get_h() const noexcept { return cfg.h; }
  [[nodiscard]] uint32_t get_nrows() const noexcept { return cfg.nrows; }
  [[nodiscard]] const lshf_sptr_t& get_lshf_sptr() const noexcept { return lshf; }
  [[nodiscard]] const Buckets& get_buckets() const noexcept { return buckets; }
  [[nodiscard]] const window_sample_t& get_windows() const noexcept { return windows; }
  [[nodiscard]] bool has_buckets() const noexcept { return loaded_buckets; }
  [[nodiscard]] bool has_windows() const noexcept { return loaded_windows; }

  // Rebuild the buckets over canonical k-mers; a no-op when already canonical.
  void canonicalize();

private:
  str rname;
  double rho = 1.0;
  uint64_t nkmers = 0;
  double card = 0.0;
  uint64_t timestamp = 0;
  uint64_t ntotal_bp = 0;
  uint64_t nvalid_bp = 0;
  bool loaded_buckets = false;
  bool loaded_windows = false;
  sketch_config_t cfg;
  lshf_sptr_t lshf;
  Buckets buckets;
  window_sample_t windows;
  std::shared_ptr<const FileMap> map; // keeps the views above alive

  friend class Container;
};

// One path per line, '#' lines skipped; "name<TAB>path" names the sketch.
void read_path_list(const std::filesystem::path& list_path, vec<str>& paths, vec<str>& names);

// True when two sketches may be compared or merged: everything but keep_seq must match.
bool compatible_configs(const sketch_config_t& a, const sketch_config_t& b);

class BaseLSH
{
public:
  void set_lshf();
  void set_nrows();
  void set_default_params()
  {
    k = 23;
    w = k;
    h = 9;
    frac = 0.5;
    canonical = true;
    nrows = uint32_t(1) << (2 * h); // recomputed by set_nrows()
  }

protected:
  uint8_t k;
  uint8_t w;
  uint8_t h;
  double frac;
  bool canonical;
  uint32_t nrows;
  lshf_sptr_t lshf = nullptr;
};

class SketchSC : public BaseLSH
{
public:
  struct params
  {
    uint64_t tau = 333;
    uint64_t sample_size = 1000;
    bool keep_seq = false;
  };

  explicit SketchSC(CLI::App& sc);
  void process();
  bool validate_configuration();

private:
  void write_file_header(std::ostream& os, uint64_t nsketches);
  scentry write_sketch(str& bytes,
                       uint64_t timestamp,
                       uint64_t ntotal_bp,
                       uint64_t nvalid_bp,
                       const Sketch& built,
                       const window_sample_t& sample);
  void write_windows(std::ostream& os, const window_sample_t& sample);
  window_sample_t sample_windows(const str& input_path, uint64_t& ntotal_bp, uint64_t& nvalid_bp);
  sketch_config_t make_config(uint64_t timestamp) const;

  vec<str> paths_v;
  vec<str> rnames_v; // optional per-input names from --input-list
  std::filesystem::path input_list_path;
  std::filesystem::path sketch_path;
  params params;
};

// Build an in-memory sketch from a FASTA/FASTQ: owned buckets plus the HLL cardinality
// estimate, with rho = retained / cardinality clamped to (0, 1].
Sketch build_sketch(const str& input_path, const sketch_config_t& cfg, str rname);

#endif
