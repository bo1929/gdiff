#ifndef _SKETCH_HPP
#define _SKETCH_HPP

#include "types.hpp"
#include "msg.hpp"
#include "hm.hpp"
#include "lshf.hpp"
#include "CLI11.hpp"
#include <filesystem>
#include <fstream>
#include <limits>

class Sketch
{
public:
  static constexpr uint32_t INVALID_BIX = std::numeric_limits<uint32_t>::max();

  Sketch(std::filesystem::path sketch_path);
  void load_from_offset(std::ifstream& stream, uint64_t offset);
  static void seek_past(std::ifstream& stream);
  // A k-mer is in the sketch iff its LSH value is below nrows
  uint32_t validate_bucket_ix(uint32_t bix) const noexcept { return bix < nrows ? bix : INVALID_BIX; }
  void prefetch_bucket_inc(uint32_t bix) const noexcept;
  void prefetch_bucket_enc(uint32_t bix) const noexcept;
  bool scan_bucket(uint32_t bix, enc_t enc_lr, uint32_t& hdist_min) const noexcept;
  void canonicalize();
  str get_rname() { return rname; }
  double get_rho() const;
  bool is_canonical() const { return canonical; }
  uint64_t get_timestamp() const { return timestamp; }
  sfhm_sptr_t get_sfhm_sptr() const { return sfhm; }
  lshf_sptr_t get_lshf_sptr() const { return lshf; }

private:
  str rname;
  uint8_t k;
  uint8_t w;
  uint8_t h;
  double rho;
  uint32_t nrows;
  bool canonical;
  uint64_t timestamp;
  lshf_sptr_t lshf = nullptr;
  sfhm_sptr_t sfhm = nullptr;
  std::filesystem::path sketch_path;
};

// Reads the sketch-file header and returns the byte offset of each sketch.
vec<uint64_t> read_sketch_offsets(const std::filesystem::path& sketch_path);

class BaseLSH
{
public:
  void set_lshf();
  void set_nrows();
  void set_sketch_defaults()
  {
    k = 27;
    w = k + 6;
    h = 11;
    frac = 1.0;
    canonical = true;
    nrows = uint32_t(1) << (2 * h); // recomputed by set_nrows()
  }

protected:
  uint8_t k;
  uint8_t w;
  uint8_t h;
  double frac; // FracMinLSHon top of minimizers: keep if LSH(x) < frac * 2^(2h)
  bool canonical;
  uint32_t nrows;
  lshf_sptr_t lshf = nullptr;
};

class SketchSC : public BaseLSH
{
public:
  SketchSC(CLI::App& sc);
  void process();
  bool validate_configuration();
  void write_header(std::ofstream& stream, uint32_t i);
  void write_config(std::ofstream& stream, uint32_t i);

private:
  std::vector<str> paths_v;
  std::vector<double> rho_v;
  std::filesystem::path sketch_path;
};

#endif
