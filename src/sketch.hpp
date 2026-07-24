#ifndef _SKETCH_HPP
#define _SKETCH_HPP

#include "types.hpp"
#include "msg.hpp"
#include "hm.hpp"
#include "lshf.hpp"

class Sketch
{
public:
  static constexpr uint32_t OFF_INVALID = std::numeric_limits<uint32_t>::max();

  Sketch(std::filesystem::path sketch_path);
  void load_from_offset(std::ifstream& stream, uint64_t offset);
  static void seek_past(std::ifstream& stream);

  // Bucket offset for a hash value: a k-mer is in the sketch iff its LSH value
  // is below nrows (the keep threshold), and the offset is the hash itself.
  uint32_t partial_offset(uint32_t rix) const noexcept { return rix < nrows ? rix : OFF_INVALID; }

  void prefetch_offset_inc(uint32_t offset) const noexcept;
  void prefetch_offset_enc(uint32_t offset) const noexcept;
  bool scan_bucket(uint32_t offset, enc_t enc_lr, uint32_t& hdist_min) const noexcept;
  sfhm_sptr_t get_sfhm_sptr();
  lshf_sptr_t get_lshf();
  void canonicalize();
  double get_rho() const;
  bool is_canonical() const { return canonical; }
  str get_rid() { return rid; }
  uint64_t get_timestamp() const { return timestamp; }

private:
  uint8_t k;
  uint8_t w;
  uint8_t h;
  double rho;
  bool canonical;
  uint32_t nrows;
  lshf_sptr_t lshf = nullptr;
  sfhm_sptr_t sfhm = nullptr;
  uint64_t timestamp;
  str rid;
  std::filesystem::path sketch_path;
};

#endif
