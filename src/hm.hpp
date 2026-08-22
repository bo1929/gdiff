#ifndef _HM_HPP
#define _HM_HPP

#include <fstream>
#include "msg.hpp"
#include "types.hpp"

class SDHM
{
  friend class SFHM;
  friend class Sketch;

public:
  void fill_table(uint32_t nrows, const rseq_sptr_t& rs);
  void make_unique();
  void sort_columns();
  uint64_t get_nmers() const;

protected:
  uint64_t nkmers = 0;
  vvec<enc_t> enc_vvec;
};

class SFHM
{
  friend class SDHM;

public:
  SFHM(const sdhm_sptr_t& source);
  SFHM() = default;
  ~SFHM();
  void save(std::ostream& sketch_stream);
  void load(std::ifstream& sketch_stream);
  // Copy SFHM payload from a memory cursor (e.g. mmap). Advances p.
  void load_mem(const char*& p, const char* end);
  // Zero-copy view into an external buffer (e.g. mmap) at cursor p. Advances p.
  // The caller must keep the underlying memory alive for the SFHM's lifetime
  // (e.g., hold a shared_ptr<Sketch2::MappedFile>). No vectors are allocated.
  void view_mem(const char*& p, const char* end);
  bool is_view() const noexcept { return view_mode_; }
  uint64_t get_nkmers() const { return nkmers; }
  uint32_t get_nrows() const { return nrows; }
  std::vector<enc_t>::const_iterator bucket_iter_start(uint32_t rix);
  std::vector<enc_t>::const_iterator bucket_iter_next(uint32_t rix);
  const enc_t* bucket_ptr_start(uint32_t rix) const noexcept;
  const enc_t* bucket_ptr_next(uint32_t rix) const noexcept;
  // Prefetch the inc_v cache line that holds the bucket boundaries for rix
  void prefetch_inc(uint32_t rix) const noexcept;
  // Requires inc_v[rix] to already be in cache (call after prefetch_inc has resolved)
  void prefetch_enc(uint32_t rix) const noexcept;
  // Mark bit i if bucket i is nonempty. bits must hold ceil(nbits/64) words.
  void fill_nonempty_bitmap(uint64_t* bits, uint32_t nbits) const noexcept;

private:
  // Inc/enc arrays are either owned (inc_v/enc_v) or a borrowed view into a
  // live external buffer (inc_ptr/enc_ptr, view_mode_==true).
  uint32_t nrows = 0;
  uint64_t nkmers = 0;
  const inc_t* inc_view = nullptr;
  const enc_t* enc_view = nullptr;
  bool view_mode_ = false;
  vec<inc_t> inc_v;
  vec<enc_t> enc_v;
};

#endif
