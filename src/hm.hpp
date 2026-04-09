#ifndef _HM_HPP
#define _HM_HPP

#include <fstream>
#include "msg.hpp"
#include "types.hpp"
#include "rqseq.hpp"

class SDHM
{
  friend class SFHM;

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
  void save(std::ofstream& sketch_stream);
  void load(std::ifstream& sketch_stream);
  std::vector<enc_t>::const_iterator bucket_iter_start(uint32_t rix);
  std::vector<enc_t>::const_iterator bucket_iter_next(uint32_t rix);
  const enc_t* bucket_ptr_start(uint32_t rix) const noexcept;
  const enc_t* bucket_ptr_next(uint32_t rix) const noexcept;
  // Prefetch the inc_v cache line that holds the bucket boundaries for rix.
  void prefetch_inc(uint32_t rix) const noexcept;
  // Requires inc_v[rix] to already be in cache (call after prefetch_inc has resolved).
  void prefetch_enc(uint32_t rix) const noexcept;

private:
  uint32_t nrows = 0;
  uint64_t nkmers = 0;
  vec<inc_t> inc_v;
  vec<enc_t> enc_v;
};

#endif
