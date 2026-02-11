#ifndef _TABLE_H
#define _TABLE_H

#include "common.hpp"
#include "rqseq.hpp"

class SDynHT
{
  friend class SFlatHT;

public:
  void make_unique();
  void sort_columns();
  void fill_table(uint32_t nrows, rseq_sptr_t rs);
  uint64_t get_nkmers() { return nkmers; }

protected:
  uint64_t nkmers = 0;
  vvec<enc_t> enc_vvec;
};

class SFlatHT
{
  friend class SDynHT;

public:
  SFlatHT(sdynht_sptr_t source);
  SFlatHT(){};
  ~SFlatHT()
  {
    inc_v.clear();
    enc_v.clear();
  }
  void save(std::ofstream& sketch_stream);
  void load(std::ifstream& sketch_stream);
  std::vector<enc_t>::const_iterator bucket_start(uint32_t rix)
  {
    if (rix) {
      return std::next(enc_v.begin(), inc_v[rix - 1]);
    } else {
      return enc_v.begin();
    }
  }
  std::vector<enc_t>::const_iterator bucket_next(uint32_t rix)
  {
    if (rix < inc_v.size()) {
      return std::next(enc_v.begin(), inc_v[rix]);
    } else {
      return enc_v.end();
    }
  }

private:
  uint32_t nrows = 0;
  uint64_t nkmers = 0;
  vec<inc_t> inc_v;
  vec<enc_t> enc_v;
};

#endif
