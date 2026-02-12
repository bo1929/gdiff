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
  void fill_table(uint32_t nrows, rseq_sptr_t rs);
  void make_unique();
  void sort_columns();
  uint64_t get_nkmers();

protected:
  uint64_t nkmers = 0;
  vvec<enc_t> enc_vvec;
};

class SFHM
{
  friend class SDHM;

public:
  SFHM(sdhm_sptr_t source);
  SFHM() {};
  ~SFHM();
  void save(std::ofstream& sketch_stream);
  void load(std::ifstream& sketch_stream);
  std::vector<enc_t>::const_iterator bucket_start(uint32_t rix);
  std::vector<enc_t>::const_iterator bucket_next(uint32_t rix);

private:
  uint32_t nrows = 0;
  uint64_t nkmers = 0;
  vec<inc_t> inc_v;
  vec<enc_t> enc_v;
};

#endif
