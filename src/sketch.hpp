#ifndef _SKETCH_HPP
#define _SKETCH_HPP

#include "types.hpp"
#include "lshf.hpp"
#include "hm.hpp"

typedef std::vector<enc_t>::const_iterator vec_enc_it;

class Sketch
{
public:
  Sketch(std::filesystem::path sketch_path);
  void load_full_sketch();
  void make_rho_partial();
  bool check_partial(uint32_t rix);
  uint32_t search_mer(uint32_t rix, enc_t enc_lr);
  std::pair<vec_enc_it, vec_enc_it> bucket_indices(uint32_t rix);
  sfhm_sptr_t get_sfhm_sptr();
  lshf_sptr_t get_lshf();
  double get_rho();

private:
  uint8_t k;
  uint8_t w;
  uint8_t h;
  bool frac;
  uint32_t r;
  uint32_t m;
  double rho;
  uint32_t nrows;
  lshf_sptr_t lshf = nullptr;
  sfhm_sptr_t sfhm = nullptr;
  std::filesystem::path sketch_path;
};

#endif
