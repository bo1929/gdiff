#ifndef _LSHF_HPP
#define _LSHF_HPP

#include "enc.hpp"
#include "msg.hpp"
#include "types.hpp"
#include "random.hpp"
#if defined(__BMI2__)
  #include <immintrin.h>
#endif

class LSHF
{
public:
  LSHF(uint8_t k, uint8_t h, uint32_t m);
  LSHF(uint32_t m, vec<uint8_t> ppos_v, vec<uint8_t> npos_v);
  void get_random_positions();
  void set_lshf();
  uint32_t compute_hash(uint64_t enc_bp);
  uint32_t drop_ppos_lr(uint64_t enc64_lr);
  uint32_t drop_ppos_bp(uint64_t enc64_bp);
  uint32_t get_npos_accdiff(uint32_t& zc, uint32_t& i);
  uint32_t get_npos_diff(uint32_t zc);
  char* npos_data();
  char* ppos_data();
  vec<uint8_t> get_npos();
  vec<uint8_t> get_ppos();
  uint8_t get_k();
  uint8_t get_h();
  uint32_t get_m();

private:
  uint8_t k;
  uint8_t h;
  uint32_t m;
  vec<uint8_t> npos_v;
  vec<uint8_t> ppos_v;
  vec<std::pair<int8_t, int8_t>> glsh_v;
  uint64_t mask_drop_lr = 0;
  uint64_t mask_drop_bp = 0;
  uint64_t mask_hash_lr = 0;
  uint64_t mask_hash_bp = 0;
  // uint64_t mask_drop_l = 0;
  // uint64_t mask_drop_r = 0;
};

#endif
