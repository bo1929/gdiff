#ifndef _LSHF_HPP
#define _LSHF_HPP

#include "enc.hpp"
#include "types.hpp"
#if defined(__BMI2__)
  #include <immintrin.h>
#endif

// Branchless bit compaction via a precomputed delta-swap network.
// (Hacker's Delight "compress", generalized to 64 bits).
// Applies the same transformation as PEXT for a fixed mask.
inline uint64_t compress_mv(uint64_t x, const arr<uint64_t, 6>& mv)
{
  uint64_t t;
  t = x & mv[0];
  x = (x ^ t) | (t >> 1);
  t = x & mv[1];
  x = (x ^ t) | (t >> 2);
  t = x & mv[2];
  x = (x ^ t) | (t >> 4);
  t = x & mv[3];
  x = (x ^ t) | (t >> 8);
  t = x & mv[4];
  x = (x ^ t) | (t >> 16);
  t = x & mv[5];
  x = (x ^ t) | (t >> 32);
  return x;
}

class LSHF
{
public:
  LSHF(uint8_t k, uint8_t h);
  LSHF(const vec<uint8_t>& ppos_v, const vec<uint8_t>& npos_v);
  void get_random_positions();
  void set_lshf();
#if defined(__BMI2__)
  uint32_t compute_hash_bp(uint64_t enc64_bp) const { return static_cast<uint32_t>(_pext_u64(enc64_bp, mask_hash_bp)); }
  uint32_t drop_ppos_lr(uint64_t enc64_lr) const { return static_cast<uint32_t>(_pext_u64(enc64_lr, mask_drop_lr)); }
  uint32_t drop_ppos_bp(uint64_t enc64_bp) const { return static_cast<uint32_t>(_pext_u64(enc64_bp, mask_drop_bp)); }
#else
  uint32_t compute_hash_bp(uint64_t enc64_bp) const
  {
    return static_cast<uint32_t>(compress_mv(enc64_bp & mask_hash_bp, mv_hash_bp));
  }
  uint32_t drop_ppos_lr(uint64_t enc64_lr) const
  {
    return static_cast<uint32_t>(compress_mv(enc64_lr & mask_drop_lr, mv_drop_lr));
  }
  uint32_t drop_ppos_bp(uint64_t enc64_bp) const
  {
    return static_cast<uint32_t>(compress_mv(enc64_bp & mask_drop_bp, mv_drop_bp));
  }
#endif
  uint64_t inv_ppos_bp(uint32_t bp);
  uint64_t inv_ppos_lr(uint32_t lr);
  uint32_t get_npos_accdiff(uint32_t& zc, uint32_t& i);
  uint32_t get_npos_diff(uint32_t zc);
  char* npos_data();
  char* ppos_data();
  vec<uint8_t> get_npos_v();
  vec<uint8_t> get_ppos_v();
  uint8_t get_k() const;
  uint8_t get_h() const;

private:
  uint8_t k;
  uint8_t h;
  vec<uint8_t> ppos_v;
  vec<uint8_t> npos_v;
  uint64_t mask_drop_lr = 0;
  uint64_t mask_drop_bp = 0;
  uint64_t mask_hash_lr = 0;
  uint64_t mask_hash_bp = 0;
  arr<uint64_t, 6> mv_hash_bp{};
  arr<uint64_t, 6> mv_drop_lr{};
  arr<uint64_t, 6> mv_drop_bp{};
};

#endif
