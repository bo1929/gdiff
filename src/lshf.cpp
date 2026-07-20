#include "lshf.hpp"

LSHF::LSHF(uint8_t k, uint8_t h)
  : k(k)
  , h(h)
{
  get_random_positions();
  set_lshf();
}

LSHF::LSHF(const vec<uint8_t>& ppos_v, const vec<uint8_t>& npos_v)
  : ppos_v(ppos_v)
  , npos_v(npos_v)
{
  k = npos_v.size() + ppos_v.size();
  h = ppos_v.size();
  set_lshf();
}

void LSHF::get_random_positions()
{
  uint8_t n;
  assert(h <= 16);
  assert(h < k);
  std::uniform_int_distribution<uint8_t> distrib(0, k - 1);
  for (uint8_t c = 0; c < h; c++) {
    n = distrib(gen);
    if (std::count(ppos_v.begin(), ppos_v.end(), n)) {
      c -= 1;
    } else {
      ppos_v.push_back(n);
    }
  }
  std::sort(ppos_v.begin(), ppos_v.end());
  uint8_t ix_pos = 0;
  for (uint8_t i = 0; i < k; ++i) {
    if (ix_pos < h && i == ppos_v[ix_pos])
      ix_pos++; // this position is a hash (positive) position; skip it
    else
      npos_v.push_back(i);
  }
  std::sort(ppos_v.begin(), ppos_v.end(), std::greater<uint8_t>());
}

// Precompute the delta-swap stage masks implementing compress for a fixed
// selection mask (see compress_staged in lshf.hpp).
static void gen_compress_mv(uint64_t m, arr<uint64_t, 6>& mv)
{
  uint64_t mk = ~m << 1;
  for (uint32_t i = 0; i < 6; ++i) {
    uint64_t mp = mk ^ (mk << 1);
    mp ^= mp << 2;
    mp ^= mp << 4;
    mp ^= mp << 8;
    mp ^= mp << 16;
    mp ^= mp << 32;
    mv[i] = mp & m;
    m = (m ^ mv[i]) | (mv[i] >> (1u << i));
    mk &= ~mp;
  }
}

void LSHF::set_lshf()
{
  for (int i = npos_v.size() - 1; i >= 0; --i) {
    mask_drop_lr += (0x0000000100000001ull << npos_v[i]);
    mask_drop_bp += (0x0000000000000003ull << (npos_v[i] * 2));
  }
  for (uint32_t i = 0; i < 16 - (k - h); ++i) {
    mask_drop_lr += 0x0000000000000001ull << (i + k);
  }
  for (int i = ppos_v.size() - 1; i >= 0; --i) {
    mask_hash_lr += (0x0000000100000001ull << ppos_v[i]);
    mask_hash_bp += (0x0000000000000003ull << (ppos_v[i] * 2));
  }
  for (uint32_t i = (2 * h) + 1; i < 32; ++i) {
    mask_hash_lr += (0x0000000000000001ull << i);
  }
  gen_compress_mv(mask_hash_bp, mv_hash);
  gen_compress_mv(mask_drop_lr, mv_drop_lr);
  gen_compress_mv(mask_drop_bp, mv_drop_bp);
}

uint32_t LSHF::get_npos_diff(uint32_t zc)
{
  uint32_t i = __builtin_ctz(zc);
  zc = zc >> (i + 1);
  return npos_v[i];
}

uint32_t LSHF::get_npos_accdiff(uint32_t& zc, uint32_t& i)
{
  uint32_t j = __builtin_ctz(zc) + 1;
  i += j;
  zc >>= j;
  return npos_v.rbegin()[i - 1];
}

uint64_t LSHF::inv_ppos_bp(uint32_t bp) { return deposit_bits(static_cast<uint64_t>(bp), mask_hash_bp); }

uint64_t LSHF::inv_ppos_lr(uint32_t lr) { return deposit_bits(static_cast<uint64_t>(lr), mask_drop_lr); }

char* LSHF::npos_data() { return reinterpret_cast<char*>(npos_v.data()); }

char* LSHF::ppos_data() { return reinterpret_cast<char*>(ppos_v.data()); }

uint8_t LSHF::get_k() const { return k; }

uint8_t LSHF::get_h() const { return h; }

vec<uint8_t> LSHF::get_npos() { return npos_v; }

vec<uint8_t> LSHF::get_ppos() { return ppos_v; }
