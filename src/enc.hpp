#ifndef _ENC_HPP
#define _ENC_HPP

#include <cstdint>
#include <array>

extern const std::array<uint8_t, 128> SEQ_NT4_TABLE;
extern const std::array<uint64_t, 4> NT4_LR_TABLE;
extern const std::array<uint64_t, 4> NT4_BP_TABLE;

inline void compute_encoding(const char* begin, const char* end, uint64_t& enc_lr, uint64_t& enc_bp)
{
  enc_lr = 0;
  enc_bp = 0;
  for (; begin < end; ++begin) {
    enc_lr <<= 1;
    enc_bp <<= 2;
    enc_bp += NT4_BP_TABLE[SEQ_NT4_TABLE[static_cast<uint8_t>(*begin)]];
    enc_lr += NT4_LR_TABLE[SEQ_NT4_TABLE[static_cast<uint8_t>(*begin)]];
  }
}

inline void update_encoding(const char* position, uint64_t& enc_lr, uint64_t& enc_bp)
{
  enc_lr <<= 1;
  enc_bp <<= 2;
  enc_lr &= 0xFFFFFFFEFFFFFFFE;
  enc_bp += NT4_BP_TABLE[SEQ_NT4_TABLE[static_cast<uint8_t>(*position)]];
  enc_lr += NT4_LR_TABLE[SEQ_NT4_TABLE[static_cast<uint8_t>(*position)]];
}

inline uint64_t revcomp_bp64(uint64_t x, uint8_t k)
{
  uint64_t res = ~x;
  res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
  res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
  res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
  res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
  res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
  return (res >> (2 * (32 - k)));
}

inline uint64_t rmoddp_bp64(uint64_t x)
{
  x = x & 0x5555555555555555;
  x = (x | (x >> 1)) & 0x3333333333333333;
  x = (x | (x >> 2)) & 0x0f0f0f0f0f0f0f0f;
  x = (x | (x >> 4)) & 0x00ff00ff00ff00ff;
  x = (x | (x >> 8)) & 0x0000ffff0000ffff;
  x = (x | (x >> 16)) & 0x00000000ffffffff;
  return x;
}

inline uint64_t bp64_to_lr64(uint64_t x) { return (rmoddp_bp64(x >> 1) << 32) | rmoddp_bp64(x); }

inline uint32_t hdist_lr64(uint64_t x, uint64_t y)
{
  uint64_t z = x ^ y;
  return __builtin_popcount(z | (z >> 32));
}

inline uint32_t hdist_lr32(uint32_t x, uint32_t y)
{
  uint32_t z = x ^ y;
  return __builtin_popcount((z | (z >> 16)) & 0x0000ffff);
}

inline uint32_t popcount_lr32(uint32_t z) { return __builtin_popcount((z | (z >> 16)) & 0x0000ffff); }

// Bit extraction for non-BMI2 architectures
template<typename Integral>
constexpr Integral extract_bits(Integral x, Integral mask)
{
  Integral res = 0;
  for (Integral bb = 1; mask != 0; bb += bb) {
    if (x & mask & -mask) {
      res |= bb;
    }
    mask &= (mask - 1);
  }
  return res;
}

#endif
