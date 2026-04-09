#ifndef _EXTHASH_HPP
#define _EXTHASH_HPP

#include <string>
#include <cstdint>
#include "MurmurHash3.hpp"

constexpr uint32_t MURMUR_SEED_0 = 0x5bd1e995;
constexpr uint32_t MURMUR_SEED_1 = 0x1b873593;

inline uint32_t ghhp(const str& str)
{
  uint32_t b = 378551;
  uint32_t a = 63689;
  uint32_t h = 0;
  for (char c : str) {
    h = (h * a) + c;
    a = a * b;
  }
  return (h & 0x7FFFFFFF);
}

inline uint32_t xhur32(uint32_t h)
{
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

inline uint64_t xhur64(uint64_t h)
{
  h ^= (h >> 33);
  h *= 0xff51afd7ed558ccdL;
  h ^= (h >> 33);
  h *= 0xc4ceb9fe1a85ec53L;
  h ^= (h >> 33);
  return h;
}

#endif
