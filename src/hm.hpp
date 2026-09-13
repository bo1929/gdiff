#ifndef _HM_HPP
#define _HM_HPP

#include "msg.hpp"
#include "types.hpp"
#include <ostream>
#include <utility>

inline uint64_t pack_key(uint32_t bix, enc_t enc) noexcept
{
  return (static_cast<uint64_t>(bix) << 32) | static_cast<uint64_t>(enc);
}

inline uint32_t key_bix(uint64_t key) noexcept { return static_cast<uint32_t>(key >> 32); }

inline enc_t key_enc(uint64_t key) noexcept { return static_cast<enc_t>(key & 0xffffffffull); }

class Buckets
{
public:
  Buckets() = default;
  Buckets(const Buckets&) = delete;
  Buckets& operator=(const Buckets&) = delete;
  Buckets(Buckets&& other) noexcept { steal(std::move(other)); }
  Buckets& operator=(Buckets&& other) noexcept
  {
    if (this != &other) steal(std::move(other));
    return *this;
  }

  void build(uint32_t nrows, vec<uint64_t>&& keys);
  void save(std::ostream& os) const;
  void view(const char*& p, const char* end, uint32_t nrows);
  static const char* skip(const char* p, const char* end);
  static uint64_t byte_size(uint64_t nkmers, uint64_t nnonempty, uint32_t nrows);

  // Entry range of bucket bix; false when out of range or empty.
  bool range(uint32_t bix, const enc_t*& beg, const enc_t*& end) const noexcept
  {
    if (bix >= nrows) return false;
    const uint64_t word = bitmap[bix >> 6];
    const uint64_t bit = uint64_t(1) << (bix & 63);
    if (!(word & bit)) return false;
    const uint32_t rank =
      blockrank[bix >> 6] + static_cast<uint32_t>(__builtin_popcountll(word & (bit - 1)));
    beg = enc + start[rank];
    end = enc + start[rank + 1];
    return true;
  }

  void prefetch(uint32_t bix) const noexcept
  {
    if (bix >= nrows) return;
    __builtin_prefetch(&bitmap[bix >> 6], 0, 3);
    __builtin_prefetch(&blockrank[bix >> 6], 0, 3);
  }

  bool nonempty(uint32_t bix) const noexcept
  {
    if (bix >= nrows) return false;
    return (bitmap[bix >> 6] >> (bix & 63)) & 1ull;
  }

  [[nodiscard]] uint64_t get_nkmers() const noexcept { return nkmers; }
  [[nodiscard]] uint64_t get_nnonempty() const noexcept { return nnonempty; }
  [[nodiscard]] uint32_t get_nrows() const noexcept { return nrows; }
  [[nodiscard]] bool empty() const noexcept { return bitmap == nullptr; }

private:
  // Point the accessors at either the owned vectors or the borrowed buffer.
  void bind_owned() noexcept;
  void steal(Buckets&& other) noexcept;

  uint32_t nrows = 0;
  uint64_t nkmers = 0;
  uint64_t nnonempty = 0;
  uint64_t nblocks = 0;
  const uint64_t* bitmap = nullptr;
  const uint32_t* blockrank = nullptr;
  const uint32_t* start = nullptr;
  const enc_t* enc = nullptr;
  // Populated by build(); empty for a borrowed view.
  vec<uint64_t> bitmap_v;
  vec<uint32_t> blockrank_v;
  vec<uint32_t> start_v;
  vec<enc_t> enc_v;
};

// Pad a stream to the next 8-byte boundary so views stay naturally aligned.
void pad_to_8(std::ostream& os);

#endif
