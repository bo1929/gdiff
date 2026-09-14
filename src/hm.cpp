#include "hm.hpp"
#include <algorithm>
#include <cstring>

namespace {

  constexpr uint64_t round_up_8(uint64_t x) noexcept { return (x + 7) & ~uint64_t(7); }

  template<typename T>
  void bind_section(const char*& p, const char* end, const T*& out, uint64_t n, const char* what)
  {
    const uint64_t bytes = round_up_8(n * sizeof(T));
    if (static_cast<uint64_t>(end - p) < bytes) {
      error_exit(concat_msg("Truncated sketch: not enough bytes for ", what));
    }
    out = n ? reinterpret_cast<const T*>(p) : nullptr;
    p += bytes;
  }

  template<typename T>
  void write_section(std::ostream& os, const vec<T>& v)
  {
    if (!v.empty()) {
      os.write(reinterpret_cast<const char*>(v.data()), static_cast<std::streamsize>(v.size() * sizeof(T)));
    }
    pad_to_8(os);
  }

} // namespace

void pad_to_8(std::ostream& os)
{
  static const char zeros[8] = {};
  const uint64_t pos = static_cast<uint64_t>(os.tellp());
  const uint64_t pad = round_up_8(pos) - pos;
  if (pad) os.write(zeros, static_cast<std::streamsize>(pad));
}

void Buckets::bind_owned() noexcept
{
  bitmap = bitmap_v.empty() ? nullptr : bitmap_v.data();
  blockrank = blockrank_v.empty() ? nullptr : blockrank_v.data();
  start = start_v.empty() ? nullptr : start_v.data();
  enc = enc_v.empty() ? nullptr : enc_v.data();
}

void Buckets::steal(Buckets&& other) noexcept
{
  const bool owned = !other.bitmap_v.empty();
  nrows = other.nrows;
  nkmers = other.nkmers;
  nnonempty = other.nnonempty;
  nblocks = other.nblocks;
  bitmap_v = std::move(other.bitmap_v);
  blockrank_v = std::move(other.blockrank_v);
  start_v = std::move(other.start_v);
  enc_v = std::move(other.enc_v);
  if (owned) {
    bind_owned();
  } else {
    bitmap = other.bitmap;
    blockrank = other.blockrank;
    start = other.start;
    enc = other.enc;
  }
  other.bitmap = nullptr;
  other.blockrank = nullptr;
  other.start = nullptr;
  other.enc = nullptr;
  other.nkmers = other.nnonempty = other.nblocks = 0;
}

void Buckets::build(uint32_t nrows_arg, vec<uint64_t>&& keys)
{
  nrows = nrows_arg;
  nblocks = (static_cast<uint64_t>(nrows) + 63) / 64;

  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  nkmers = keys.size();

  bitmap_v.assign(nblocks, 0);
  blockrank_v.assign(nblocks, 0);
  enc_v.resize(nkmers);
  // One offset per nonempty bucket plus a terminator; runs of equal key_bix.
  start_v.clear();
  start_v.reserve(nkmers ? 1024 : 1);
  for (uint64_t i = 0; i < nkmers;) {
    const uint32_t bix = key_bix(keys[i]);
    bitmap_v[bix >> 6] |= uint64_t(1) << (bix & 63);
    start_v.push_back(static_cast<uint32_t>(i));
    do {
      enc_v[i] = key_enc(keys[i]);
      ++i;
    } while (i < nkmers && key_bix(keys[i]) == bix);
  }
  nnonempty = start_v.size();
  start_v.push_back(static_cast<uint32_t>(nkmers));
  start_v.shrink_to_fit();

  uint32_t acc = 0;
  for (uint64_t b = 0; b < nblocks; ++b) {
    blockrank_v[b] = acc;
    acc += static_cast<uint32_t>(__builtin_popcountll(bitmap_v[b]));
  }
  bind_owned();
  vec<uint64_t>().swap(keys);
}

uint64_t Buckets::byte_size(uint64_t nkmers, uint64_t nnonempty, uint32_t nrows)
{
  const uint64_t nblocks = (static_cast<uint64_t>(nrows) + 63) / 64;
  return 24 + nblocks * 8 + round_up_8(nblocks * 4) + round_up_8((nnonempty + 1) * 4) + round_up_8(nkmers * 4);
}

void Buckets::save(std::ostream& os) const
{
  pad_to_8(os);
  const uint64_t hdr[3] = {nkmers, nnonempty, nblocks};
  os.write(reinterpret_cast<const char*>(hdr), sizeof(hdr));
  write_section(os, bitmap_v);
  write_section(os, blockrank_v);
  write_section(os, start_v);
  write_section(os, enc_v);
}

void Buckets::view(const char*& p, const char* end, uint32_t nrows_arg)
{
  if (end - p < 24) error_exit("Truncated sketch: bucket section header");
  uint64_t hdr[3];
  std::memcpy(hdr, p, sizeof(hdr));
  p += sizeof(hdr);
  nkmers = hdr[0];
  nnonempty = hdr[1];
  nblocks = hdr[2];
  nrows = nrows_arg;
  if (nblocks != (static_cast<uint64_t>(nrows) + 63) / 64) {
    error_exit("Sketch bucket section does not match the configured nrows");
  }
  bind_section(p, end, bitmap, nblocks, "bucket bitmap");
  bind_section(p, end, blockrank, nblocks, "bucket ranks");
  bind_section(p, end, start, nnonempty + 1, "bucket offsets");
  bind_section(p, end, enc, nkmers, "bucket entries");
  bitmap_v.clear();
  blockrank_v.clear();
  start_v.clear();
  enc_v.clear();
}

const char* Buckets::skip(const char* p, const char* end)
{
  if (end - p < 24) error_exit("Truncated sketch: bucket section header");
  uint64_t hdr[3];
  std::memcpy(hdr, p, sizeof(hdr));
  return p + Buckets::byte_size(hdr[0], hdr[1], static_cast<uint32_t>(hdr[2] * 64));
}
