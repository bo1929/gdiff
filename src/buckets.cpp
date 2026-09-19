#include "buckets.hpp"
#include <algorithm>
#include <cstring>

void Buckets::bind() noexcept
{
  bitmap = bitmap_v.empty() ? nullptr : bitmap_v.data();
  brank = brank_v.empty() ? nullptr : brank_v.data();
  start = start_v.empty() ? nullptr : start_v.data();
  enc = enc_v.empty() ? nullptr : enc_v.data();
}

void Buckets::swap(Buckets& other) noexcept
{
  using std::swap;
  swap(nrows, other.nrows);
  swap(nkmers, other.nkmers);
  swap(nnonempty, other.nnonempty);
  swap(nblocks, other.nblocks);
  swap(bitmap, other.bitmap);
  swap(brank, other.brank);
  swap(start, other.start);
  swap(enc, other.enc);
  swap(bitmap_v, other.bitmap_v);
  swap(brank_v, other.brank_v);
  swap(start_v, other.start_v);
  swap(enc_v, other.enc_v);
}

void Buckets::build(uint32_t nrows, vec<uint64_t>&& keys_v)
{
  set_nrows(nrows);
  nblocks = (static_cast<uint64_t>(nrows) + 63) / 64;

  std::sort(keys_v.begin(), keys_v.end());
  keys_v.erase(std::unique(keys_v.begin(), keys_v.end()), keys_v.end());
  nkmers = keys_v.size();

  bitmap_v.assign(nblocks, 0);
  brank_v.assign(nblocks, 0);
  enc_v.resize(nkmers);
  // One offset per nonempty bucket plus a terminator; runs of equal key_bix.
  start_v.clear();
  start_v.reserve(nkmers ? 1024 : 1);
  for (uint64_t i = 0; i < nkmers;) {
    const uint32_t bix = key_bix(keys_v[i]);
    bitmap_v[bix >> 6] |= uint64_t(1) << (bix & 63);
    start_v.push_back(static_cast<uint32_t>(i));
    do {
      enc_v[i] = key_enc(keys_v[i]);
      ++i;
    } while (i < nkmers && key_bix(keys_v[i]) == bix);
  }
  nnonempty = start_v.size();
  start_v.push_back(static_cast<uint32_t>(nkmers));
  start_v.shrink_to_fit();

  uint32_t acc = 0;
  for (uint64_t b = 0; b < nblocks; ++b) {
    brank_v[b] = acc;
    acc += static_cast<uint32_t>(__builtin_popcountll(bitmap_v[b]));
  }
  bind();
  vec<uint64_t>().swap(keys_v);
}

uint64_t Buckets::get_byte_size(uint64_t nkmers, uint64_t nnonempty, uint32_t nrows)
{
  const uint64_t nblocks = (static_cast<uint64_t>(nrows) + 63) / 64;
  return 24 + nblocks * 8 + round_up_word(nblocks * 4) + round_up_word((nnonempty + 1) * 4) + round_up_word(nkmers * 4);
}

void Buckets::save(std::ostream& os) const
{
  pad_to_word(os);
  const uint64_t nnn[3] = {nkmers, nnonempty, nblocks};
  os.write(reinterpret_cast<const char*>(nnn), sizeof(nnn));
  write_array(os, bitmap_v.data(), bitmap_v.size());
  write_array(os, brank_v.data(), brank_v.size());
  write_array(os, start_v.data(), start_v.size());
  write_array(os, enc_v.data(), enc_v.size());
}

void Buckets::view(const char*& p, const char* end, uint32_t nrows)
{
  if (end - p < 24) error_exit("Truncated sketch: bucket section header");
  uint64_t nnn[3];
  std::memcpy(nnn, p, sizeof(nnn));
  p += sizeof(nnn);
  nkmers = nnn[0];
  nnonempty = nnn[1];
  nblocks = nnn[2];
  set_nrows(nrows);
  if (nblocks != (static_cast<uint64_t>(nrows) + 63) / 64) {
    error_exit("Sketch bucket section does not match the configured nrows");
  }
  bitmap = view_array<uint64_t>(p, end, nblocks, "bucket bitmap");
  brank = view_array<uint32_t>(p, end, nblocks, "bucket ranks");
  start = view_array<uint32_t>(p, end, nnonempty + 1, "bucket offsets");
  enc = view_array<enc_t>(p, end, nkmers, "bucket entries");
  bitmap_v.clear();
  brank_v.clear();
  start_v.clear();
  enc_v.clear();
}

const char* Buckets::skip(const char* p, const char* end)
{
  if (end - p < 24) error_exit("Truncated sketch: bucket section header");
  uint64_t nnn[3];
  std::memcpy(nnn, p, sizeof(nnn));
  return p + Buckets::get_byte_size(nnn[0], nnn[1], static_cast<uint32_t>(nnn[2] * 64));
}
