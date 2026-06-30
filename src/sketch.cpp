#include "sketch.hpp"

#include <utility>

Sketch::Sketch(std::filesystem::path sketch_path)
  : sketch_path(std::move(sketch_path))
{
}

void Sketch::load_from_offset(std::ifstream& stream, uint64_t offset)
{
  if (offset > 0) {
    stream.seekg(offset);
  }

  uint64_t rid_len;
  stream.read(reinterpret_cast<char*>(&rid_len), sizeof(uint64_t));
  rid.resize(rid_len);
  stream.read(&rid[0], rid_len);
  stream.read(reinterpret_cast<char*>(&timestamp), sizeof(uint64_t));

  stream.read(reinterpret_cast<char*>(&k), sizeof(uint8_t));
  stream.read(reinterpret_cast<char*>(&w), sizeof(uint8_t));
  stream.read(reinterpret_cast<char*>(&h), sizeof(uint8_t));
  stream.read(reinterpret_cast<char*>(&m), sizeof(uint32_t));
  stream.read(reinterpret_cast<char*>(&r), sizeof(uint32_t));
  stream.read(reinterpret_cast<char*>(&frac), sizeof(bool));
  stream.read(reinterpret_cast<char*>(&canonical), sizeof(bool));
  stream.read(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));

  vec<uint8_t> ppos_v(h), npos_v(k - h);
  stream.read(reinterpret_cast<char*>(ppos_v.data()), h * sizeof(uint8_t));
  stream.read(reinterpret_cast<char*>(npos_v.data()), (k - h) * sizeof(uint8_t));

  lshf = std::make_shared<LSHF>(m, ppos_v, npos_v);

  stream.read(reinterpret_cast<char*>(&rho), sizeof(double));

  sfhm = std::make_shared<SFHM>();
  sfhm->load(stream);

  // check_fstream(stream, "Failed to read the sketch file!", sketch_path);
}

void Sketch::seek_past(std::ifstream& stream)
{
  uint64_t rid_len;
  stream.read(reinterpret_cast<char*>(&rid_len), sizeof(uint64_t));
  // Skip rid (rid_len bytes) + timestamp (8 bytes)
  stream.seekg(static_cast<std::streamoff>(rid_len) + static_cast<std::streamoff>(sizeof(uint64_t)), std::ios::cur);

  uint8_t k, w, h; // Read k, skip one byte for w and then read h
  stream.read(reinterpret_cast<char*>(&k), sizeof(uint8_t));
  stream.seekg(sizeof(uint8_t), std::ios::cur);
  stream.read(reinterpret_cast<char*>(&h), sizeof(uint8_t));
  // skip: m (4), r (4), frac (1) = 9 bytes
  stream.seekg(sizeof(uint32_t) + sizeof(uint32_t) + sizeof(bool), std::ios::cur);
  // read canonical (1 byte)
  bool canonical = true;
  stream.read(reinterpret_cast<char*>(&canonical), sizeof(bool));
  // skip: nrows (4 bytes)
  stream.seekg(sizeof(uint32_t), std::ios::cur);
  // skip: ppos_v (h bytes) + npos_v ((k-h) bytes) = k bytes total
  stream.seekg(static_cast<std::streamoff>(k), std::ios::cur);
  // skip: rho (8 bytes)
  stream.seekg(static_cast<std::streamoff>(sizeof(double)), std::ios::cur);

  uint64_t nkmers;
  stream.read(reinterpret_cast<char*>(&nkmers), sizeof(uint64_t));
  stream.seekg(static_cast<std::streamoff>(nkmers) * static_cast<std::streamoff>(sizeof(enc_t)), std::ios::cur);
  uint32_t sfhm_nrows;
  stream.read(reinterpret_cast<char*>(&sfhm_nrows), sizeof(uint32_t));
  stream.seekg(static_cast<std::streamoff>(sfhm_nrows) * static_cast<std::streamoff>(sizeof(inc_t)), std::ios::cur);
}

void Sketch::make_rho_partial()
{
  if (frac) {
    rho *= (static_cast<double>(r) + 1.0) / static_cast<double>(m);
  } else {
    rho *= 1.0 / static_cast<double>(m);
  }
}

sfhm_sptr_t Sketch::get_sfhm_sptr() { return sfhm; }

lshf_sptr_t Sketch::get_lshf() { return lshf; }

double Sketch::get_rho() const { return rho; }

bool Sketch::check_partial(uint32_t rix) const
{
  uint32_t rix_res = rix % m;
  return (frac && (rix_res <= r)) || (rix_res == r);
}

uint32_t Sketch::search_mer(uint32_t rix, enc_t enc_lr)
{
  uint32_t offset;
  if (frac) {
    offset = ((rix / m) * (r + 1)) + (rix % m);
  } else {
    offset = rix / m;
  }
  vec_enc_it ix1 = sfhm->bucket_iter_start(offset); // TODO: Use pointers?
  vec_enc_it ix2 = sfhm->bucket_iter_next(offset);
  uint32_t hdist_min = std::numeric_limits<uint32_t>::max();
  for (; ix1 < ix2; ++ix1) {
    const uint32_t hdist_curr = popcount_lr32((*ix1) ^ enc_lr);
    if (hdist_curr < hdist_min) {
      hdist_min = hdist_curr;
      if (hdist_curr == 0) break; // Perfect match, can't improve further
    }
  }
  return hdist_min;
}

bool Sketch::search_mer_partial(uint32_t rix, enc_t enc_lr, uint32_t& hdist_min)
{
  const uint32_t rix_res = rix % m;
  // if (__builtin_expect((frac ? (rix_res > r) : (rix_res != r)), 0)) return false;
  if (frac ? (rix_res > r) : (rix_res != r)) return false;

  const uint32_t offset = frac ? ((rix / m) * (r + 1)) + rix_res : rix / m;
  vec_enc_it ix1 = sfhm->bucket_iter_start(offset);
  vec_enc_it ix2 = sfhm->bucket_iter_next(offset);
  hdist_min = std::numeric_limits<uint32_t>::max();
  for (; ix1 < ix2; ++ix1) {
    const uint32_t hdist_curr = popcount_lr32((*ix1) ^ enc_lr);
    if (hdist_curr < hdist_min) {
      hdist_min = hdist_curr;
      if (hdist_curr == 0) break; // HD=0, can't improve further
    }
  }
  return true;
}

uint32_t Sketch::partial_offset(uint32_t rix) const noexcept
{
  const uint32_t rix_res = rix % m;
  if (frac ? (rix_res > r) : (rix_res != r)) {
    return std::numeric_limits<uint32_t>::max();
  }
  return frac ? ((rix / m) * (r + 1)) + rix_res : rix / m;
}

void Sketch::prefetch_offset_inc(uint32_t offset) const noexcept
{
  if (offset != std::numeric_limits<uint32_t>::max()) {
    sfhm->prefetch_inc(offset);
  }
}

void Sketch::prefetch_offset_enc(uint32_t offset) const noexcept
{
  // inc_v must already be in cache for this to be effective
  if (offset != std::numeric_limits<uint32_t>::max()) {
    sfhm->prefetch_enc(offset);
  }
}

bool Sketch::scan_bucket(uint32_t offset, enc_t enc_lr, uint32_t& hdist_min) const noexcept
{
  if (offset == std::numeric_limits<uint32_t>::max()) return false;
  const enc_t* ix1 = sfhm->bucket_ptr_start(offset);
  const enc_t* ix2 = sfhm->bucket_ptr_next(offset);
  hdist_min = std::numeric_limits<uint32_t>::max();
  for (; ix1 < ix2; ++ix1) {
    const uint32_t hd = popcount_lr32((*ix1) ^ enc_lr);
    if (hd < hdist_min) {
      hdist_min = hd;
      if (hd == 0) break; // perfect match, can't improve further
    }
  }
  return true;
}

std::pair<vec_enc_it, vec_enc_it> Sketch::bucket_indices(uint32_t rix)
{
  uint32_t offset;
  if (frac) {
    offset = ((rix / m) * (r + 1)) + (rix % m);
  } else {
    offset = rix / m;
  }
  return std::make_pair(sfhm->bucket_iter_start(offset), sfhm->bucket_iter_next(offset));
}

void Sketch::canonicalize()
{
  if (canonical) {
    return;
  }
  const uint64_t mask_bp = std::numeric_limits<uint64_t>::max() >> ((32 - k) * 2);
  sdhm_sptr_t sdhm = std::make_shared<SDHM>();
  sdhm->enc_vvec.resize(nrows);
  for (uint32_t off = 0; off < nrows; ++off) {
    const enc_t* ix1 = sfhm->bucket_ptr_start(off);
    const enc_t* ix2 = sfhm->bucket_ptr_next(off);
    const uint32_t fw_rix = frac ? (off / (r + 1)) * m + (off % (r + 1)) : off * m + r;
    for (; ix1 < ix2; ++ix1) {
      const uint64_t bp_ppos = lshf->inv_compute_hash(fw_rix);
      const uint64_t bp_npos = lr64_to_bp64(lshf->inv_drop_ppos_lr(*ix1));
      const uint64_t fw_bp = (bp_ppos | bp_npos) & mask_bp;
      const uint64_t rc_bp = revcomp_bp64(fw_bp, k);
      const uint64_t can_bp = std::max(fw_bp, rc_bp);
      const uint32_t rixn = lshf->compute_hash(can_bp);
      const uint32_t resn = rixn % m;
      const bool keep = frac ? (resn <= r) : (resn == r);
      if (!keep) continue;
      const uint32_t new_off = frac ? (rixn / m) * (r + 1) + resn : rixn / m;
      const enc_t new_enc = lshf->drop_ppos_lr(bp64_to_lr64(can_bp));
      sdhm->enc_vvec[new_off].push_back(new_enc);
    }
  }
  sdhm->sort_columns();
  sdhm->make_unique();
  sfhm = std::make_shared<SFHM>(sdhm);
  canonical = true;
}
