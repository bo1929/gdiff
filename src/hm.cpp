#include "hm.hpp"

SFHM::SFHM(sdhm_sptr_t source)
{
  nkmers = source->nkmers;
  nrows = source->enc_vvec.size();
  inc_v.resize(nrows);
  enc_v.reserve(nkmers);
  inc_t limit_inc = std::numeric_limits<inc_t>::max();
  inc_t cpinc;
  inc_t lix = 0;
  for (uint32_t rix = 0; rix < nrows; ++rix) {
    cpinc = std::min(limit_inc, static_cast<inc_t>(source->enc_vvec[rix].size()));
    for (inc_t i = 0; i < cpinc; ++i) {
      enc_v.push_back(source->enc_vvec[rix][i]);
    }
    lix += cpinc;
    inc_v[rix] = lix;
    source->enc_vvec[rix].clear();
  }
}

SFHM::~SFHM()
{
  inc_v.clear();
  enc_v.clear();
}

void SFHM::load(std::ifstream& sketch_stream)
{
  sketch_stream.read(reinterpret_cast<char*>(&nkmers), sizeof(uint64_t));
  enc_v.resize(nkmers);
  sketch_stream.read(reinterpret_cast<char*>(enc_v.data()), nkmers * sizeof(enc_t));
  assert(nkmers == enc_v.size());
  sketch_stream.read(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));
  inc_v.resize(nrows);
  sketch_stream.read(reinterpret_cast<char*>(inc_v.data()), nrows * sizeof(inc_t));
  assert(nrows == inc_v.size());
}

void SFHM::save(std::ofstream& sketch_stream)
{
  sketch_stream.write(reinterpret_cast<const char*>(&nkmers), sizeof(uint64_t));
  sketch_stream.write(reinterpret_cast<const char*>(enc_v.data()), sizeof(enc_t) * nkmers);
  sketch_stream.write(reinterpret_cast<const char*>(&nrows), sizeof(uint32_t));
  sketch_stream.write(reinterpret_cast<const char*>(inc_v.data()), sizeof(inc_t) * nrows);
}

std::vector<enc_t>::const_iterator SFHM::bucket_start(uint32_t rix)
{
  if ((rix != 0) && (rix <= inc_v.size())) {
    return std::next(enc_v.begin(), inc_v[rix - 1]);
  } else {
    return enc_v.begin();
  }
}

std::vector<enc_t>::const_iterator SFHM::bucket_next(uint32_t rix)
{
  if (rix < inc_v.size()) {
    return std::next(enc_v.begin(), inc_v[rix]);
  } else {
    return enc_v.end();
  }
}

void SDHM::sort_columns()
{
  for (uint32_t i = 0; i < enc_vvec.size(); ++i) {
    if (!enc_vvec[i].empty()) {
      std::sort(enc_vvec[i].begin(), enc_vvec[i].end());
    }
  }
}

void SDHM::make_unique()
{
  nkmers = 0;
  for (uint32_t i = 0; i < enc_vvec.size(); ++i) {
    if (!enc_vvec[i].empty()) {
      enc_vvec[i].erase(std::unique(enc_vvec[i].begin(), enc_vvec[i].end()), enc_vvec[i].end());
    }
    nkmers += enc_vvec[i].size();
  }
}

void SDHM::fill_table(uint32_t nrows, rseq_sptr_t rs)
{
  enc_vvec.resize(nrows);
  while (rs->read_next_seq()) {
    if (rs->set_curr_seq()) {
      rs->extract_mers(enc_vvec);
    }
  }
  rs->compute_rho();
  sort_columns();
  make_unique();
}

uint64_t SDHM::get_nmers() { return nkmers; }
