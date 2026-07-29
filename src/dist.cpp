#include "dist.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <simde/x86/avx512.h>

#include "common.hpp"
#include "enc.hpp"
#include "msg.hpp"
#include "random.hpp"

extern uint32_t num_threads;

namespace {

  double linear_quantile(const vec<double>& v, const double p)
  {
    if (v.empty()) return nanx();
    if (v.size() == 1) return v.front();
    const double ix = p * static_cast<double>(v.size() - 1);
    const size_t lo = static_cast<size_t>(std::floor(ix));
    const size_t hi = static_cast<size_t>(std::ceil(ix));
    return v[lo] + (ix - static_cast<double>(lo)) * (v[hi] - v[lo]);
  }

  // Branchless minimum Hamming distance of enc against a bucket; identical
  // result to the early-exit scan but without data-dependent branches.
  inline uint32_t bucket_hdist_min(const enc_t* ix1, const enc_t* ix2, const enc_t enc)
  {
    uint32_t hmin = std::numeric_limits<uint32_t>::max();
    for (; ix1 < ix2; ++ix1) {
      const uint32_t hd = popcount_lr32((*ix1) ^ enc);
      hmin = hd < hmin ? hd : hmin;
    }
    return hmin;
  }

  struct scan_ctx_t
  {
    const Sketch* sketch;
    const LSHF* lshf;
    const SFHM* sfhm;
    uint32_t k;
    uint64_t mask_bp;
    uint64_t mask_lr;
    uint64_t bin_shift;
    uint32_t hdist_th;
  };

  // Block size of the software pipeline: hashes/offsets for scan_blk mers are
  // resolved first while their prefetches are in flight, hiding memory latency.
  constexpr size_t scan_blk = 64;

  // Scans mers whose start position j lies in [j0, j1) (j1 <= enmers).
  // agg(bin_j, hdist, is_rc) is invoked for every observed k-mer (valid LSH bucket).
  // STRAND_AWARE: query both strands; otherwise query only the canonical one.
  template<bool STRAND_AWARE, typename Agg>
  inline void scan_mers_range(const scan_ctx_t& ctx, const char* cseq, const uint64_t j0, const uint64_t j1, Agg&& agg)
  {
    const uint32_t k = ctx.k;
    const uint64_t i1 = j1 + k - 1;
    const LSHF* lshf = ctx.lshf;
    const Sketch* sketch = ctx.sketch;
    const SFHM* sfhm = ctx.sfhm;
    uint64_t l = 0;
    uint64_t enc_lr = 0, enc_bp = 0;
    size_t n = 0;
    uint64_t b_bin[scan_blk];
    uint32_t b_off[2][scan_blk];
    enc_t b_enc[2][scan_blk];

    auto flush = [&]() {
      // Phase A: prefetch the bucket-boundary lines for the whole block.
      for (size_t s = 0; s < n; ++s) {
        if (b_off[0][s] != Sketch::OFF_INVALID) sfhm->prefetch_inc(b_off[0][s]);
        if constexpr (STRAND_AWARE) {
          if (b_off[1][s] != Sketch::OFF_INVALID) sfhm->prefetch_inc(b_off[1][s]);
        }
      }
      // Phase B: resolve bucket bounds (now cache-warm) and prefetch enc lines.
      const enc_t* beg[2][scan_blk];
      const enc_t* end[2][scan_blk];
      for (size_t s = 0; s < n; ++s) {
        const uint32_t off0 = b_off[0][s];
        if (off0 == Sketch::OFF_INVALID) {
          beg[0][s] = end[0][s] = nullptr;
        } else {
          beg[0][s] = sfhm->bucket_ptr_start(off0);
          end[0][s] = sfhm->bucket_ptr_next(off0);
          if (beg[0][s] < end[0][s]) __builtin_prefetch(beg[0][s], 0, 0);
        }
        if constexpr (STRAND_AWARE) {
          const uint32_t off1 = b_off[1][s];
          if (off1 == Sketch::OFF_INVALID) {
            beg[1][s] = end[1][s] = nullptr;
          } else {
            beg[1][s] = sfhm->bucket_ptr_start(off1);
            end[1][s] = sfhm->bucket_ptr_next(off1);
            if (beg[1][s] < end[1][s]) __builtin_prefetch(beg[1][s], 0, 0);
          }
        }
      }
      // Phase C: scan buckets and aggregate observed k-mers (hits and misses).
      for (size_t s = 0; s < n; ++s) {
        if (beg[0][s] != nullptr) {
          const uint32_t hd = bucket_hdist_min(beg[0][s], end[0][s], b_enc[0][s]);
          agg(b_bin[s], hd, false);
        }
        if constexpr (STRAND_AWARE) {
          if (beg[1][s] != nullptr) {
            const uint32_t hd = bucket_hdist_min(beg[1][s], end[1][s], b_enc[1][s]);
            agg(b_bin[s], hd, true);
          }
        }
      }
      n = 0;
    };

    for (uint64_t i = j0; i < i1; ++i) {
      if (__builtin_expect(SEQ_NT4_TABLE[static_cast<uint8_t>(cseq[i])] >= 4, 0)) {
        l = 0;
        continue;
      }
      ++l;
      if (l < k) continue;
      const uint64_t j = i - k + 1;
      if (l == k) {
        compute_encoding(cseq + j, cseq + i + 1, enc_lr, enc_bp);
      } else {
        update_encoding(cseq + i, enc_lr, enc_bp);
      }
      enc_bp &= ctx.mask_bp;
      enc_lr &= ctx.mask_lr;
      if (__builtin_expect(j >= j1, 0)) break;
      const uint64_t rc_bp = revcomp_bp64(enc_bp, k);
      if constexpr (STRAND_AWARE) {
        b_off[0][n] = sketch->partial_offset(lshf->compute_hash(enc_bp));
        b_off[1][n] = sketch->partial_offset(lshf->compute_hash(rc_bp));
        b_enc[0][n] = lshf->drop_ppos_lr(enc_lr);
        b_enc[1][n] = lshf->drop_ppos_lr(bp64_to_lr64(rc_bp));
      } else {
        if (rc_bp < enc_bp) {
          b_off[0][n] = sketch->partial_offset(lshf->compute_hash(enc_bp));
          b_enc[0][n] = lshf->drop_ppos_lr(enc_lr);
        } else {
          b_off[0][n] = sketch->partial_offset(lshf->compute_hash(rc_bp));
          b_enc[0][n] = lshf->drop_ppos_lr(bp64_to_lr64(rc_bp));
        }
      }
      b_bin[n] = j >> ctx.bin_shift;
      if (++n == scan_blk) flush();
    }
    if (n) flush();
  }

} // namespace

HDHist::HDHist(const uint64_t nbins, const uint32_t hdist_th, const uint64_t bin_shift, bool zero)
  : nbins(nbins)
  , hdist_th(hdist_th)
  , bin_shift(bin_shift)
  , hdisthist_v(new uint64_t[(nbins + 1) * (hdist_th + 1)])
  , miss_v(new uint64_t[nbins + 1])
{
  if (zero) {
    std::fill_n(hdisthist_v.get(), (nbins + 1) * (hdist_th + 1), uint64_t(0));
    std::fill_n(miss_v.get(), nbins + 1, uint64_t(0));
  }
}

void HDHist::zero_range(const uint64_t r0, const uint64_t r1)
{
  std::fill_n(hdisthist_v.get() + r0 * (hdist_th + 1), (r1 - r0) * (hdist_th + 1), uint64_t(0));
  std::fill_n(miss_v.get() + r0, r1 - r0, uint64_t(0));
}

void HDHist::aggregate_mer(const uint32_t hdist_min, const uint64_t i)
{
  if (i >= nbins) return;
  if (hdist_min <= hdist_th) {
    ++hdisthist_v[((i + 1) * (hdist_th + 1)) + hdist_min];
  } else {
    ++miss_v[i + 1];
  }
}

void HDHist::aggregate_mer_atomic(const uint32_t hdist_min, const uint64_t i)
{
  if (i >= nbins) return;
  if (hdist_min <= hdist_th) {
    __atomic_add_fetch(&hdisthist_v[((i + 1) * (hdist_th + 1)) + hdist_min], 1, __ATOMIC_RELAXED);
  } else {
    __atomic_add_fetch(&miss_v[i + 1], 1, __ATOMIC_RELAXED);
  }
}

void HDHist::compute_prefhistsum()
{
  const uint32_t W = hdist_th + 1;
  for (uint64_t i = 0; i < nbins; ++i) {
    for (uint32_t d = 0; d < W; ++d) {
      hdisthist_v[((i + 1) * W) + d] += hdisthist_v[(i * W) + d];
    }
    miss_v[i + 1] += miss_v[i];
  }
}

void HDHist::compute_prefhistsum_parallel(ThreadPool& pool, uint32_t nchunks)
{
  const uint32_t W = hdist_th + 1;
  const uint64_t n_rows = nbins + 1;
  if (nchunks <= 1 || n_rows < (uint64_t(1) << 16)) {
    compute_prefhistsum();
    return;
  }
  nchunks = std::min<uint32_t>(nchunks, static_cast<uint32_t>((n_rows + 4095) / 4096));
  if (nchunks <= 1) {
    compute_prefhistsum();
    return;
  }
  const uint64_t rows_per = (n_rows + nchunks - 1) / nchunks;
  uint64_t* h = hdisthist_v.get();
  // Phase A: exclusive-local prefix sums inside each chunk.
  pool.parallel_for(nchunks, 1, [&](uint64_t c) {
    const uint64_t c0 = c * rows_per;
    const uint64_t c1 = std::min(c0 + rows_per, n_rows);
    for (uint64_t r = c0 + 1; r < c1; ++r) {
      for (uint32_t d = 0; d < W; ++d)
        h[r * W + d] += h[(r - 1) * W + d];
    }
  });
  // Serial combine of chunk bases (few chunks).
  vec<uint64_t> base(static_cast<uint64_t>(nchunks) * W, 0);
  for (uint32_t c = 1; c < nchunks; ++c) {
    const uint64_t last = std::min((c * rows_per), n_rows) - 1;
    for (uint32_t d = 0; d < W; ++d)
      base[c * W + d] = base[(c - 1) * W + d] + h[last * W + d];
  }
  // Phase B: add each chunk's base to all of its rows.
  pool.parallel_for(nchunks - 1, 1, [&](uint64_t ci) {
    const uint64_t c = ci + 1;
    const uint64_t c0 = c * rows_per;
    const uint64_t c1 = std::min(c0 + rows_per, n_rows);
    const uint64_t* b = &base[c * W];
    for (uint64_t r = c0; r < c1; ++r) {
      for (uint32_t d = 0; d < W; ++d)
        h[r * W + d] += b[d];
    }
  });
  for (uint64_t i = 0; i < nbins; ++i)
    miss_v[i + 1] += miss_v[i];
}

void HDHist::extract_histogram(const uint64_t a, const uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const
{
  assert(a <= b && b <= nbins);
  const uint32_t W = hdist_th + 1;
  v.resize(hdist_bound + 1);
  const simde__mmask8 mask = static_cast<simde__mmask8>((1u << W) - 1);
  const simde__m512i vb = simde_mm512_maskz_loadu_epi64(mask, &hdisthist_v[b * W]);
  const simde__m512i va = simde_mm512_maskz_loadu_epi64(mask, &hdisthist_v[a * W]);
  const simde__m512i vd = simde_mm512_sub_epi64(vb, va);
  simde_mm512_storeu_si512(v.data(), vd);
  const simde__m256i lend = simde_mm512_castsi512_si256(vd);
  const simde__m256i rend = simde_mm512_extracti64x4_epi64(vd, 1);
  const simde__m256i s4 = simde_mm256_add_epi64(lend, rend);
  const simde__m128i s4_lend = simde_mm256_castsi256_si128(s4);
  const simde__m128i s4_rend = simde_mm256_extracti128_si256(s4, 1);
  const simde__m128i s2 = simde_mm_add_epi64(s4_lend, s4_rend);
  t = static_cast<uint64_t>(simde_mm_extract_epi64(s2, 0) + simde_mm_extract_epi64(s2, 1));
  u = miss_v[b] - miss_v[a];
}

dist_summary_t summarize_distances(vec<double> d_v, const vec<uint8_t>* nomatch_v)
{
  dist_summary_t summary;
  assert(nomatch_v == nullptr || nomatch_v->size() == d_v.size());
  vec<double> rel_v; // distances of matched (non-ceiling) windows
  size_t wi = 0;
  for (size_t i = 0; i < d_v.size(); ++i) {
    if (!std::isfinite(d_v[i])) continue;
    d_v[wi++] = d_v[i];
    if (nomatch_v != nullptr) {
      if ((*nomatch_v)[i])
        ++summary.n_nomatch;
      else
        rel_v.push_back(d_v[i]);
    }
  }
  d_v.resize(wi);
  if (d_v.empty()) {
    summary.quantiles.fill(nanx());
    return summary;
  }
  std::sort(d_v.begin(), d_v.end());
  summary.n = d_v.size();
  summary.mean = std::accumulate(d_v.begin(), d_v.end(), 0.0) / static_cast<double>(summary.n);
  if (summary.n > 1) {
    double sum_sq = 0.0;
    for (const double d : d_v) {
      const double delta = d - summary.mean;
      sum_sq += delta * delta;
    }
    summary.sd = std::sqrt(sum_sq / static_cast<double>(summary.n - 1));
  } else {
    summary.sd = 0.0;
  }
  constexpr arr<double, 7> probs{0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99};
  for (size_t i = 0; i < probs.size(); ++i)
    summary.quantiles[i] = linear_quantile(d_v, probs[i]);
  if (nomatch_v != nullptr) {
    std::sort(rel_v.begin(), rel_v.end());
    summary.capped_median = linear_quantile(rel_v, 0.5); // NaN when all capped
  }
  return summary;
}

std::pair<double, char> select_strand_distance(const double d_fw, const double d_rc)
{
  const bool fw_valid = std::isfinite(d_fw);
  const bool rc_valid = std::isfinite(d_rc);
  if (!fw_valid && !rc_valid) return {nanx(), '.'};
  if (!rc_valid || (fw_valid && d_fw <= d_rc)) return {d_fw, '+'};
  return {d_rc, '-'};
}

vec<uint64_t> sample_region_starts(const uint64_t npos, const uint64_t n_samples, std::mt19937& rng)
{
  assert(npos >= 1);
  vec<uint64_t> starts;
  starts.reserve(n_samples);
  std::uniform_int_distribution<uint64_t> rstart(0, npos - 1);
  for (uint64_t i = 0; i < n_samples; ++i)
    starts.push_back(rstart(rng));
  return starts;
}

vec<size_t> select_reservoir_slots(const uint64_t offset, const uint64_t size, const uint64_t total, std::mt19937& rng)
{
  vec<size_t> slots;
  if (size == 0 || total == 0) return slots;
  const uint64_t cum = offset + size;
  slots.reserve((total * size) / cum + 1);
  std::uniform_real_distribution<double> runif(0.0, 1.0);
  for (size_t i = 0; i < total; ++i) {
    if (runif(rng) * static_cast<double>(cum) < static_cast<double>(size)) slots.push_back(i);
  }
  return slots;
}

DistSC::DistSC(CLI::App& sc)
{
  sc.add_option("-q,--query-path", query_path, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible)")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("-i,--sketch-path", sketch_path, "Sketch file at <path> to query")->required()->check(CLI::ExistingFile);
  sc.add_option("-l,--length", tau, "Length of sampled query regions")->required()->check(CLI::PositiveNumber);
  sc.add_option("-b,--bin-shift", bin_shift, "Group consecutive k-mers into bins of size 2^b [0]")->check(CLI::Range(0, 62));
  sc.add_option("--sample-size", sample_size, "Regions sampled across the query file per reference [200]")
    ->check(CLI::PositiveNumber);
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match [4]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("-o,--output-path", output_path, "Write summary output to a file at <path> [stdout]");
  sc.add_option("--samples-output", samples_output_path, "Write sampled regions and distances to a TSV file");
  sc.callback([&]() {
    if (!validate_configuration()) {
      error_exit("Invalid configuration!");
    }
    if (!output_path.empty()) {
      output_file.open(output_path);
      check_fstream(output_file, "Cannot open output file", output_path.string());
      output_stream = &output_file;
    }
    if (!samples_output_path.empty()) {
      samples_output_file.open(samples_output_path);
      check_fstream(samples_output_file, "Cannot open samples output file", samples_output_path.string());
      samples_output_stream = &samples_output_file;
    }
  });
}

bool DistSC::validate_configuration()
{
  bool is_invalid = false;
  if (bin_shift >= 63) {
    is_invalid = true;
    cerr_msg("--bin-shift must be less than 63; got ", bin_shift);
  }
  if (!output_path.empty() && !samples_output_path.empty() && output_path == samples_output_path) {
    is_invalid = true;
    cerr_msg("--output-path and --samples-output must be different files");
  }
  return !is_invalid;
}

void DistSC::dist()
{
  qseq_sptr_t qs = std::make_shared<QSeq>(query_path);
  bool cont_reading;
  while ((cont_reading = qs->read_next_batch())) {
  }

  std::ifstream sketch_stream(sketch_path, std::ifstream::binary);
  check_fstream(sketch_stream, "Cannot open sketch file", sketch_path.string());

  uint32_t nsketches = 0;
  sketch_stream.read(reinterpret_cast<char*>(&nsketches), sizeof(uint32_t));

  vec<uint64_t> sketch_offsets(nsketches);
  for (uint32_t i = 0; i < nsketches; ++i) {
    sketch_offsets[i] = static_cast<uint64_t>(sketch_stream.tellg());
    Sketch::seek_past(sketch_stream);
  }
  sketch_stream.close();

  ThreadPool pool(num_threads);
  cerr_msg("Processing ", nsketches, " sketches w/ ", pool.size(), " thread(s)...");

  init_thread_rng(1);

  const auto& seq_batch = qs->get_seq_batch();
  const auto& qid_batch = qs->get_qid_batch();

  std::ifstream sin(sketch_path, std::ifstream::binary);
  check_fstream(sin, "Cannot open sketch file for reading", sketch_path.string());

  for (uint32_t i = 0; i < nsketches; ++i) {
    sketch_sptr_t sketch = std::make_shared<Sketch>(sketch_path);
    sketch->load_from_offset(sin, sketch_offsets[i]);

    strstream sout;
    strstream samples_sout;
    sample_sequences(sketch, seq_batch, qid_batch, sout, samples_output_stream ? &samples_sout : nullptr, pool);

    if (sout.tellp() > 0) *output_stream << sout.rdbuf();
    if (samples_output_stream && samples_sout.tellp() > 0) *samples_output_stream << samples_sout.rdbuf();

    std::cerr << "\rProcessed sketch " << i + 1 << "/" << nsketches << "..." << std::flush;
    if (i + 1 == nsketches) std::cerr << std::endl;
  }
}

namespace {

  // Per-sequence sampling plan; all randomness is resolved before the parallel
  // phases begin, so results are independent of scheduling.
  struct seq_plan_t
  {
    uint64_t bix = 0;
    uint64_t n_samples = 0;
    bool full_pass = false;
    uint64_t enmers = 0;
    uint64_t nbins = 0;
    bool big = false;
    vec<uint64_t> starts;
    vec<double> d_res;
    vec<char> strand_res;
    vec<uint64_t> hits_res; // matching k-mers (t) of the selected strand; 0 = ceiling
    std::unique_ptr<HDHist> hist_fw;
    std::unique_ptr<HDHist> hist_rc;
  };

  enum task_kind_t : uint8_t
  {
    TASK_HIST_CHUNK, // mer range of a big full-pass plan
    TASK_SAMPLES,    // sample range (window plans, or eval for big full-pass)
    TASK_FUSED       // whole small full-pass sequence: build + prefix + eval
  };

  struct task_t
  {
    uint32_t plan;
    uint64_t a, b;
    task_kind_t kind;
  };

  struct hist_agg_t
  {
    HDHist* fw;
    HDHist* rc;
    bool atomic;
    inline void operator()(uint64_t bin, uint32_t hd, bool is_rc) const
    {
      HDHist* h = is_rc ? rc : fw;
      if (atomic)
        h->aggregate_mer_atomic(hd, bin);
      else
        h->aggregate_mer(hd, bin);
    }
  };

  struct window_agg_t
  {
    uint64_t* fw;
    uint64_t* rc;
    uint64_t* u_fw;
    uint64_t* u_rc;
    uint32_t hdist_th;
    inline void operator()(uint64_t, uint32_t hd, bool is_rc) const
    {
      if (hd <= hdist_th)
        ++(is_rc ? rc[hd] : fw[hd]);
      else
        ++(is_rc ? *u_rc : *u_fw);
    }
  };

} // namespace

void DistSC::sample_sequences(const sketch_sptr_t& sketch,
                              const vec<str>& seq_batch,
                              const vec<str>& qid_batch,
                              strstream& sout,
                              strstream* samples_sout,
                              ThreadPool& pool)
{
  const lshf_sptr_t lshf = sketch->get_lshf();
  const uint32_t k = lshf->get_k();
  if (tau < k) {
    error_exit(concat_msg("--length must be at least the sketch k-mer length (length=", tau, ", k=", k, ")"));
  }
  const uint64_t bin_size = uint64_t(1) << bin_shift;
  if (bin_size > tau) {
    error_exit(concat_msg("--bin-shift gives bin_size=", bin_size, ", which exceeds --length=", tau));
  }
  const uint64_t region_nmers = tau - k + 1;
  const uint64_t tau_bin = std::max<uint64_t>(1, (region_nmers + bin_size - 1) >> bin_shift);
  const bool canonical = sketch->is_canonical();
  const uint32_t nworkers = pool.size();

  uint64_t total_len = 0;
  vec<uint64_t> cum_lens(seq_batch.size());
  for (size_t bix = 0; bix < seq_batch.size(); ++bix) {
    if (seq_batch[bix].size() >= tau) total_len += seq_batch[bix].size();
    cum_lens[bix] = total_len;
  }

  if (total_len == 0) {
    const dist_summary_t summary = summarize_distances({});
    sout << std::setprecision(8);
    write_tsv(sout, query_path, sketch->get_rid(), summary.n, summary.mean, summary.sd);
    for (const double q : summary.quantiles)
      sout << '\t' << q;
    sout << '\t' << summary.capped_median << '\t' << summary.n_nomatch << '\n';
    return;
  }

  const bool prof = std::getenv("GDIFF_PROF") != nullptr;
  auto t0 = std::chrono::steady_clock::now();

  // File-wide uniform positions, then per-sequence sample counts (identical
  // RNG consumption to a sequential implementation).
  vec<uint64_t> positions;
  positions.reserve(sample_size);
  std::uniform_int_distribution<uint64_t> rpos(0, total_len - 1);
  for (uint64_t i = 0; i < sample_size; ++i)
    positions.push_back(rpos(gen));
  std::sort(positions.begin(), positions.end());

  // Build per-sequence plans and pre-draw region starts.
  vec<seq_plan_t> plans;
  plans.reserve(seq_batch.size());
  size_t pidx = 0;
  uint32_t n_fullpass = 0;
  for (size_t bix = 0; bix < seq_batch.size() && pidx < positions.size(); ++bix) {
    const uint64_t L = seq_batch[bix].size();
    if (L < tau) continue;
    const uint64_t seq_end = cum_lens[bix];
    uint64_t n_for_seq = 0;
    while (pidx < positions.size() && positions[pidx] < seq_end) {
      ++n_for_seq;
      ++pidx;
    }
    if (n_for_seq == 0) continue;
    const uint64_t enmers = L - k + 1;
    const uint64_t nbins = (enmers + bin_size - 1) >> bin_shift;
    if (nbins < tau_bin) continue; // defensive; cannot happen when L >= tau
    seq_plan_t p;
    p.bix = bix;
    p.n_samples = n_for_seq;
    p.enmers = enmers;
    p.nbins = nbins;
    p.full_pass = n_for_seq * region_nmers >= enmers;
    n_fullpass += p.full_pass;
    const uint64_t npos = nbins - tau_bin + 1;
    p.starts = sample_region_starts(npos, n_for_seq, gen);
    p.d_res.assign(n_for_seq, nanx());
    p.strand_res.assign(n_for_seq, '+');
    p.hits_res.assign(n_for_seq, 0);
    plans.push_back(std::move(p));
  }
  // Full-pass sequences become "big" (parallel chunked build) only when there
  // aren't enough of them to keep workers busy; others use fused per-sequence tasks.
  const uint64_t big_thresh = uint64_t(1) << 21;
  uint32_t n_big = 0;
  for (auto& p : plans) {
    if (p.full_pass && p.enmers >= big_thresh && n_fullpass < nworkers * 2) {
      p.big = true;
      ++n_big;
    }
  }
  for (auto& p : plans) {
    if (!p.big) continue;
    p.hist_fw = std::make_unique<HDHist>(p.nbins, hdist_th, bin_shift, /*zero=*/false);
    if (!canonical) p.hist_rc = std::make_unique<HDHist>(p.nbins, hdist_th, bin_shift, /*zero=*/false);
  }

  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  const scan_ctx_t ctx{sketch.get(),
                       lshf.get(),
                       sketch->get_sfhm_sptr().get(),
                       k,
                       u64m >> ((32 - k) * 2),
                       ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k)),
                       bin_shift,
                       hdist_th};
  const LLH<double> llhf(k, lshf->get_h(), sketch->get_rho(), hdist_th, 0.0, false);

  using clock = std::chrono::steady_clock;
  auto toc = [&](const char* name, clock::time_point t0) {
    if (prof) std::cerr << "    [" << name << "] " << std::chrono::duration<double>(clock::now() - t0).count() << " s\n";
  };
  if (prof) std::cerr << "    [plan+alloc] " << std::chrono::duration<double>(clock::now() - t0).count() << " s\n";
  auto tpz = clock::now();

  // Zero the big histogram backing stores in parallel (page faults and memsets
  // are spread across workers).
  if (n_big > 0) {
    struct zero_task_t
    {
      uint32_t plan;
      uint64_t r0, r1;
    };
    vec<zero_task_t> zero_tasks;
    for (uint32_t pi = 0; pi < plans.size(); ++pi) {
      const seq_plan_t& p = plans[pi];
      if (!p.big) continue;
      const uint64_t rows = p.hist_fw->storage_rows();
      const uint64_t step = std::max<uint64_t>(rows / (nworkers * 4) + 1, 1u << 16);
      for (uint64_t r0 = 0; r0 < rows; r0 += step)
        zero_tasks.push_back({pi, r0, std::min(r0 + step, rows)});
    }
    pool.parallel_for(zero_tasks.size(), 1, [&](const uint64_t ti) {
      const zero_task_t& z = zero_tasks[ti];
      seq_plan_t& p = plans[z.plan];
      p.hist_fw->zero_range(z.r0, z.r1);
      if (p.hist_rc) p.hist_rc->zero_range(z.r0, z.r1);
    });
  }
  if (prof) toc("zero hists", tpz);
  auto tp0 = clock::now();

  // Per-sample evaluation helper for full-pass histograms (used by fused tasks
  // and by phase 3 for big plans).
  auto eval_fullpass =
    [&](seq_plan_t& p, const HDHist& hist_fw, const HDHist* hist_rc, const uint64_t s0, const uint64_t s1) {
      vec<uint64_t> v_scratch(hdist_bound + 1, 0);
      for (uint64_t s = s0; s < s1; ++s) {
        const uint64_t a_bin = p.starts[s];
        const uint64_t b_bin = a_bin + tau_bin;
        uint64_t u_fw = 0, t_fw = 0;
        hist_fw.extract_histogram(a_bin, b_bin, v_scratch, u_fw, t_fw);
        const double d_fw = validate_distance(llhf.mle_at(v_scratch.data(), u_fw));
        double d = d_fw;
        char strand = '+';
        uint64_t t_sel = t_fw;
        if (!canonical) {
          uint64_t u_rc = 0, t_rc = 0;
          hist_rc->extract_histogram(a_bin, b_bin, v_scratch, u_rc, t_rc);
          const double d_rc = validate_distance(llhf.mle_at(v_scratch.data(), u_rc));
          std::tie(d, strand) = select_strand_distance(d_fw, d_rc);
          t_sel = (strand == '-') ? t_rc : t_fw;
        }
        p.d_res[s] = d;
        p.strand_res[s] = strand;
        p.hits_res[s] = t_sel;
      }
    };

  // Phase 1: window scans, fused small full-pass sequences, and chunked builds
  // for big full-pass sequences.
  {
    vec<task_t> tasks;
    // Fused plans first, largest first, to avoid stragglers.
    vec<uint32_t> fused;
    for (uint32_t pi = 0; pi < plans.size(); ++pi)
      if (plans[pi].full_pass && !plans[pi].big) fused.push_back(pi);
    std::sort(fused.begin(), fused.end(), [&](uint32_t x, uint32_t y) { return plans[x].enmers > plans[y].enmers; });
    for (const uint32_t pi : fused)
      tasks.push_back({pi, 0, 0, TASK_FUSED});
    for (uint32_t pi = 0; pi < plans.size(); ++pi) {
      const seq_plan_t& p = plans[pi];
      if (p.big) {
        const uint64_t chunk_mers = std::max<uint64_t>(uint64_t(1) << 18, (p.enmers + nworkers * 2 - 1) / (nworkers * 2));
        for (uint64_t c0 = 0; c0 < p.enmers; c0 += chunk_mers)
          tasks.push_back({pi, c0, std::min(c0 + chunk_mers, p.enmers), TASK_HIST_CHUNK});
      }
    }
    for (uint32_t pi = 0; pi < plans.size(); ++pi) {
      const seq_plan_t& p = plans[pi];
      if (p.full_pass) continue;
      const uint64_t n_per = std::max<uint64_t>(4, (p.n_samples + nworkers * 4 - 1) / (nworkers * 4));
      for (uint64_t s0 = 0; s0 < p.n_samples; s0 += n_per)
        tasks.push_back({pi, s0, std::min(s0 + n_per, p.n_samples), TASK_SAMPLES});
    }
    pool.parallel_for(tasks.size(), 1, [&](const uint64_t ti) {
      const task_t& t = tasks[ti];
      seq_plan_t& p = plans[t.plan];
      const char* cseq = seq_batch[p.bix].data();
      if (t.kind == TASK_HIST_CHUNK) {
        hist_agg_t agg{p.hist_fw.get(), p.hist_rc.get(), /*atomic=*/true};
        if (canonical)
          scan_mers_range<false>(ctx, cseq, t.a, t.b, agg);
        else
          scan_mers_range<true>(ctx, cseq, t.a, t.b, agg);
      } else if (t.kind == TASK_FUSED) {
        HDHist hist_fw(p.nbins, hdist_th, bin_shift);
        std::unique_ptr<HDHist> hist_rc;
        if (!canonical) hist_rc = std::make_unique<HDHist>(p.nbins, hdist_th, bin_shift);
        hist_agg_t agg{&hist_fw, hist_rc.get(), /*atomic=*/false};
        if (canonical)
          scan_mers_range<false>(ctx, cseq, 0, p.enmers, agg);
        else
          scan_mers_range<true>(ctx, cseq, 0, p.enmers, agg);
        hist_fw.compute_prefhistsum();
        if (hist_rc) hist_rc->compute_prefhistsum();
        eval_fullpass(p, hist_fw, hist_rc.get(), 0, p.n_samples);
      } else {
        uint64_t v_fw[hdist_bound + 1];
        uint64_t v_rc[hdist_bound + 1];
        for (uint64_t s = t.a; s < t.b; ++s) {
          const uint64_t a_bin = p.starts[s];
          const uint64_t j0 = a_bin << bin_shift;
          const uint64_t j1 = std::min(j0 + region_nmers, p.enmers);
          std::fill(v_fw, v_fw + hdist_bound + 1, 0);
          std::fill(v_rc, v_rc + hdist_bound + 1, 0);
          uint64_t u_fw = 0, u_rc = 0;
          window_agg_t agg{v_fw, v_rc, &u_fw, &u_rc, hdist_th};
          if (canonical)
            scan_mers_range<false>(ctx, cseq, j0, j1, agg);
          else
            scan_mers_range<true>(ctx, cseq, j0, j1, agg);
          const uint64_t t_fw = std::accumulate(v_fw, v_fw + hdist_bound + 1, uint64_t(0));
          const double d_fw = validate_distance(llhf.mle_at(v_fw, u_fw));
          double d = d_fw;
          char strand = '+';
          uint64_t t_sel = t_fw;
          if (!canonical) {
            const uint64_t t_rc = std::accumulate(v_rc, v_rc + hdist_bound + 1, uint64_t(0));
            const double d_rc = validate_distance(llhf.mle_at(v_rc, u_rc));
            std::tie(d, strand) = select_strand_distance(d_fw, d_rc);
            t_sel = (strand == '-') ? t_rc : t_fw;
          }
          p.d_res[s] = d;
          p.strand_res[s] = strand;
          p.hits_res[s] = t_sel;
        }
      }
    });
  }

  // Phase 2: prefix sums over the big per-bin histograms.
  toc("phase1 scan", tp0);
  tp0 = clock::now();
  if (n_big > 1) {
    vec<HDHist*> big_hists;
    for (auto& p : plans) {
      if (!p.big) continue;
      big_hists.push_back(p.hist_fw.get());
      if (p.hist_rc) big_hists.push_back(p.hist_rc.get());
    }
    pool.parallel_for(big_hists.size(), 1, [&](uint64_t j) { big_hists[j]->compute_prefhistsum(); });
  } else {
    for (auto& p : plans) {
      if (!p.big) continue;
      p.hist_fw->compute_prefhistsum_parallel(pool, nworkers);
      if (p.hist_rc) p.hist_rc->compute_prefhistsum_parallel(pool, nworkers);
    }
  }
  toc("phase2 prefix", tp0);
  tp0 = clock::now();

  // Phase 3: window histogram extraction + MLE for big full-pass plans.
  if (n_big > 0) {
    vec<task_t> tasks;
    for (uint32_t pi = 0; pi < plans.size(); ++pi) {
      const seq_plan_t& p = plans[pi];
      if (!p.big) continue;
      const uint64_t n_per = std::max<uint64_t>(4, (p.n_samples + nworkers * 4 - 1) / (nworkers * 4));
      for (uint64_t s0 = 0; s0 < p.n_samples; s0 += n_per)
        tasks.push_back({pi, s0, std::min(s0 + n_per, p.n_samples), TASK_SAMPLES});
    }
    pool.parallel_for(tasks.size(), 1, [&](const uint64_t ti) {
      const task_t& t = tasks[ti];
      seq_plan_t& p = plans[t.plan];
      eval_fullpass(p, *p.hist_fw, p.hist_rc.get(), t.a, t.b);
    });
  }

  // Merge per-sequence results in sequence order (deterministic).
  toc("phase3 eval", tp0);
  tp0 = clock::now();
  vec<double> d_v;
  vec<uint8_t> nomatch_v;
  d_v.reserve(sample_size);
  nomatch_v.reserve(sample_size);
  if (samples_sout) *samples_sout << std::setprecision(8);
  for (const auto& p : plans) {
    const str& qid = qid_batch[p.bix];
    for (uint64_t s = 0; s < p.n_samples; ++s) {
      const double d = p.d_res[s];
      if (!std::isfinite(d)) continue;
      d_v.push_back(d);
      nomatch_v.push_back(p.hits_res[s] == 0);
      if (samples_sout) {
        const uint64_t a = p.starts[s] << bin_shift;
        write_tsv(*samples_sout, qid, p.enmers + k - 1, a + 1, a + tau, p.strand_res[s], sketch->get_rid(), d)
          << '\t' << p.hits_res[s] << '\n';
      }
    }
  }

  const dist_summary_t summary = summarize_distances(std::move(d_v), &nomatch_v);
  toc("merge+summarize", tp0);
  sout << std::setprecision(8);
  write_tsv(sout, query_path, sketch->get_rid(), summary.n, summary.mean, summary.sd);
  for (const double q : summary.quantiles)
    sout << '\t' << q;
  sout << '\t' << summary.capped_median << '\t' << summary.n_nomatch << '\n';
}
