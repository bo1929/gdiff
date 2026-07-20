#include "dist.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <mutex>
#include <numeric>
#include <simde/x86/avx512.h>
#include <thread>

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

} // namespace

HDHist::HDHist(const uint64_t nbins, const uint64_t nmers, const uint32_t hdist_th, const uint64_t bin_shift)
  : nbins(nbins)
  , nmers(nmers)
  , hdist_th(hdist_th)
  , bin_shift(bin_shift)
  , hdisthist_v((nbins + 1) * (hdist_th + 1), 0)
{
}

void HDHist::aggregate_mer(const uint32_t hdist_min, const uint64_t i)
{
  if (hdist_min <= hdist_th && i < nbins) {
    ++hdisthist_v[((i + 1) * (hdist_th + 1)) + hdist_min];
  }
}

void HDHist::compute_prefhistsum()
{
  const uint32_t W = hdist_th + 1;
  for (uint64_t i = 0; i < nbins; ++i) {
    for (uint32_t d = 0; d < W; ++d) {
      hdisthist_v[((i + 1) * W) + d] += hdisthist_v[(i * W) + d];
    }
  }
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
  const uint64_t mers_b = std::min(b << bin_shift, nmers);
  const uint64_t mers_a = std::min(a << bin_shift, nmers);
  u = (mers_b - mers_a) - t;
}

dist_summary_t summarize_distances(vec<double> d_v)
{
  dist_summary_t summary;
  d_v.erase(std::remove_if(d_v.begin(), d_v.end(), [](const double d) { return !std::isfinite(d); }), d_v.end());
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
  const uint32_t nthreads = std::max(1u, std::min(num_threads, nsketches));
  cerr_msg("Processing ", nsketches, " sketches w/ ", nthreads, " thread(s)...");

  vec<uint64_t> sketch_offsets(nsketches);
  for (uint32_t i = 0; i < nsketches; ++i) {
    sketch_offsets[i] = static_cast<uint64_t>(sketch_stream.tellg());
    Sketch::seek_past(sketch_stream);
  }
  sketch_stream.close();

  vec<strstream> results(nsketches);
  vec<strstream> samples_results(samples_output_stream ? nsketches : 0);
  std::atomic<uint32_t> next_idx{0};
  std::atomic<uint32_t> count_p{0};
  std::mutex cerr_mtx;

  const auto& seq_batch = qs->get_seq_batch();
  const auto& qid_batch = qs->get_qid_batch();

  auto worker = [&](const uint32_t tseed) {
    init_thread_rng(tseed);
    uint32_t i;
    while ((i = next_idx.fetch_add(1, std::memory_order_relaxed)) < nsketches) {
      std::ifstream sketch_stream(sketch_path, std::ifstream::binary);
      sketch_sptr_t sketch = std::make_shared<Sketch>(sketch_path);
      sketch->load_from_offset(sketch_stream, sketch_offsets[i]);
      sketch_stream.close();
      sketch->make_rho_partial();

      strstream* samples_sout = samples_output_stream ? &samples_results[i] : nullptr;
      sample_sequences(sketch, seq_batch, qid_batch, results[i], samples_sout);

      const uint32_t num_p = count_p.fetch_add(1, std::memory_order_relaxed) + 1;
      {
        std::lock_guard<std::mutex> lock(cerr_mtx);
        std::cerr << "\rProcessed sketch " << num_p << "/" << nsketches << "..." << std::flush;
        if (num_p == nsketches) std::cerr << std::endl;
      }
    }
  };

  vec<std::thread> threads;
  threads.reserve(nthreads);
  for (uint32_t t = 0; t < nthreads; ++t)
    threads.emplace_back([&, t]() { worker(t + 1); });
  for (auto& t : threads)
    t.join();

  for (uint32_t i = 0; i < nsketches; ++i) {
    if (results[i].tellp() > 0) *output_stream << results[i].rdbuf();
  }
  if (samples_output_stream) {
    for (uint32_t i = 0; i < nsketches; ++i) {
      if (samples_results[i].tellp() > 0) *samples_output_stream << samples_results[i].rdbuf();
    }
  }
}

void DistSC::sample_sequences(const sketch_sptr_t& sketch,
                              const vec<str>& seq_batch,
                              const vec<str>& qid_batch,
                              strstream& sout,
                              strstream* samples_sout)
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

  uint64_t total_len = 0;
  vec<uint64_t> cum_lens(seq_batch.size());
  for (size_t bix = 0; bix < seq_batch.size(); ++bix) {
    if (seq_batch[bix].size() >= tau) total_len += seq_batch[bix].size();
    cum_lens[bix] = total_len;
  }

  vec<double> d_v;
  d_v.reserve(sample_size);
  vec<dist_sample_t> samples_v;
  vec<dist_sample_t>* samples_ptr = nullptr;
  if (samples_sout) {
    samples_v.reserve(sample_size);
    samples_ptr = &samples_v;
  }

  if (total_len == 0) {
    const dist_summary_t summary = summarize_distances(std::move(d_v));
    sout << std::setprecision(8);
    write_tsv(sout, query_path, sketch->get_rid(), summary.n, summary.mean, summary.sd);
    for (const double q : summary.quantiles)
      sout << '\t' << q;
    sout << '\n';
    return;
  }

  vec<uint64_t> positions;
  positions.reserve(sample_size);
  std::uniform_int_distribution<uint64_t> rpos(0, total_len - 1);
  for (uint64_t i = 0; i < sample_size; ++i)
    positions.push_back(rpos(gen));
  std::sort(positions.begin(), positions.end());

  size_t pidx = 0;
  for (size_t bix = 0; bix < seq_batch.size() && pidx < positions.size(); ++bix) {
    const str& seq = seq_batch[bix];
    if (seq.size() < tau) continue;
    const uint64_t seq_end = cum_lens[bix];
    uint64_t n_for_seq = 0;
    while (pidx < positions.size() && positions[pidx] < seq_end) {
      ++n_for_seq;
      ++pidx;
    }
    if (n_for_seq > 0) sample_sequence(sketch, seq, qid_batch[bix], n_for_seq, d_v, samples_ptr);
  }

  if (samples_ptr) {
    *samples_sout << std::setprecision(8);
    for (const auto& s : samples_v) {
      if (!std::isfinite(s.d)) continue;
      write_tsv(*samples_sout, s.qid, s.L, s.a + 1, s.a + tau, s.strand, sketch->get_rid(), s.d) << '\n';
    }
  }

  const dist_summary_t summary = summarize_distances(std::move(d_v));
  sout << std::setprecision(8);
  write_tsv(sout, query_path, sketch->get_rid(), summary.n, summary.mean, summary.sd);
  for (const double q : summary.quantiles)
    sout << '\t' << q;
  sout << '\n';
}

void DistSC::sample_sequence(const sketch_sptr_t& sketch,
                             const str& seq,
                             const str& qid,
                             uint64_t n_samples,
                             vec<double>& d_v,
                             vec<dist_sample_t>* samples_v)
{
  const lshf_sptr_t lshf = sketch->get_lshf();
  const uint32_t k = lshf->get_k();
  const uint64_t L = seq.size();
  const uint64_t enmers = L - k + 1;
  const uint64_t region_nmers = tau - k + 1;
  const uint64_t bin_size = uint64_t(1) << bin_shift;
  const uint64_t nbins = (enmers + bin_size - 1) >> bin_shift;
  const uint64_t tau_bin = std::max<uint64_t>(1, (region_nmers + bin_size - 1) >> bin_shift);
  if (nbins < tau_bin) return;

  const uint64_t npos = nbins - tau_bin + 1;
  std::uniform_int_distribution<uint64_t> rvstart(0, npos - 1);
  LLH<double> llhf(k, lshf->get_h(), sketch->get_rho(), hdist_th, 0.0, false);
  const bool canonical = sketch->is_canonical();

  const bool full_pass = n_samples * region_nmers >= enmers;
  std::unique_ptr<HDHist> hist_fw;
  std::unique_ptr<HDHist> hist_rc;
  if (full_pass) {
    if (canonical) {
      hist_fw = std::make_unique<HDHist>(nbins, enmers, hdist_th, bin_shift);
      search_mers(sketch, seq.data(), L, *hist_fw);
      hist_fw->compute_prefhistsum();
    } else {
      hist_fw = std::make_unique<HDHist>(nbins, enmers, hdist_th, bin_shift);
      hist_rc = std::make_unique<HDHist>(nbins, enmers, hdist_th, bin_shift);
      search_mers(sketch, seq.data(), L, *hist_fw, *hist_rc);
      hist_fw->compute_prefhistsum();
      hist_rc->compute_prefhistsum();
    }
  }

  vec<uint64_t> v_scratch(hdist_bound + 1, 0);
  vec<uint64_t> v_rc(hdist_bound + 1, 0);

  for (uint64_t i = 0; i < n_samples; ++i) {
    const uint64_t a_bin = rvstart(gen);
    const uint64_t b_bin = a_bin + tau_bin;
    uint64_t u = 0, t = 0;
    double d_fw = nanx(), d_rc = nanx();

    if (full_pass) {
      hist_fw->extract_histogram(a_bin, b_bin, v_scratch, u, t);
      d_fw = validate_distance(llhf.mle(v_scratch.data(), u));
      if (hist_rc) {
        hist_rc->extract_histogram(a_bin, b_bin, v_scratch, u, t);
        d_rc = validate_distance(llhf.mle(v_scratch.data(), u));
      }
    } else {
      const uint64_t j0 = a_bin << bin_shift;
      const uint64_t j1 = std::min(j0 + region_nmers, enmers);
      v_scratch.assign(hdist_bound + 1, 0);
      if (!canonical) v_rc.assign(hdist_bound + 1, 0);
      search_window(sketch, seq.data(), L, j0, j1, v_scratch, canonical ? nullptr : &v_rc);
      t = 0;
      for (uint32_t d = 0; d <= hdist_th; ++d)
        t += v_scratch[d];
      u = (j1 - j0) - t;
      d_fw = validate_distance(llhf.mle(v_scratch.data(), u));
      if (!canonical) {
        t = 0;
        for (uint32_t d = 0; d <= hdist_th; ++d)
          t += v_rc[d];
        u = (j1 - j0) - t;
        d_rc = validate_distance(llhf.mle(v_rc.data(), u));
      }
    }

    double d = d_fw;
    char strand = '+';
    if (!canonical) std::tie(d, strand) = select_strand_distance(d_fw, d_rc);
    if (!std::isfinite(d)) continue;

    d_v.push_back(d);
    if (samples_v) {
      const uint64_t a = a_bin << bin_shift;
      samples_v->push_back({qid, L, a, strand, d});
    }
  }
}

void DistSC::search_mers(const sketch_sptr_t& sketch, const char* cseq, const uint64_t len, HDHist& hist) const
{
  const lshf_sptr_t lshf = sketch->get_lshf();
  const uint32_t k = lshf->get_k();
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  const uint64_t mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  const uint64_t mask_bp = u64m >> ((32 - k) * 2);

  uint64_t i = 0, j = 0, l = 0;
  uint32_t orrix, rcrix;
  uint64_t orenc64_bp, orenc64_lr, rcenc64_bp;
  for (; i < len; ++i) {
    if (__builtin_expect(SEQ_NT4_TABLE[cseq[i]] >= 4, 0)) {
      l = 0;
      continue;
    }
    ++l;
    if (l < k) continue;

    j = i - k + 1;
    if (l == k) {
      compute_encoding(cseq + j, cseq + i + 1, orenc64_lr, orenc64_bp);
    } else {
      update_encoding(cseq + i, orenc64_lr, orenc64_bp);
    }
    orenc64_bp &= mask_bp;
    orenc64_lr &= mask_lr;
    rcenc64_bp = revcomp_bp64(orenc64_bp, k);
    const uint64_t bin_j = j >> bin_shift;

    if (rcenc64_bp < orenc64_bp) {
      orrix = lshf->compute_hash(orenc64_bp);
      const uint32_t off_fw = sketch->partial_offset(orrix);
      sketch->prefetch_offset_inc(off_fw);
      const enc_t enc_lr_fw = lshf->drop_ppos_lr(orenc64_lr);
      sketch->prefetch_offset_enc(off_fw);
      uint32_t hdist_fw;
      if (sketch->scan_bucket(off_fw, enc_lr_fw, hdist_fw)) {
        hist.aggregate_mer(hdist_fw, bin_j);
      }
    } else {
      rcrix = lshf->compute_hash(rcenc64_bp);
      const uint32_t off_rc = sketch->partial_offset(rcrix);
      sketch->prefetch_offset_inc(off_rc);
      const enc_t enc_lr_rc = lshf->drop_ppos_lr(bp64_to_lr64(rcenc64_bp));
      sketch->prefetch_offset_enc(off_rc);
      uint32_t hdist_rc;
      if (sketch->scan_bucket(off_rc, enc_lr_rc, hdist_rc)) {
        hist.aggregate_mer(hdist_rc, bin_j);
      }
    }
  }
}

void DistSC::search_mers(const sketch_sptr_t& sketch,
                         const char* cseq,
                         const uint64_t len,
                         HDHist& hist_fw,
                         HDHist& hist_rc) const
{
  const lshf_sptr_t lshf = sketch->get_lshf();
  const uint32_t k = lshf->get_k();
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  const uint64_t mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  const uint64_t mask_bp = u64m >> ((32 - k) * 2);

  uint64_t i = 0, j = 0, l = 0;
  uint32_t orrix, rcrix;
  uint64_t orenc64_bp, orenc64_lr, rcenc64_bp;
  for (; i < len; ++i) {
    if (__builtin_expect(SEQ_NT4_TABLE[cseq[i]] >= 4, 0)) {
      l = 0;
      continue;
    }
    ++l;
    if (l < k) continue;

    j = i - k + 1;
    if (l == k) {
      compute_encoding(cseq + j, cseq + i + 1, orenc64_lr, orenc64_bp);
    } else {
      update_encoding(cseq + i, orenc64_lr, orenc64_bp);
    }
    orenc64_bp &= mask_bp;
    orenc64_lr &= mask_lr;
    rcenc64_bp = revcomp_bp64(orenc64_bp, k);
    const uint64_t bin_j = j >> bin_shift;

    orrix = lshf->compute_hash(orenc64_bp);
    rcrix = lshf->compute_hash(rcenc64_bp);
    const uint32_t off_fw = sketch->partial_offset(orrix);
    const uint32_t off_rc = sketch->partial_offset(rcrix);
    sketch->prefetch_offset_inc(off_fw);
    sketch->prefetch_offset_inc(off_rc);
    const enc_t enc_lr_fw = lshf->drop_ppos_lr(orenc64_lr);
    const enc_t enc_lr_rc = lshf->drop_ppos_lr(bp64_to_lr64(rcenc64_bp));
    sketch->prefetch_offset_enc(off_fw);
    sketch->prefetch_offset_enc(off_rc);
    uint32_t hdist_fw;
    if (sketch->scan_bucket(off_fw, enc_lr_fw, hdist_fw)) {
      hist_fw.aggregate_mer(hdist_fw, bin_j);
    }
    uint32_t hdist_rc;
    if (sketch->scan_bucket(off_rc, enc_lr_rc, hdist_rc)) {
      hist_rc.aggregate_mer(hdist_rc, bin_j);
    }
  }
}

void DistSC::search_window(const sketch_sptr_t& sketch,
                           const char* cseq,
                           const uint64_t len,
                           const uint64_t j0,
                           const uint64_t j1,
                           vec<uint64_t>& v_fw,
                           vec<uint64_t>* v_rc) const
{
  if (j0 >= j1) return;
  const lshf_sptr_t lshf = sketch->get_lshf();
  const uint32_t k = lshf->get_k();
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  const uint64_t mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  const uint64_t mask_bp = u64m >> ((32 - k) * 2);

  if (j0 + k > len) return;
  uint64_t orenc64_bp = 0, orenc64_lr = 0;
  compute_encoding(cseq + j0, cseq + j0 + k, orenc64_lr, orenc64_bp);
  orenc64_bp &= mask_bp;
  orenc64_lr &= mask_lr;

  for (uint64_t j = j0;;) {
    const uint64_t rcenc64_bp = revcomp_bp64(orenc64_bp, k);
    if (v_rc) {
      uint32_t orrix = lshf->compute_hash(orenc64_bp);
      uint32_t rcrix = lshf->compute_hash(rcenc64_bp);
      const uint32_t off_fw = sketch->partial_offset(orrix);
      const uint32_t off_rc = sketch->partial_offset(rcrix);
      sketch->prefetch_offset_inc(off_fw);
      sketch->prefetch_offset_inc(off_rc);
      const enc_t enc_lr_fw = lshf->drop_ppos_lr(orenc64_lr);
      const enc_t enc_lr_rc = lshf->drop_ppos_lr(bp64_to_lr64(rcenc64_bp));
      sketch->prefetch_offset_enc(off_fw);
      sketch->prefetch_offset_enc(off_rc);
      uint32_t hdist_fw;
      if (sketch->scan_bucket(off_fw, enc_lr_fw, hdist_fw) && hdist_fw <= hdist_th) ++v_fw[hdist_fw];
      uint32_t hdist_rc;
      if (sketch->scan_bucket(off_rc, enc_lr_rc, hdist_rc) && hdist_rc <= hdist_th) ++(*v_rc)[hdist_rc];
    } else if (rcenc64_bp < orenc64_bp) {
      const uint32_t orrix = lshf->compute_hash(orenc64_bp);
      const uint32_t off_fw = sketch->partial_offset(orrix);
      sketch->prefetch_offset_inc(off_fw);
      const enc_t enc_lr_fw = lshf->drop_ppos_lr(orenc64_lr);
      sketch->prefetch_offset_enc(off_fw);
      uint32_t hdist_fw;
      if (sketch->scan_bucket(off_fw, enc_lr_fw, hdist_fw) && hdist_fw <= hdist_th) ++v_fw[hdist_fw];
    } else {
      const uint32_t rcrix = lshf->compute_hash(rcenc64_bp);
      const uint32_t off_rc = sketch->partial_offset(rcrix);
      sketch->prefetch_offset_inc(off_rc);
      const enc_t enc_lr_rc = lshf->drop_ppos_lr(bp64_to_lr64(rcenc64_bp));
      sketch->prefetch_offset_enc(off_rc);
      uint32_t hdist_rc;
      if (sketch->scan_bucket(off_rc, enc_lr_rc, hdist_rc) && hdist_rc <= hdist_th) ++v_fw[hdist_rc];
    }

    ++j;
    if (j >= j1) break;
    const uint64_t i = j + k - 1;
    if (__builtin_expect(SEQ_NT4_TABLE[cseq[i]] >= 4, 0)) {
      uint64_t nxt = i + 1;
      while (nxt + k <= len) {
        bool ok = true;
        for (uint64_t p = nxt; p < nxt + k; ++p) {
          if (SEQ_NT4_TABLE[cseq[p]] >= 4) {
            nxt = p + 1;
            ok = false;
            break;
          }
        }
        if (!ok) continue;
        if (nxt >= j1) return;
        j = nxt;
        compute_encoding(cseq + j, cseq + j + k, orenc64_lr, orenc64_bp);
        orenc64_bp &= mask_bp;
        orenc64_lr &= mask_lr;
        break;
      }
      if (j >= j1 || j + k > len) return;
      continue;
    }
    update_encoding(cseq + i, orenc64_lr, orenc64_bp);
    orenc64_bp &= mask_bp;
    orenc64_lr &= mask_lr;
  }
}
