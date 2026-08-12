#include "dist.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_set>

#include "common.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "scan.hpp"

extern uint32_t num_threads;

static LLH<double> make_llhf(const sketch_sptr_t& sketch, uint32_t hdist_th)
{
  const auto& lshf = *sketch->get_lshf_sptr();
  return {lshf.get_k(), lshf.get_h(), sketch->get_rho(), hdist_th, 0.0, false};
}

double linear_quantile(const vec<double>& v, const double p)
{
  if (v.empty()) return nanx();
  if (v.size() == 1) return v.front();
  const double ix = p * static_cast<double>(v.size() - 1);
  const size_t lo = static_cast<size_t>(std::floor(ix));
  const size_t hi = static_cast<size_t>(std::ceil(ix));
  return v[lo] + (ix - static_cast<double>(lo)) * (v[hi] - v[lo]);
}

xy_t bracket_distance(const double d, const vec<double>& th_v)
{
  xy_t d_range{d_eps, d_ub};
  if (!is_valid_distance(d)) return d_range;
  const auto it = std::lower_bound(th_v.begin(), th_v.end(), d);
  if (it != th_v.begin()) d_range.first = *(it - 1);
  if (it != th_v.end()) d_range.second = *it;
  return d_range;
}

std::pair<double, char> select_strand_distance(const double d_fw, const double d_rc)
{
  const bool fw_valid = is_valid_distance(d_fw);
  const bool rc_valid = is_valid_distance(d_rc);
  if (!fw_valid && !rc_valid) return {nanx(), '.'};
  if (!rc_valid || (fw_valid && d_fw <= d_rc)) return {d_fw, '+'};
  return {d_rc, '-'};
}

static vec<uint64_t> sample_random_coordinates(const uint64_t npos, const uint64_t nsamples, std::mt19937& rng)
{
  assert(npos >= 1);
  const uint64_t n = std::min(npos, nsamples);
  // The full-shuffle draw is O(npos) memory, which can dominate for long
  // queries. For sparse draws, rejection-sample distinct positions instead:
  // O(n) memory and O(n) expected time while n << npos.
  if (n < npos / 2) {
    vec<uint64_t> starts_v;
    starts_v.reserve(n);
    std::unordered_set<uint64_t> seen;
    seen.reserve(n);
    std::uniform_int_distribution<uint64_t> pick(0, npos - 1);
    while (starts_v.size() < n) {
      const uint64_t x = pick(rng);
      if (seen.insert(x).second) starts_v.push_back(x);
    }
    return starts_v;
  }
  vec<uint64_t> starts_v(npos);
  std::iota(starts_v.begin(), starts_v.end(), uint64_t(0));
  std::shuffle(starts_v.begin(), starts_v.end(), rng);
  starts_v.resize(n);
  return starts_v;
}

DistanceSampler::DistanceSampler(const sketch_sptr_t& sketch,
                                 const vec<qseq_t>& batch_v,
                                 uint64_t tau,
                                 uint64_t bin_shift,
                                 uint32_t hdist_th)
  : sketch(sketch)
  , batch_v(batch_v)
  , hdist_th(hdist_th)
  , tau(tau)
  , bin_shift(bin_shift)
  , llhf(make_llhf(sketch, hdist_th))
  , k(llhf.k)
{
  canonical = sketch->is_canonical();
  bin_size = uint64_t(1) << bin_shift;
  tau_bin = std::max<uint64_t>(1, (this->tau + bin_size - 1) >> bin_shift);
  nwinmers = tau_bin << bin_shift;
}

void DistanceSampler::run_for_all(uint64_t sample_size, bool keep_counts, ThreadPool& pool)
{
  build_for_all(sample_size, keep_counts);
  evaluate(pool);
}

void DistanceSampler::run_per_sequence(uint64_t sample_size, bool keep_counts, ThreadPool& pool)
{
  build_per_sequence(sample_size, keep_counts);
  evaluate(pool);
}

void DistanceSampler::build_for_all(uint64_t sample_size, bool keep_counts)
{
  schemes_v.clear();
  const uint64_t xtau = nwinmers + k - 1;
  uint64_t total_npos = 0;
  vec<uint64_t> lenc_v(batch_v.size());
  for (size_t bix = 0; bix < batch_v.size(); ++bix) {
    const uint64_t L = batch_v[bix].seq.size();
    if (L >= xtau) {
      const uint64_t enmers = L - k + 1;
      total_npos += (enmers - nwinmers) / bin_size + 1;
    }
    lenc_v[bix] = total_npos;
  }
  if (total_npos == 0) return;

  vec<uint64_t> positions_v = sample_random_coordinates(total_npos, sample_size, gen);
  std::sort(positions_v.begin(), positions_v.end());

  schemes_v.reserve(batch_v.size());
  size_t pidx = 0;
  for (size_t bix = 0; bix < batch_v.size() && pidx < positions_v.size(); ++bix) {
    const uint64_t L = batch_v[bix].seq.size();
    if (L < xtau) continue;
    const uint64_t rend = lenc_v[bix];
    const uint64_t roff = bix == 0 ? 0 : lenc_v[bix - 1];
    vec<uint64_t> starts_v;
    while (pidx < positions_v.size() && positions_v[pidx] < rend) {
      starts_v.push_back(positions_v[pidx] - roff);
      ++pidx;
    }
    if (starts_v.empty()) continue;
    const uint64_t enmers = L - k + 1;
    const uint64_t nbins = (enmers + bin_size - 1) >> bin_shift;
    schemes_v.emplace_back(bix, starts_v.size(), enmers, nbins, std::move(starts_v), keep_counts);
  }
}

void DistanceSampler::build_per_sequence(const uint64_t sample_size, const bool keep_counts)
{
  schemes_v.clear();
  const uint64_t xtau = nwinmers + k - 1;
  schemes_v.reserve(batch_v.size());
  for (size_t bix = 0; bix < batch_v.size(); ++bix) {
    const uint64_t L = batch_v[bix].seq.size();
    if (L < xtau) continue;
    const uint64_t enmers = L - k + 1;
    const uint64_t nbins = (enmers + bin_size - 1) >> bin_shift;
    const uint64_t npos = (enmers - nwinmers) / bin_size + 1;
    vec<uint64_t> starts_v = sample_random_coordinates(npos, sample_size, gen);
    schemes_v.emplace_back(bix, starts_v.size(), enmers, nbins, std::move(starts_v), keep_counts);
  }
}

void DistanceSampler::evaluate(ThreadPool& pool)
{
  const uint32_t nworkers = pool.size();
  const scan_ctx_t ctx = make_scan_ctx(*sketch, bin_shift, hdist_th);

  struct task_t
  {
    uint32_t six;
    uint64_t a, b;
  };
  vec<task_t> tasks_v;
  for (uint32_t six = 0; six < schemes_v.size(); ++six) {
    const scheme_t& scheme = schemes_v[six];
    const uint64_t nper = std::max<uint64_t>(4, (scheme.nsamples + nworkers * 4 - 1) / (nworkers * 4));
    for (uint64_t s0 = 0; s0 < scheme.nsamples; s0 += nper)
      tasks_v.push_back({six, s0, std::min(s0 + nper, scheme.nsamples)});
  }

  pool.parallel_for(tasks_v.size(), 1, [&](const uint64_t ti) {
    const task_t& t = tasks_v[ti];
    scheme_t& scheme = schemes_v[t.six];
    const char* cseq = batch_v[scheme.bix].seq.data();
    if (canonical) {
      window_counts_t agg(hdist_th);
      for (uint64_t s = t.a; s < t.b; ++s) {
        const uint64_t a_bin = scheme.starts_v[s];
        const uint64_t jx = a_bin << bin_shift;
        const uint64_t jy = std::min(jx + nwinmers, scheme.enmers);
        agg.clear();
        scan_mers_range<false>(ctx, cseq, jx, jy, agg);
        const double d = llhf.mle(agg.hist(), agg.u);
        scheme.d_v[s] = d;
        scheme.strand_v[s] = '.';
        if (scheme.keep_counts && is_valid_distance(d)) {
          std::copy(agg.hist(), agg.hist() + hdist_bound + 1, scheme.hist_v.data() + s * (hdist_bound + 1));
          scheme.u_v[s] = agg.u;
        }
      }
    } else {
      swindow_counts_t agg(hdist_th);
      for (uint64_t s = t.a; s < t.b; ++s) {
        const uint64_t a_bin = scheme.starts_v[s];
        const uint64_t jx = a_bin << bin_shift;
        const uint64_t jy = std::min(jx + nwinmers, scheme.enmers);
        agg.clear();
        scan_mers_range<true>(ctx, cseq, jx, jy, agg);
        const double d_fw = llhf.mle(agg.hist_fw(), agg.u_fw);
        const double d_rc = llhf.mle(agg.hist_rc(), agg.u_rc);
        const auto [d, strand] = select_strand_distance(d_fw, d_rc);
        scheme.d_v[s] = d;
        scheme.strand_v[s] = strand;
        if (scheme.keep_counts && is_valid_distance(d)) {
          if (strand == '-') {
            std::copy(agg.hist_rc(), agg.hist_rc() + hdist_bound + 1, scheme.hist_v.data() + s * (hdist_bound + 1));
            scheme.u_v[s] = agg.u_rc;
          } else {
            std::copy(agg.hist_fw(), agg.hist_fw() + hdist_bound + 1, scheme.hist_v.data() + s * (hdist_bound + 1));
            scheme.u_v[s] = agg.u_fw;
          }
        }
      }
    }
  });
}

void DistanceSampler::collect_distances(vec<double>& d_v) const
{
  for (const auto& scheme : schemes_v) {
    for (const double d : scheme.d_v) {
      if (is_valid_distance(d)) d_v.push_back(d);
    }
  }
}

void DistanceSampler::collect_distances(vec<vec<double>>& d_vvec) const
{
  for (const auto& scheme : schemes_v) {
    auto& dst = d_vvec[scheme.bix];
    for (const double d : scheme.d_v) {
      if (is_valid_distance(d)) dst.push_back(d);
    }
  }
}

uint64_t DistanceSampler::get_nsamples() const
{
  uint64_t n = 0;
  for (const auto& scheme : schemes_v)
    n += scheme.nsamples;
  return n;
}

void DistSC::sample_distances(const sketch_sptr_t& sketch, const vec<qseq_t>& batch_v, strstream& sout, ThreadPool& pool)
{
  DistanceSampler sampler(sketch, batch_v, tau, bin_shift, hdist_th);
  sampler.run_for_all(sample_size, true, pool);

  const uint64_t nwinmers = sampler.get_nwinmers();
  const uint32_t k = sketch->get_lshf_sptr()->get_k();
  const LLH<double>& llhf = sampler.get_llhf();

  vec<double> d_v;
  d_v.reserve(sample_size);
  sampler.collect_distances(d_v);
  std::sort(d_v.begin(), d_v.end());
  const double d_median = linear_quantile(d_v, 0.5);

  set_precision(sout, output_samples ? 5 : 8);

  if (output_samples) {
    sampler.for_each_sample_counts(
      [&](uint64_t bix, uint64_t enmers, uint64_t start_bin, double d, char strand, const uint64_t* hist, uint64_t u) {
        double lr_bg = nanx();
        double lr_ub = nanx();
        if (hist && is_valid_distance(d)) {
          uint64_t n_total = u;
          for (uint32_t di = 0; di <= llhf.hdist_th; ++di)
            n_total += hist[di];
          lr_ub = compute_lr_ub(llhf, d, n_total);
          if (is_valid_distance(d_median))
            lr_bg = likelihood_ratio_statistic(llhf.nll(d_median, hist, u), llhf.nll(d, hist, u));
        }
        const uint64_t jx = start_bin << bin_shift;
        const uint64_t jy = std::min(jx + nwinmers, enmers);
        write_tsv(sout, batch_v[bix].qid, jx + 1, jy + k - 1, strand, sketch->get_rname(), d, lr_bg, lr_ub) << '\n';
      });
    return;
  }

  const uint64_t n = d_v.size();
  const uint64_t nwinu = sampler.get_nsamples() - n;
  if (nwinu > 0) {
    cerr_msg(
      "[", sketch->get_rname(), "] unmapped sampled windows (no k-mer hits): ", nwinu, " (excluded from the summary)");
  }

  write_tsv(sout, target_path, sketch->get_rname(), n, d_median) << '\n';
}

void DistSC::dist()
{
  qseq_sptr_t qs = std::make_shared<QSeq>(target_path);
  while (qs->read_next_batch()) {
  }

  const vec<uint64_t> sketch_offsets = read_sketch_offsets(sketch_path);
  const uint32_t nsketches = static_cast<uint32_t>(sketch_offsets.size());

  ThreadPool pool(num_threads);
  cerr_msg("Processing ", nsketches, " sketches w/ ", pool.size(), " thread(s)...");

  init_thread_rng(1);

  const auto& batch_v = qs->get_batch_v();

  std::ifstream sin(sketch_path, std::ifstream::binary);
  check_fstream(sin, "Cannot open sketch file for reading", sketch_path.string());

  for (uint32_t i = 0; i < nsketches; ++i) {
    sketch_sptr_t sketch = std::make_shared<Sketch>(sketch_path);
    sketch->load_from_offset(sin, sketch_offsets[i]);

    strstream sout;
    sample_distances(sketch, batch_v, sout, pool);
    if (sout.tellp() > 0) *output_stream << sout.rdbuf();

    std::cerr << "\rProcessed sketch " << i + 1 << "/" << nsketches << "..." << std::flush;
    if (i + 1 == nsketches) std::cerr << std::endl;
  }
}

DistSC::DistSC(CLI::App& sc)
{
  sc.add_option("target-path", target_path, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible)")
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("sketch-path", sketch_path, "Reference sketch file <path>")->check(CLI::ExistingFile);
  sc.add_option("-l", tau, "Length of sampled query regions in k-mers")->required()->check(CLI::PositiveNumber);
  sc.add_option("-b,--bin-shift", bin_shift, "Group consecutive k-mers into bins of size 2^b [0]")->check(CLI::Range(0, 16));
  sc.add_option("--sample-size", sample_size, "Regions sampled across the query file per reference [200]")
    ->check(CLI::PositiveNumber);
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match [4]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.add_flag("--output-samples", output_samples, "Write per-sample output instead of per-reference summary");
  sc.callback([&]() {
    if (!validate_configuration()) {
      error_exit("Invalid configuration!");
    }
    if (!output_path.empty()) {
      output_file.open(output_path);
      check_fstream(output_file, "Cannot open output file", output_path.string());
      output_stream = &output_file;
    }
  });
}

bool validate_binning(const uint64_t bin_shift, const uint64_t tau)
{
  bool is_invalid = false;
  if (bin_shift > 16) {
    is_invalid = true;
    cerr_msg("--bin-shift must be less than or equal to 16; got ", bin_shift);
  }
  const uint64_t bin_size = (bin_shift <= 16) ? (uint64_t(1) << bin_shift) : 0;
  if (bin_size > tau) {
    is_invalid = true;
    cerr_msg("--bin-shift gives bin_size=", bin_size, ", which exceeds -l=", tau);
  }
  return !is_invalid;
}

bool DistSC::validate_configuration()
{
  bool is_invalid = false;
  if (!validate_binning(bin_shift, tau)) {
    is_invalid = true;
  }
  return !is_invalid;
}
