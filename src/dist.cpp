#include "dist.hpp"

#include <algorithm>
#include <cmath>

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
  if (!std::isfinite(d)) return d_range;
  const auto it = std::lower_bound(th_v.begin(), th_v.end(), d);
  if (it != th_v.begin()) d_range.first = *(it - 1);
  if (it != th_v.end()) d_range.second = *it;
  return d_range;
}

std::pair<double, char> select_strand_distance(const double d_fw, const double d_rc)
{
  const bool fw_valid = std::isfinite(d_fw);
  const bool rc_valid = std::isfinite(d_rc);
  if (!fw_valid && !rc_valid) return {nanx(), '.'};
  if (!rc_valid || (fw_valid && d_fw <= d_rc)) return {d_fw, '+'};
  return {d_rc, '-'};
}

static vec<uint64_t> sample_random_coordinates(const uint64_t npos, const uint64_t nsamples, std::mt19937& rng)
{
  assert(npos >= 1);
  vec<uint64_t> starts_v;
  starts_v.reserve(nsamples);
  std::uniform_int_distribution<uint64_t> rstart(0, npos - 1);
  for (uint64_t i = 0; i < nsamples; ++i)
    starts_v.push_back(rstart(rng));
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
  , bin_shift(bin_shift)
  , tau(tau)
  , llhf(make_llhf(sketch, hdist_th))
  , k(llhf.k)
{
  canonical = sketch->is_canonical();
  bin_size = uint64_t(1) << bin_shift;
  tau_bin = std::max<uint64_t>(1, (tau + bin_size - 1) >> bin_shift);
  nwinmers = tau_bin << bin_shift;
}

void DistanceSampler::run(uint64_t sample_size, bool keep_counts, ThreadPool& pool)
{
  build(sample_size, keep_counts);
  evaluate(pool);
}

void DistanceSampler::build(uint64_t sample_size, bool keep_counts)
{
  schemes_v.clear();
  const uint64_t xtau = nwinmers + k - 1;
  uint64_t total_len = 0;
  vec<uint64_t> lenc_v(batch_v.size());
  for (size_t bix = 0; bix < batch_v.size(); ++bix) {
    if (batch_v[bix].seq.size() >= xtau) {
      total_len += batch_v[bix].seq.size();
    }
    lenc_v[bix] = total_len;
  }
  if (total_len == 0) return;

  vec<uint64_t> positions_v;
  positions_v.reserve(sample_size);
  std::uniform_int_distribution<uint64_t> rpos(0, total_len - 1);
  for (uint64_t i = 0; i < sample_size; ++i) {
    positions_v.push_back(rpos(gen));
  }
  std::sort(positions_v.begin(), positions_v.end());

  schemes_v.reserve(batch_v.size());
  size_t pidx = 0;
  for (size_t bix = 0; bix < batch_v.size() && pidx < positions_v.size(); ++bix) {
    const uint64_t L = batch_v[bix].seq.size();
    if (L < xtau) continue;
    const uint64_t rend = lenc_v[bix];
    uint64_t nsamples = 0;
    while (pidx < positions_v.size() && positions_v[pidx] < rend) {
      ++nsamples;
      ++pidx;
    }
    if (nsamples == 0) continue;
    const uint64_t enmers = L - k + 1;
    if (enmers < nwinmers) continue;
    const uint64_t nbins = (enmers + bin_size - 1) >> bin_shift;
    const uint64_t npos = (enmers - nwinmers) / bin_size + 1;
    schemes_v.emplace_back(bix, nsamples, enmers, nbins, sample_random_coordinates(npos, nsamples, gen), keep_counts);
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
        const double d = llhf.mle(agg.hist(), agg.u());
        scheme.d_v[s] = d;
        scheme.strand_v[s] = '.';
        if (scheme.keep_counts && std::isfinite(d)) {
          std::copy(agg.hist(), agg.hist() + hdist_bound + 1, scheme.hist_v.data() + s * (hdist_bound + 1));
          scheme.u_v[s] = agg.u();
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
        const double d_fw = llhf.mle(agg.hist_fw(), agg.u_fw());
        const double d_rc = llhf.mle(agg.hist_rc(), agg.u_rc());
        const auto [d, strand] = select_strand_distance(d_fw, d_rc);
        scheme.d_v[s] = d;
        scheme.strand_v[s] = strand;
        if (scheme.keep_counts && std::isfinite(d)) {
          if (strand == '-') {
            std::copy(agg.hist_rc(), agg.hist_rc() + hdist_bound + 1, scheme.hist_v.data() + s * (hdist_bound + 1));
            scheme.u_v[s] = agg.u_rc();
          } else {
            std::copy(agg.hist_fw(), agg.hist_fw() + hdist_bound + 1, scheme.hist_v.data() + s * (hdist_bound + 1));
            scheme.u_v[s] = agg.u_fw();
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
      if (std::isfinite(d)) d_v.push_back(d);
    }
  }
}

void DistanceSampler::collect_distances(vec<vec<double>>& d_vvec) const
{
  for (const auto& scheme : schemes_v) {
    auto& dst = d_vvec[scheme.bix];
    for (const double d : scheme.d_v) {
      if (std::isfinite(d)) dst.push_back(d);
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
  sampler.run(sample_size, output_samples, pool);

  if (output_samples) {
    const uint64_t nwinmers = sampler.get_nwinmers();
    const uint32_t k = sketch->get_lshf_sptr()->get_k();
    const LLH<double>& llhf = sampler.get_llhf();
    set_precision(sout, 5);

    sampler.for_each_sample([&](uint64_t bix, uint64_t enmers, uint64_t start_bin, double d, char strand) {
      const uint64_t jx = start_bin << bin_shift;
      const uint64_t jy = std::min(jx + nwinmers, enmers);
      write_tsv(sout, batch_v[bix].qid, jx + 1, jy + k - 1, strand, sketch->get_rname(), d) << '\n';
    });
    return;
  } else {
    vec<double> d_v;
    d_v.reserve(sample_size);
    sampler.collect_distances(d_v);
  }
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
