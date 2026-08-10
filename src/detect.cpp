#include "detect.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <optional>
#include <unordered_map>

#include "common.hpp"
#include "dist.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "scan.hpp"

extern uint32_t num_threads;

namespace {

  struct pair_hash
  {
    size_t operator()(const std::pair<uint64_t, uint64_t>& p) const noexcept
    {
      return std::hash<uint64_t>{}(p.first) ^ (std::hash<uint64_t>{}(p.second) << 1);
    }
  };

  // Counts of the sampled background window closest to a target distance;
  // represents a typical window's information content for the threshold screen.
  struct median_window_t
  {
    std::array<uint64_t, hdist_bound + 1> hist{};
    uint64_t u = 0;
    double dev = std::numeric_limits<double>::infinity();
    bool valid = false;
  };

  // Single pass over the samples: per query with a finite dmed_v[bix], copy the
  // counts of the window whose distance is closest to that query's median.
  vec<median_window_t> find_median_windows(const DistanceSampler& sampler, const vec<double>& dmed_v)
  {
    vec<median_window_t> out_v(dmed_v.size());
    sampler.for_each_sample_counts([&](uint64_t bix, uint64_t, uint64_t, double d, char, const uint64_t* hist, uint64_t u) {
      if (bix >= dmed_v.size() || !is_valid_distance(dmed_v[bix])) return;
      if (hist == nullptr || u == 0 || !is_valid_distance(d)) return;
      const double dev = std::abs(d - dmed_v[bix]);
      median_window_t& out = out_v[bix];
      if (dev < out.dev) {
        std::copy(hist, hist + hdist_bound + 1, out.hist.begin());
        out.u = u;
        out.dev = dev;
        out.valid = true;
      }
    });
    return out_v;
  }

} // namespace

DetectSC::DetectSC(CLI::App& sc)
{
  sc.add_option("target-path", target_path, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible)")
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("sketch-path", sketch_path, "Reference sketch file <path>")->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.add_option("-l", tau, "Minimum interval length in k-mers; also the background window length")
    ->required()
    ->check(CLI::PositiveNumber);
  sc.add_option("-b,--bin-shift", bin_shift, "Group consecutive k-mers into bins of size 2^b [0]")->check(CLI::Range(0, 16));
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match [4]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("--chisq", chisq, "Chi-square threshold for merging overlapping intervals [33.00051]")
    ->check(CLI::NonNegativeNumber);
  sc.add_option("--sample-size", sample_size, "Background windows sampled per fit [1000]")->check(CLI::PositiveNumber);
  sc.add_option("--levels", levels, "Two-sided confidence level(s) - provide exactly 1 or 4 values [0.05 0.01 0.001 0.0001]")
    ->expected(1, 4);
  sc.add_flag("--per-sequence", per_sequence, "Fit and sample the background per query sequence [per-sketch]");
  sc.add_option("--fit-quantiles", fit_quantiles, "Three central quantile probabilities for the gamma fit [0.2 0.4 0.6]")
    ->expected(3);
  sc.add_option("--verbosity", verbosity, "Statistics and report detail level (0-3) [1]")->check(CLI::Range(0, 3));
  sc.callback([&]() {
    if (levels.empty()) levels = {0.05, 0.01, 0.001, 0.0001};
    if (fit_quantiles.empty()) fit_quantiles = {0.2, 0.4, 0.6};
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

bool DetectSC::validate_configuration()
{
  bool is_invalid = false;
  if (levels.size() != 1 && levels.size() != 4) {
    is_invalid = true;
    cerr_msg("--levels requires exactly 1 or 4 levels; got ", levels.size());
  }
  for (size_t i = 0; i < levels.size(); ++i) {
    if (!std::isfinite(levels[i]) || levels[i] <= 1e-8 || levels[i] >= 0.5) {
      is_invalid = true;
      cerr_msg("--levels[", i, "] must be in (1e-8, 0.5): ", levels[i]);
    }
  }
  {
    auto slevels = levels;
    std::sort(slevels.begin(), slevels.end());
    if (const auto it = std::adjacent_find(slevels.begin(), slevels.end()); it != slevels.end()) {
      is_invalid = true;
      cerr_msg("--levels values must be unique; duplicate: ", *it);
    }
  }
  if (!std::is_sorted(fit_quantiles.begin(), fit_quantiles.end()) ||
      std::adjacent_find(fit_quantiles.begin(), fit_quantiles.end()) != fit_quantiles.end()) {
    is_invalid = true;
    cerr_msg("--fit-quantiles must be strictly increasing");
  }
  for (size_t i = 0; i < fit_quantiles.size(); ++i) {
    if (!std::isfinite(fit_quantiles[i]) || fit_quantiles[i] <= 0.0 || fit_quantiles[i] >= 1.0) {
      is_invalid = true;
      cerr_msg("--fit-quantiles[", i, "] must be in (0, 1): ", fit_quantiles[i]);
    }
  }
  if (!validate_binning(bin_shift, tau)) {
    is_invalid = true;
  }
  const uint64_t bin_size = (bin_shift <= 16) ? (uint64_t(1) << bin_shift) : 0;
  const uint64_t tau_bin = (bin_size > 0) ? ((tau + bin_size - 1) >> bin_shift) : 0;
  if (tau_bin < 2) {
    is_invalid = true;
    cerr_msg("-l must span at least two bins after binning ", "(tau=", tau, ", bin_size=", bin_size, ")");
  }
  return !is_invalid;
}

DetectSC::background_fit_t DetectSC::fit_background(const vec<double>& d_v) const
{
  background_fit_t fit;
  const auto prepared = GammaModel::prepare_samples(d_v, d_eps);
  fit.ndropped = prepared.ndropped;
  fit.nfloored = prepared.nfloored;
  fit.nsamples = prepared.x.size();
  if (prepared.x.size() < GammaModel::min_nsamples) return fit;
  if (fit.nfloored * 20 > fit.nsamples) {
    warn_msg(concat_msg(
      "more than 5% floored windows (", fit.nfloored, "/", fit.nsamples, "); low-side thresholds may be meaningless"));
  }

  GammaModel::Config cfg{};
  cfg.quantile_probs = {fit_quantiles[0], fit_quantiles[1], fit_quantiles[2]};
  fit.params = GammaModel::fit_from_samples(prepared.x, cfg, &fit.objective, &fit.niter);
  if (!GammaModel::validate_params(fit.params)) {
    fit.params = GammaModel::moments_estimate(prepared.x);
    if (!GammaModel::validate_params(fit.params)) return fit;
    // Diagnostics from the failed Nelder-Mead run do not describe this estimate.
    fit.objective = nanx();
    fit.niter = 0;
    warn_msg("Nelder-Mead gamma fit failed; falling back to the moment estimate");
  }
  if (fit.params.shape > 1e6) {
    fit.params.scale *= fit.params.shape / 1e6;
    fit.params.shape = 1e6;
    warn_msg("near-degenerate background (all samples nearly identical); clamping gamma shape while preserving its mean");
  }
  const double mass_above = 1.0 - GammaModel::cdf(d_ub, fit.params.shape, fit.params.scale);
  if (mass_above > 1e-3) {
    warn_msg(concat_msg(
      "fitted gamma puts ", mass_above * 100.0, "% of its mass above the maximum distance; tail thresholds are biased low"));
  }
  fit.ok = true;
  return fit;
}

bool DetectSC::build_lanes(const background_fit_t& fit,
                           lane_config_t& lanes,
                           const LLH<double>& llhf,
                           const uint64_t* win_hist,
                           const uint64_t win_u,
                           const double d_med) const
{
  lanes = lane_config_t{};
  vec<double> lev = levels;
  std::sort(lev.begin(), lev.end(), std::greater<double>());

  vec<double> used; // all thresholds so far, for the uniqueness requirement
  for (const double alpha : lev) {
    double t_lo = GammaModel::quantile(alpha / 2.0, fit.params.shape, fit.params.scale);
    double t_hi = GammaModel::quantile(1.0 - alpha / 2.0, fit.params.shape, fit.params.scale);
    t_lo = std::clamp(t_lo, d_eps, d_ub - d_eps);
    t_hi = std::clamp(t_hi, d_eps, d_ub - d_eps);
    if (!std::isfinite(t_lo) || !std::isfinite(t_hi) || !(t_lo < t_hi)) {
      warn_msg(concat_msg("degenerate thresholds at level ", alpha, "; dropping the level"));
      continue;
    }
    if (win_hist != nullptr && win_u > 0 && is_valid_distance(d_med)) {
      // Screen: test H0 "D = threshold" against "D = median" on a representative
      // background window (LR, chi-square(1) at this level's own alpha). A
      // threshold a typical window cannot tell apart from the median only
      // re-labels background noise, so detection at such a level is dropped.
      const double crit = GammaModel::quantile(1.0 - alpha, 0.5, 2.0); // chi-square(1) == Gamma(1/2, 2)
      const double nll_med = llhf.nll(d_med, win_hist, win_u);
      const double lr_lo = likelihood_ratio_statistic(llhf.nll(t_lo, win_hist, win_u), nll_med);
      const double lr_hi = likelihood_ratio_statistic(llhf.nll(t_hi, win_hist, win_u), nll_med);
      if (!std::isfinite(lr_lo) || !std::isfinite(lr_hi) || lr_lo < crit || lr_hi < crit) {
        warn_msg(concat_msg("level ",
                            alpha,
                            " threshold(s) indistinguishable from the background median (LR: lo=",
                            lr_lo,
                            ", hi=",
                            lr_hi,
                            ", chi2 crit=",
                            crit,
                            "); dropping the level"));
        continue;
      }
    }
    const bool dup = std::any_of(
      used.begin(), used.end(), [&](double t) { return std::abs(t - t_lo) < 1e-12 || std::abs(t - t_hi) < 1e-12; });
    if (dup) {
      warn_msg(concat_msg("duplicate thresholds at level ", alpha, " after clipping; dropping the level"));
      continue;
    }
    used.push_back(t_lo);
    used.push_back(t_hi);
    lanes.thresholds.push_back({alpha, t_lo, t_hi});
  }
  if (lanes.thresholds.empty()) return false;

  lanes.nlevels = lanes.thresholds.size();
  lanes.nlanes = 2 * lanes.nlevels;
  lanes.extrema.fill(0.5);
  lanes.alpha_v.fill(0.0);
  for (size_t j = 0; j < lanes.nlevels; ++j) {
    lanes.extrema[j] = -lanes.thresholds[j].t_hi;                // high side: d > t_hi
    lanes.extrema[lanes.nlevels + j] = lanes.thresholds[j].t_lo; // low side: d < t_lo
    lanes.alpha_v[j] = lanes.thresholds[j].alpha;
    lanes.alpha_v[lanes.nlevels + j] = lanes.thresholds[j].alpha;
  }
  return true;
}

void DetectSC::detect_queries(const sketch_sptr_t& sketch,
                              const vec<qseq_t>& batch_v,
                              const std::optional<lane_config_t>& lanes_pooled,
                              const vec<lane_config_t>* lanes_per_seq,
                              vec<detect_level_stats_t>& stats,
                              uint64_t& unmapped_iv,
                              uint64_t& unmapped_bp,
                              strstream& sout,
                              ThreadPool& pool) const
{
  const lshf_sptr_t lshf = sketch->get_lshf_sptr();
  const uint32_t k = lshf->get_k();
  const bool canonical = sketch->is_canonical();
  const str rname = sketch->get_rname();
  const uint64_t bin_size = uint64_t(1) << bin_shift;
  const uint64_t tau_bin = std::max<uint64_t>(1, (tau + bin_size - 1) >> bin_shift);

  const scan_ctx_t ctx = make_scan_ctx(*sketch, bin_shift, hdist_th);
  const LLH<double> llhf_d(k, lshf->get_h(), sketch->get_rho(), hdist_th, 0.0, false);

  const size_t nseq = batch_v.size();
  vec<strstream> out_v(nseq);
  vvec<detect_level_stats_t> stats_v(nseq);
  vec<uint64_t> unmapped_iv_v(nseq, 0);
  vec<uint64_t> unmapped_bp_v(nseq, 0);

  // Pooled mode shares one lane configuration across all sequences: build the
  // SIMD params/LLH once (LLH is const after construction; safe to share).
  std::optional<params_t<cm512_t>> pooled_params;
  llh_sptr_t<cm512_t> pooled_llhf;
  if (lanes_pooled) {
    pooled_params.emplace(lanes_pooled->extrema, hdist_th, tau, chisq, bin_shift, 0, canonical, false);
    pooled_llhf = std::make_shared<LLH<cm512_t>>(k, lshf->get_h(), sketch->get_rho(), hdist_th, lanes_pooled->extrema);
  }

  auto process = [&](const uint64_t bix) {
    const lane_config_t& lanes = lanes_pooled ? *lanes_pooled : (*lanes_per_seq)[bix];
    if (lanes.nlanes == 0) return;

    const char* cseq = batch_v[bix].seq.data();
    const uint64_t len = batch_v[bix].seq.size();
    if (len < static_cast<uint64_t>(k)) {
      warn_pmsg(batch_v[bix].qid, "skipped: sequence shorter than k-mer length ", "(len=", len, ", k=", k, ")");
      return;
    }
    const uint64_t enmers = len - k + 1;
    const uint64_t nbins = (enmers + bin_size - 1) >> bin_shift;
    if (nbins < 2) {
      warn_pmsg(batch_v[bix].qid, "skipped: fewer than two bins after binning ", "(len=", len, ", bin_size=", bin_size, ")");
      return;
    }
    const uint64_t tau_eff = std::min(tau_bin, nbins) - 1;

    // Per-sequence mode: lanes differ per query, so params/LLH are built here.
    std::optional<params_t<cm512_t>> seq_params;
    llh_sptr_t<cm512_t> seq_llhf;
    if (!pooled_llhf) {
      seq_params.emplace(lanes.extrema, hdist_th, tau, chisq, bin_shift, 0, canonical, false);
      seq_llhf = std::make_shared<LLH<cm512_t>>(k, lshf->get_h(), sketch->get_rho(), hdist_th, lanes.extrema);
    }
    const params_t<cm512_t>& params = pooled_params ? *pooled_params : *seq_params;
    const llh_sptr_t<cm512_t>& llhf = pooled_llhf ? pooled_llhf : seq_llhf;

    // Per-side ascending thresholds for distance bracketing.
    vec<double> th_hi, th_lo;
    for (const auto& th : lanes.thresholds) {
      th_hi.push_back(th.t_hi); // ascending: alpha is sorted descending
      th_lo.push_back(th.t_lo); // descending; reversed below
    }
    std::reverse(th_lo.begin(), th_lo.end());

    strstream& os = out_v[bix];
    set_precision(os, 5);
    stats_v[bix].assign(lanes.nlanes, {});
    vec<uint64_t> v_scratch;
    uint64_t u = 0, t = 0;

    // Likelihood estimates are memoized because nested lanes can emit the same
    // interval and differ only in mask and level.
    std::unordered_map<std::pair<uint64_t, uint64_t>, likelihood_estimate_t, pair_hash> estimates[2];

    auto emit = [&](DIM<cm512_t>& dim,
                    const uint64_t a_bin,
                    const uint64_t b_bin,
                    const size_t lane,
                    const bool is_rc,
                    const double d_q,
                    const double d_diff) {
      const interval_t seq_iv = get_coordinates({a_bin, b_bin}, bin_shift, enmers, k);
      auto& cache = estimates[is_rc ? 1 : 0];
      const auto key = std::make_pair(a_bin, b_bin);
      auto it = cache.find(key);
      if (it == cache.end()) {
        dim.extract_histogram(a_bin - 1, b_bin - 1, v_scratch, u, t);
        it = cache.emplace(key, compute_likelihood_estimate(llhf_d, v_scratch.data(), u, t, d_q)).first;
        if (!it->second.has_hits) {
          // Count unmapped intervals once per unique (strand, a_bin, b_bin):
          // nested lanes re-emit the same interval at each level.
          unmapped_iv_v[bix] += 1;
          unmapped_bp_v[bix] += seq_iv.b - seq_iv.a + 1;
        }
      }
      const likelihood_estimate_t& est = it->second;
      const bool unmapped = !est.has_hits;
      const double d = est.d;
      const double info = est.I, lr_bg = est.lr_bg, lr_ub = est.lr_ub;
      const bool high_side = lane < lanes.nlevels;
      const auto d_range = bracket_distance(d, high_side ? th_hi : th_lo);
      std::ostringstream d_bin;
      d_bin.flags(os.flags());
      d_bin.precision(os.precision());
      d_bin << '(' << d_range.first << ", " << d_range.second << ')';

      if (!unmapped) {
        stats_v[bix][lane].nintervals += 1;
        stats_v[bix][lane].bp_covered += seq_iv.b - seq_iv.a + 1;
      }

      const uint32_t mask = static_cast<uint32_t>(1u << lane);
      const char* side = unmapped ? "unmapped" : (high_side ? "high" : "low");
      const uint64_t L = enmers + k - 1;
      if (canonical) {
        write_tsv(os,
                  batch_v[bix].qid,
                  L,
                  seq_iv.a,
                  seq_iv.b,
                  rname,
                  d,
                  side,
                  mask,
                  lanes.alpha_v[lane],
                  d_bin.str(),
                  d_q,
                  info,
                  lr_bg,
                  lr_ub)
          << '\n';
      } else {
        write_tsv(os,
                  batch_v[bix].qid,
                  L,
                  seq_iv.a,
                  seq_iv.b,
                  report_strand(is_rc, d_diff),
                  static_cast<uint32_t>(is_rc),
                  rname,
                  d,
                  side,
                  mask,
                  lanes.alpha_v[lane],
                  d_bin.str(),
                  d_q,
                  d_diff,
                  info,
                  lr_bg,
                  lr_ub)
          << '\n';
      }
    };

    DIM<cm512_t> dim_fw(params, llhf, nbins, enmers);
    std::optional<DIM<cm512_t>> dim_rc;
    if (!canonical) dim_rc.emplace(params, llhf, nbins, enmers);
    dim_agg_t<cm512_t> agg{dim_fw, canonical ? nullptr : &*dim_rc};
    if (canonical)
      scan_mers_range<false>(ctx, cseq, 0, enmers, agg);
    else
      scan_mers_range<true>(ctx, cseq, 0, enmers, agg);

    // Uniform per-strand pipeline: prefix sums, whole-query distance, then
    // per-lane extraction. Canonical mode runs a single forward "strand".
    struct strand_t
    {
      DIM<cm512_t>* dim;
      bool is_rc;
      double d_q;
    };
    vec<strand_t> strands;
    strands.push_back({&dim_fw, false, nanx()});
    if (dim_rc) strands.push_back({&*dim_rc, true, nanx()});
    for (auto& s : strands) {
      s.dim->inclusive_scan();
      s.dim->compute_prefhistsum();
      s.dim->extrema_scan();
      vec<uint64_t> v_q;
      uint64_t u_q = 0, t_q = 0;
      s.dim->total_histogram(v_q, u_q, t_q);
      s.d_q = llhf_d.mle(v_q.data(), u_q);
    }
    const double d_diff = canonical ? nanx() : strand_diff(strands[0].d_q, strands[1].d_q);

    for (size_t ix = 0; ix < lanes.nlanes; ++ix) {
      for (auto& s : strands) {
        s.dim->extract_intervals_mx(tau_eff, 1, nbins, ix);
        s.dim->expand_intervals(params.chisq, ix);
        for (const auto& iv : s.dim->get_intervals_v(ix))
          emit(*s.dim, iv.a, iv.b + 1, ix, s.is_rc, s.d_q, d_diff);
      }
    }
  };

  pool.parallel_for(nseq, 1, process);

  unmapped_iv = 0;
  unmapped_bp = 0;
  for (size_t bix = 0; bix < nseq; ++bix) {
    if (out_v[bix].tellp() > 0) sout << out_v[bix].rdbuf();
    for (size_t lane = 0; lane < stats_v[bix].size(); ++lane) {
      stats[lane].nintervals += stats_v[bix][lane].nintervals;
      stats[lane].bp_covered += stats_v[bix][lane].bp_covered;
    }
    unmapped_iv += unmapped_iv_v[bix];
    unmapped_bp += unmapped_bp_v[bix];
  }
}

void DetectSC::report_fit(const str& rname,
                          const background_fit_t& fit,
                          const lane_config_t& lanes,
                          const double mean,
                          const double sd,
                          const uint64_t nunmapped,
                          const uint64_t nmapped) const
{
  cerr_msg("[",
           rname,
           "] background windows: n=",
           fit.nsamples,
           " (floored=",
           fit.nfloored,
           ", dropped=",
           fit.ndropped,
           ", unmapped=",
           nunmapped,
           ", hit=",
           nmapped,
           ")",
           " mean=",
           mean,
           " sd=",
           sd);
  if (verbosity >= 2) {
    cerr_msg("[",
             rname,
             "] fit: shape=",
             fit.params.shape,
             " scale=",
             fit.params.scale,
             " objective=",
             fit.objective,
             " niter=",
             fit.niter,
             " median=",
             GammaModel::quantile(0.5, fit.params.shape, fit.params.scale));
  } else {
    cerr_msg("[", rname, "] fit: shape=", fit.params.shape, " scale=", fit.params.scale);
  }
  for (const auto& th : lanes.thresholds) {
    cerr_msg("[", rname, "] level ", th.alpha, ": t_lo=", th.t_lo, " t_hi=", th.t_hi);
  }
}

void DetectSC::report_stats(const str& rname,
                            const lane_config_t& lanes,
                            const vec<detect_level_stats_t>& stats,
                            const uint64_t unmapped_iv,
                            const uint64_t unmapped_bp) const
{
  for (size_t j = 0; j < lanes.nlevels; ++j) {
    cerr_msg("[",
             rname,
             "] level ",
             lanes.thresholds[j].alpha,
             ": high: ",
             stats[j].nintervals,
             " interval(s), ",
             stats[j].bp_covered,
             " bp | low: ",
             stats[lanes.nlevels + j].nintervals,
             " interval(s), ",
             stats[lanes.nlevels + j].bp_covered,
             " bp");
  }
  cerr_msg("[", rname, "] unmapped (no k-mer hits): ", unmapped_iv, " interval(s), ", unmapped_bp, " bp");
}

void DetectSC::detect()
{
  set_precision(*output_stream, 5);

  qseq_sptr_t qs = std::make_shared<QSeq>(target_path);
  while (qs->read_next_batch()) {
  }
  const auto& batch_v = qs->get_batch_v();

  const vec<uint64_t> sketch_offsets = read_sketch_offsets(sketch_path);
  const uint32_t nsketches = static_cast<uint32_t>(sketch_offsets.size());

  ThreadPool pool(num_threads);
  cerr_msg("Processing ", nsketches, " sketches w/ ", pool.size(), " thread(s)...");

  init_thread_rng(1);

  std::ifstream sin(sketch_path, std::ifstream::binary);
  check_fstream(sin, "Cannot open sketch file for reading", sketch_path.string());

  for (uint32_t i = 0; i < nsketches; ++i) {
    sketch_sptr_t sketch = std::make_shared<Sketch>(sketch_path);
    sketch->load_from_offset(sin, sketch_offsets[i]);
    const str rname = sketch->get_rname();
    const lshf_sptr_t lshf = sketch->get_lshf_sptr();
    const uint64_t bin_size = uint64_t(1) << bin_shift;
    const uint64_t nwinmers = ((tau + bin_size - 1) >> bin_shift) << bin_shift;
    const uint64_t min_len = nwinmers + lshf->get_k() - 1;
    // Scalar LLH for the threshold screen (same configuration as in detect_queries).
    const LLH<double> llhf_d(lshf->get_k(), lshf->get_h(), sketch->get_rho(), hdist_th, 0.0, false);

    // Pass 1: sample background windows. Counts are kept so that a representative
    // window per fit can feed the threshold screen in build_lanes.
    DistanceSampler sampler(sketch, batch_v, tau, bin_shift, hdist_th);
    if (per_sequence)
      sampler.run_per_sequence(sample_size, true, pool);
    else
      sampler.run_for_all(sample_size, true, pool);
    vec<vec<double>> d_per_seq(batch_v.size());
    sampler.collect_distances(d_per_seq);
    const uint64_t nsamples = sampler.get_nsamples();
    uint64_t nmapped = 0;
    for (const auto& dv : d_per_seq)
      nmapped += dv.size();
    const uint64_t nunmapped = nsamples - nmapped;

    // Fit and build lane configurations.
    std::optional<lane_config_t> lanes_pooled;
    vec<lane_config_t> lanes_per_seq;
    const vec<lane_config_t>* lanes_seq_ptr = nullptr;
    size_t nlanes_max = 0;

    if (!per_sequence) {
      vec<double> d_all;
      d_all.reserve(nmapped);
      for (const auto& dv : d_per_seq)
        d_all.insert(d_all.end(), dv.begin(), dv.end());
      const background_fit_t fit = fit_background(d_all);
      if (!fit.ok) {
        error_exit(concat_msg("gamma fit failed for sketch ",
                              rname,
                              " (",
                              fit.nsamples,
                              " usable samples; need at least ",
                              GammaModel::min_nsamples,
                              ")"));
      }
      const double d_med = GammaModel::quantile(0.5, fit.params.shape, fit.params.scale);
      // Representative background window: closest to the fitted median across all queries.
      median_window_t win;
      for (auto& w : find_median_windows(sampler, vec<double>(batch_v.size(), d_med))) {
        if (w.valid && w.dev < win.dev) win = std::move(w);
      }
      lane_config_t lanes;
      if (!build_lanes(fit, lanes, llhf_d, win.valid ? win.hist.data() : nullptr, win.u, d_med)) {
        error_exit(concat_msg("no usable confidence levels for sketch ", rname));
      }
      if (verbosity >= 1) {
        double mean = nanx(), sd = nanx();
        if (!d_all.empty()) {
          mean = std::accumulate(d_all.begin(), d_all.end(), 0.0) / static_cast<double>(d_all.size());
          if (d_all.size() > 1) {
            double sum_sq = 0.0;
            for (const double d : d_all) {
              const double delta = d - mean;
              sum_sq += delta * delta;
            }
            sd = std::sqrt(sum_sq / static_cast<double>(d_all.size() - 1));
          } else {
            sd = 0.0;
          }
        }
        report_fit(rname, fit, lanes, mean, sd, nunmapped, nmapped);
      }
      nlanes_max = lanes.nlanes;
      lanes_pooled = std::move(lanes);
    } else {
      lanes_per_seq.resize(batch_v.size());
      // Fit all sequences first so every median is known, then gather each
      // sequence's representative (median-closest) window in one pass.
      vec<background_fit_t> fits(batch_v.size());
      vec<double> dmed_v(batch_v.size(), nanx());
      for (size_t bix = 0; bix < batch_v.size(); ++bix) {
        fits[bix] = fit_background(d_per_seq[bix]);
        if (fits[bix].ok) dmed_v[bix] = GammaModel::quantile(0.5, fits[bix].params.shape, fits[bix].params.scale);
      }
      const vec<median_window_t> win_v = find_median_windows(sampler, dmed_v);
      uint64_t nfit = 0;
      uint64_t nskipped = 0;
      for (size_t bix = 0; bix < batch_v.size(); ++bix) {
        lane_config_t lanes;
        if (!fits[bix].ok ||
            !build_lanes(
              fits[bix], lanes, llhf_d, win_v[bix].valid ? win_v[bix].hist.data() : nullptr, win_v[bix].u, dmed_v[bix]) ||
            lanes.nlevels != levels.size()) {
          if (batch_v[bix].seq.size() >= min_len) {
            warn_pmsg(batch_v[bix].qid, "gamma fit failed or degenerate; detection skipped for this sequence");
            ++nskipped;
          }
          continue;
        }
        lanes_per_seq[bix] = std::move(lanes);
        ++nfit;
      }
      if (verbosity >= 1) {
        cerr_msg("[",
                 rname,
                 "] per-sequence fits: ",
                 nfit,
                 " ok, ",
                 nskipped,
                 " skipped",
                 " (unmapped windows: ",
                 nunmapped,
                 ", hit: ",
                 nmapped,
                 ")");
      }
      lanes_seq_ptr = &lanes_per_seq;
      nlanes_max = 2 * levels.size();
    }

    // Pass 2: extract outlier intervals per lane.
    vec<detect_level_stats_t> stats(nlanes_max);
    uint64_t unmapped_iv = 0, unmapped_bp = 0;
    strstream sout;
    set_precision(sout, 5);
    detect_queries(sketch, batch_v, lanes_pooled, lanes_seq_ptr, stats, unmapped_iv, unmapped_bp, sout, pool);
    if (sout.tellp() > 0) *output_stream << sout.rdbuf();

    if (verbosity >= 1) {
      if (lanes_pooled) report_stats(rname, *lanes_pooled, stats, unmapped_iv, unmapped_bp);
      std::cerr << "\rProcessed sketch " << i + 1 << "/" << nsketches << "..." << std::flush;
      if (i + 1 == nsketches) std::cerr << std::endl;
    }
  }
}
