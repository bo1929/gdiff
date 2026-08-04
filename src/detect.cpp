#include "detect.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <numeric>
#include <optional>

#include "common.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "scan.hpp"

extern uint32_t num_threads;

namespace {

  // Per-sequence sampling plan; all randomness is resolved before the parallel
  // phase begins, so results are independent of scheduling.
  struct dplan_t
  {
    uint64_t bix = 0;
    uint64_t n_samples = 0;
    uint64_t enmers = 0;
    uint64_t nbins = 0;
    vec<uint64_t> starts;
    vec<double> d_res;
  };

  // Routes observed k-mers into the per-bin DIM accumulators.
  struct dim_agg_t
  {
    DIM<cm512_t>& fw;
    DIM<cm512_t>* rc; // null in canonical mode
    inline void operator()(uint64_t bin, uint32_t hd, bool is_rc) const { (is_rc ? *rc : fw).aggregate_mer(hd, bin); }
  };

  // Brackets a distance within one side's thresholds (ascending), like
  // QIE::get_distance_bin but per side. th_v must be sorted ascending.
  xy_t bracket_distance(double d, const vec<double>& th_v)
  {
    xy_t d_range{d_eps, d_ub};
    if (!std::isfinite(d)) return d_range;
    const auto it = std::lower_bound(th_v.begin(), th_v.end(), d);
    if (it != th_v.begin()) d_range.first = *(it - 1);
    if (it != th_v.end()) d_range.second = *it;
    return d_range;
  }

} // namespace

DetectSC::DetectSC(CLI::App& sc)
{
  sc.add_option("-q,--query-path", query_path, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible)")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("-i,--sketch-path", sketch_path, "Sketch file at <path> to query")->required()->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.add_option("-l,--min-length", tau, "Minimum interval length in k-mers; also the null window length")
    ->required()
    ->check(CLI::PositiveNumber);
  sc.add_option("-b,--bin-shift", bin_shift, "Group consecutive k-mers into bins of size 2^b [0]")->check(CLI::Range(0, 62));
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match [4]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("--chisq", chisq, "Chi-square threshold [33.00051]")->check(CLI::NonNegativeNumber);
  sc.add_option("--sample-size", sample_size, "Null windows sampled per fit [1000]")->check(CLI::PositiveNumber);
  sc.add_option("--levels", levels, "Two-sided confidence level(s) - provide exactly 1 or 4 values [0.05 0.01 0.001 0.0001]")
    ->expected(1, 4);
  sc.add_option("--fit-scope", fit_scope, "Fit the null per-sketch (pooled across queries) or per-query [per-sketch]")
    ->check(CLI::IsMember({"per-sketch", "per-query"}));
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
    if (levels[i] <= 1e-8 || levels[i] >= 0.5) {
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
  if (!std::is_sorted(fit_quantiles.begin(), fit_quantiles.end())) {
    is_invalid = true;
    cerr_msg("--fit-quantiles must be sorted ascending");
  }
  for (size_t i = 0; i < fit_quantiles.size(); ++i) {
    if (fit_quantiles[i] <= 0.0 || fit_quantiles[i] >= 1.0) {
      is_invalid = true;
      cerr_msg("--fit-quantiles[", i, "] must be in (0, 1): ", fit_quantiles[i]);
    }
  }
  if (bin_shift >= 63) {
    is_invalid = true;
    cerr_msg("--bin-shift must be less than 63; got ", bin_shift);
  }
  const uint64_t bin_size = (bin_shift < 63) ? (uint64_t(1) << bin_shift) : 0;
  if (bin_size > tau) {
    is_invalid = true;
    cerr_msg("--bin-shift gives bin_size=", bin_size, ", which exceeds --min-length=", tau);
  }
  const uint64_t tau_bin = (bin_size > 0) ? ((tau + bin_size - 1) >> bin_shift) : 0;
  if (tau_bin < 2) {
    is_invalid = true;
    cerr_msg("--min-length must span at least two bins after binning ", "(tau=", tau, ", bin_size=", bin_size, ")");
  }
  return !is_invalid;
}

gamma_fit_t DetectSC::fit_null(const vec<double>& d_v) const
{
  gamma_fit_t fit;
  vec<double> x;
  x.reserve(d_v.size());
  for (const double d : d_v) {
    if (!std::isfinite(d)) {
      ++fit.n_dropped;
    } else if (d < d_eps) {
      ++fit.n_zeros;
      x.push_back(d_eps); // exact zeros break the gamma left tail; floor them
    } else {
      x.push_back(d);
    }
  }
  fit.n_samples = x.size();
  if (x.size() < GammaModel::min_nsamples) return fit;
  if (fit.n_zeros * 20 > fit.n_samples) {
    warn_msg(concat_msg(
      "more than 5% exact-zero windows (", fit.n_zeros, "/", fit.n_samples, "); low-side thresholds may be meaningless"));
  }

  GammaModel::Config cfg{};
  cfg.quantile_probs = {fit_quantiles[0], fit_quantiles[1], fit_quantiles[2]};
  fit.params = GammaModel::fit_from_samples(x, cfg, &fit.objective, &fit.niter);
  if (!GammaModel::validate_params(fit.params)) {
    fit.params = GammaModel::moments_estimate(x);
    if (!GammaModel::validate_params(fit.params)) return fit;
    warn_msg("Nelder-Mead gamma fit failed; falling back to the moment estimate");
  }
  if (fit.params.shape > 1e6) {
    fit.params.shape = 1e6;
    warn_msg("near-degenerate null (all samples nearly identical); clamping gamma shape to 1e6");
  }
  const double mass_above = 1.0 - GammaModel::cdf(d_ub, fit.params.shape, fit.params.scale);
  if (mass_above > 1e-3) {
    warn_msg(concat_msg(
      "fitted gamma puts ", mass_above * 100.0, "% of its mass above the maximum distance; tail thresholds are biased low"));
  }
  fit.ok = true;
  return fit;
}

bool DetectSC::build_lanes(const gamma_fit_t& fit, lane_config_t& lanes) const
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
    bool dup = false;
    for (const double t : used) {
      if (std::abs(t - t_lo) < 1e-12 || std::abs(t - t_hi) < 1e-12) {
        dup = true;
        break;
      }
    }
    if (dup) {
      warn_msg(concat_msg("duplicate thresholds at level ", alpha, " after clipping; dropping the level"));
      continue;
    }
    used.push_back(t_lo);
    used.push_back(t_hi);
    lanes.thresholds.push_back({alpha, t_lo, t_hi});
  }
  if (lanes.thresholds.empty()) return false;

  lanes.n_levels = lanes.thresholds.size();
  lanes.n_lanes = 2 * lanes.n_levels;
  lanes.extrema.fill(0.5);
  lanes.alpha_v.fill(0.0);
  for (size_t j = 0; j < lanes.n_levels; ++j) {
    lanes.extrema[j] = -lanes.thresholds[j].t_hi;                 // high side: d > t_hi
    lanes.extrema[lanes.n_levels + j] = lanes.thresholds[j].t_lo; // low side: d < t_lo
    lanes.alpha_v[j] = lanes.thresholds[j].alpha;
    lanes.alpha_v[lanes.n_levels + j] = lanes.thresholds[j].alpha;
  }
  return true;
}

void DetectSC::sample_null(const sketch_sptr_t& sketch,
                           const vec<str>& seq_batch,
                           vec<vec<double>>& d_per_seq,
                           uint64_t& n_unmapped,
                           ThreadPool& pool) const
{
  const lshf_sptr_t lshf = sketch->get_lshf_sptr();
  const uint32_t k = lshf->get_k();
  const uint64_t bin_size = uint64_t(1) << bin_shift;
  // tau counts k-mers: a null window spans tau_bin bins, exactly the minimum
  // detected interval length and the same window length as map/dist.
  const uint64_t tau_bin = std::max<uint64_t>(1, (tau + bin_size - 1) >> bin_shift);
  const uint64_t win_mers = tau_bin << bin_shift; // window span in mer starts (last bin may be partial)
  const bool canonical = sketch->is_canonical();
  const uint32_t nworkers = pool.size();

  uint64_t total_len = 0;
  vec<uint64_t> cum_lens(seq_batch.size());
  for (size_t bix = 0; bix < seq_batch.size(); ++bix) {
    if (seq_batch[bix].size() >= tau + k - 1) total_len += seq_batch[bix].size();
    cum_lens[bix] = total_len;
  }
  if (total_len == 0) return;

  // File-wide uniform positions, then per-sequence sample counts (identical
  // RNG consumption to a sequential implementation).
  vec<uint64_t> positions;
  positions.reserve(sample_size);
  std::uniform_int_distribution<uint64_t> rpos(0, total_len - 1);
  for (uint64_t i = 0; i < sample_size; ++i)
    positions.push_back(rpos(gen));
  std::sort(positions.begin(), positions.end());

  vec<dplan_t> plans;
  plans.reserve(seq_batch.size());
  size_t pidx = 0;
  for (size_t bix = 0; bix < seq_batch.size() && pidx < positions.size(); ++bix) {
    const uint64_t L = seq_batch[bix].size();
    if (L < tau + k - 1) continue;
    const uint64_t seq_end = cum_lens[bix];
    uint64_t n_for_seq = 0;
    while (pidx < positions.size() && positions[pidx] < seq_end) {
      ++n_for_seq;
      ++pidx;
    }
    if (n_for_seq == 0) continue;
    const uint64_t enmers = L - k + 1;
    const uint64_t nbins = (enmers + bin_size - 1) >> bin_shift;
    if (nbins < tau_bin) continue; // defensive; cannot happen when L >= tau + k - 1
    dplan_t p;
    p.bix = bix;
    p.n_samples = n_for_seq;
    p.enmers = enmers;
    p.nbins = nbins;
    const uint64_t npos = nbins - tau_bin + 1;
    p.starts = sample_region_starts(npos, n_for_seq, gen);
    p.d_res.assign(n_for_seq, nanx());
    plans.push_back(std::move(p));
  }

  const scan_ctx_t ctx = make_scan_ctx(*sketch, bin_shift, hdist_th);
  const LLH<double> llhf(k, lshf->get_h(), sketch->get_rho(), hdist_th, 0.0, false);

  // Window scans, chunked across workers.
  vec<std::pair<uint32_t, uint64_t>> tasks; // (plan, first sample)
  for (uint32_t pi = 0; pi < plans.size(); ++pi) {
    const dplan_t& p = plans[pi];
    const uint64_t n_per = std::max<uint64_t>(4, (p.n_samples + nworkers * 4 - 1) / (nworkers * 4));
    for (uint64_t s0 = 0; s0 < p.n_samples; s0 += n_per)
      tasks.push_back({pi, s0});
  }
  pool.parallel_for(tasks.size(), 1, [&](const uint64_t ti) {
    const auto [pi, s0] = tasks[ti];
    dplan_t& p = plans[pi];
    const uint64_t s1 = std::min(s0 + std::max<uint64_t>(4, (p.n_samples + nworkers * 4 - 1) / (nworkers * 4)), p.n_samples);
    const char* cseq = seq_batch[p.bix].data();
    uint64_t v_fw[hdist_bound + 1];
    uint64_t v_rc[hdist_bound + 1];
    for (uint64_t s = s0; s < s1; ++s) {
      const uint64_t a_bin = p.starts[s];
      const uint64_t j0 = a_bin << bin_shift;
      const uint64_t j1 = std::min(j0 + win_mers, p.enmers);
      std::fill(v_fw, v_fw + hdist_bound + 1, 0);
      std::fill(v_rc, v_rc + hdist_bound + 1, 0);
      uint64_t u_fw = 0, u_rc = 0;
      window_agg_t agg{v_fw, v_rc, u_fw, u_rc, hdist_th};
      if (canonical)
        scan_mers_range<false>(ctx, cseq, j0, j1, agg);
      else
        scan_mers_range<true>(ctx, cseq, j0, j1, agg);
      // t == 0 (no k-mer hits): mle is undefined (NaN), not a plateau MLE.
      const double d_fw = llhf.mle(v_fw, u_fw);
      double d = d_fw;
      if (!canonical) {
        const double d_rc = llhf.mle(v_rc, u_rc);
        d = select_strand_distance(d_fw, d_rc).first;
      }
      p.d_res[s] = d;
    }
  });

  // Merge per-sequence results in sequence order (deterministic).
  n_unmapped = 0;
  for (const auto& p : plans) {
    for (uint64_t s = 0; s < p.n_samples; ++s) {
      const double d = p.d_res[s];
      if (!std::isfinite(d)) {
        ++n_unmapped;
        continue;
      }
      d_per_seq[p.bix].push_back(d);
    }
  }
}

void DetectSC::detect_queries(const sketch_sptr_t& sketch,
                              const vec<str>& seq_batch,
                              const vec<str>& qid_batch,
                              const lane_config_t* lanes_pooled,
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
  const str rid = sketch->get_rid();
  const uint64_t bin_size = uint64_t(1) << bin_shift;
  const uint64_t tau_bin = std::max<uint64_t>(1, (tau + bin_size - 1) >> bin_shift);

  const scan_ctx_t ctx = make_scan_ctx(*sketch, bin_shift, hdist_th);
  const LLH<double> llhf_d(k, lshf->get_h(), sketch->get_rho(), hdist_th, 0.0, false);

  const size_t nseq = seq_batch.size();
  vec<strstream> out_v(nseq);
  vvec<detect_level_stats_t> stats_v(nseq);
  vec<uint64_t> unmapped_iv_v(nseq, 0);
  vec<uint64_t> unmapped_bp_v(nseq, 0);

  auto process = [&](const uint64_t bix) {
    const lane_config_t& lanes = lanes_pooled ? *lanes_pooled : (*lanes_per_seq)[bix];
    if (lanes.n_lanes == 0) return;

    const char* cseq = seq_batch[bix].data();
    const uint64_t len = seq_batch[bix].size();
    if (len < static_cast<uint64_t>(k)) {
      warn_pmsg(qid_batch[bix], "skipped: sequence shorter than k-mer length ", "(len=", len, ", k=", k, ")");
      return;
    }
    const uint64_t enmers = len - k + 1;
    const uint64_t nbins = (enmers + bin_size - 1) >> bin_shift;
    if (nbins < 2) {
      warn_pmsg(qid_batch[bix], "skipped: fewer than two bins after binning ", "(len=", len, ", bin_size=", bin_size, ")");
      return;
    }
    const uint64_t tau_eff = std::min(tau_bin, nbins) - 1;

    params_t<cm512_t> params(lanes.n_lanes, lanes.extrema, hdist_th, tau, chisq, bin_shift, 0, canonical, false);
    llh_sptr_t<cm512_t> llhf = std::make_shared<LLH<cm512_t>>(k, lshf->get_h(), sketch->get_rho(), hdist_th, lanes.extrema);

    // Per-side ascending thresholds for distance bracketing.
    vec<double> th_hi, th_lo;
    for (const auto& th : lanes.thresholds) {
      th_hi.push_back(th.t_hi); // ascending: alpha is sorted descending
      th_lo.push_back(th.t_lo); // descending; reversed below
    }
    std::reverse(th_lo.begin(), th_lo.end());

    strstream& os = out_v[bix];
    os << std::setprecision(5);
    stats_v[bix].assign(lanes.n_lanes, {});
    vec<uint64_t> v_scratch;
    uint64_t u = 0, t = 0;

    // Window scores are memoized per interval: nested confidence lanes often
    // re-emit the same interval, and the MLE + Fisher + deviance work is
    // identical across lanes (only the mask/level columns differ).
    std::map<std::pair<uint64_t, uint64_t>, window_score_t> score_cache[2];

    // Emits one detected interval: MLE distance + record row + stats.
    auto emit = [&](DIM<cm512_t>& dim,
                    const uint64_t a_bin,
                    const uint64_t b_bin,
                    const size_t lane,
                    const bool is_rc,
                    const double d_q,
                    const double d_diff) {
      auto& cache = score_cache[is_rc ? 1 : 0];
      const auto key = std::make_pair(a_bin, b_bin);
      auto it = cache.find(key);
      if (it == cache.end()) {
        dim.extract_histogram(a_bin - 1, b_bin - 1, v_scratch, u, t);
        it = cache.emplace(key, score_window(llhf_d, v_scratch.data(), u, t, d_q)).first;
      }
      const window_score_t& ws = it->second;
      // t == 0 (no k-mer hits): distance is undefined; classify as unmapped.
      const bool unmapped = ws.unmapped;
      const double d = ws.d;
      const double info = ws.info, lr_bg = ws.lr_bg, lr_ub = ws.lr_ub;
      const interval_t seq_iv = get_coordinates({a_bin, b_bin}, bin_shift, enmers, k);
      const bool high_side = lane < lanes.n_levels;
      const auto d_range = bracket_distance(d, high_side ? th_hi : th_lo);
      std::ostringstream d_bin;
      d_bin.flags(os.flags());
      d_bin.precision(os.precision());
      d_bin << '(' << d_range.first << ", " << d_range.second << ')';

      if (unmapped) {
        unmapped_iv_v[bix] += 1;
        unmapped_bp_v[bix] += seq_iv.b - seq_iv.a + 1;
      } else {
        stats_v[bix][lane].n_intervals += 1;
        stats_v[bix][lane].bp_covered += seq_iv.b - seq_iv.a + 1;
      }

      const uint32_t mask = static_cast<uint32_t>(1u << lane);
      const char* side = unmapped ? "unmapped" : (high_side ? "high" : "low");
      const uint64_t L = enmers + k - 1;
      if (canonical) {
        write_tsv(os,
                  qid_batch[bix],
                  L,
                  seq_iv.a,
                  seq_iv.b,
                  rid,
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
                  qid_batch[bix],
                  L,
                  seq_iv.a,
                  seq_iv.b,
                  report_strand(is_rc, d_diff),
                  static_cast<uint32_t>(is_rc),
                  rid,
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
    dim_agg_t agg{dim_fw, canonical ? nullptr : &*dim_rc};
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
      // t_q == 0 (no k-mer hits at all): whole-query distance is undefined.
      s.d_q = llhf_d.mle(v_q.data(), u_q);
    }
    const double d_diff = canonical ? nanx() : strand_diff(strands[0].d_q, strands[1].d_q);

    for (size_t ix = 0; ix < lanes.n_lanes; ++ix) {
      for (auto& s : strands) {
        s.dim->extract_intervals_mx(tau_eff, 1, nbins, ix);
        s.dim->expand_intervals(params.chisq, ix);
        for (const auto& iv : s.dim->get_intervals(ix))
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
      stats[lane].n_intervals += stats_v[bix][lane].n_intervals;
      stats[lane].bp_covered += stats_v[bix][lane].bp_covered;
    }
    unmapped_iv += unmapped_iv_v[bix];
    unmapped_bp += unmapped_bp_v[bix];
  }
}

void DetectSC::report_fit(const str& rid,
                          const dist_summary_t& summary,
                          const gamma_fit_t& fit,
                          const lane_config_t& lanes,
                          const uint64_t n_unmapped) const
{
  cerr_msg("[",
           rid,
           "] null windows: n=",
           fit.n_samples,
           " (zeros=",
           fit.n_zeros,
           ", dropped=",
           fit.n_dropped,
           ", unmapped=",
           n_unmapped,
           ")",
           " mean=",
           summary.mean,
           " sd=",
           summary.sd);
  if (verbosity >= 2) {
    cerr_msg("[",
             rid,
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
    cerr_msg("[", rid, "] fit: shape=", fit.params.shape, " scale=", fit.params.scale);
  }
  for (const auto& th : lanes.thresholds) {
    cerr_msg("[", rid, "] level ", th.alpha, ": t_lo=", th.t_lo, " t_hi=", th.t_hi);
  }
}

void DetectSC::report_stats(const str& rid,
                            const lane_config_t& lanes,
                            const vec<detect_level_stats_t>& stats,
                            const uint64_t unmapped_iv,
                            const uint64_t unmapped_bp) const
{
  for (size_t j = 0; j < lanes.n_levels; ++j) {
    cerr_msg("[",
             rid,
             "] level ",
             lanes.thresholds[j].alpha,
             ": high: ",
             stats[j].n_intervals,
             " interval(s), ",
             stats[j].bp_covered,
             " bp | low: ",
             stats[lanes.n_levels + j].n_intervals,
             " interval(s), ",
             stats[lanes.n_levels + j].bp_covered,
             " bp");
  }
  cerr_msg("[", rid, "] unmapped (no k-mer hits): ", unmapped_iv, " interval(s), ", unmapped_bp, " bp");
}

void DetectSC::detect()
{
  *output_stream << std::setprecision(5);

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
    const str rid = sketch->get_rid();

    // Pass 1: sample null windows.
    vec<vec<double>> d_per_seq(seq_batch.size());
    uint64_t n_unmapped = 0;
    sample_null(sketch, seq_batch, d_per_seq, n_unmapped, pool);

    // Fit and build lane configurations.
    lane_config_t lanes_pooled;
    vec<lane_config_t> lanes_per_seq;
    const lane_config_t* lanes_ptr = nullptr;
    const vec<lane_config_t>* lanes_seq_ptr = nullptr;
    size_t n_lanes_max = 0;

    if (fit_scope == "per-sketch") {
      vec<double> d_all;
      d_all.reserve(sample_size);
      for (const auto& dv : d_per_seq)
        d_all.insert(d_all.end(), dv.begin(), dv.end());
      const gamma_fit_t fit = fit_null(d_all);
      if (!fit.ok) {
        error_exit(concat_msg("gamma fit failed for sketch ",
                              rid,
                              " (",
                              fit.n_samples,
                              " usable samples; need at least ",
                              GammaModel::min_nsamples,
                              ")"));
      }
      if (!build_lanes(fit, lanes_pooled)) {
        error_exit(concat_msg("no usable confidence levels for sketch ", rid));
      }
      if (verbosity >= 1) {
        const dist_summary_t summary = summarize_distances(std::move(d_all));
        report_fit(rid, summary, fit, lanes_pooled, n_unmapped);
      }
      lanes_ptr = &lanes_pooled;
      n_lanes_max = lanes_pooled.n_lanes;
    } else {
      lanes_per_seq.resize(seq_batch.size());
      uint64_t n_skipped = 0;
      for (size_t bix = 0; bix < seq_batch.size(); ++bix) {
        const gamma_fit_t fit = fit_null(d_per_seq[bix]);
        lane_config_t lanes;
        if (!fit.ok || !build_lanes(fit, lanes) || lanes.n_levels != levels.size()) {
          if (seq_batch[bix].size() >= tau) {
            warn_pmsg(qid_batch[bix], "gamma fit failed or degenerate; detection skipped for this sequence");
            ++n_skipped;
          }
          continue;
        }
        lanes_per_seq[bix] = std::move(lanes);
      }
      if (verbosity >= 1) {
        cerr_msg("[",
                 rid,
                 "] per-query fits: ",
                 seq_batch.size() - n_skipped,
                 " ok, ",
                 n_skipped,
                 " skipped",
                 " (unmapped windows: ",
                 n_unmapped,
                 ")");
      }
      lanes_seq_ptr = &lanes_per_seq;
      n_lanes_max = 2 * levels.size();
    }

    // Pass 2: extract outlier intervals per lane.
    vec<detect_level_stats_t> stats(n_lanes_max);
    uint64_t unmapped_iv = 0, unmapped_bp = 0;
    strstream sout;
    sout << std::setprecision(5);
    detect_queries(sketch, seq_batch, qid_batch, lanes_ptr, lanes_seq_ptr, stats, unmapped_iv, unmapped_bp, sout, pool);
    if (sout.tellp() > 0) *output_stream << sout.rdbuf();

    if (verbosity >= 1) {
      if (lanes_ptr) report_stats(rid, lanes_pooled, stats, unmapped_iv, unmapped_bp);
      std::cerr << "\rProcessed sketch " << i + 1 << "/" << nsketches << "..." << std::flush;
      if (i + 1 == nsketches) std::cerr << std::endl;
    }
  }
}
