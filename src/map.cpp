#include "map.hpp"
#include "dist.hpp"
#include "gamma.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "scan.hpp"
#include <algorithm>
#include <atomic>
#include <mutex>
#include <numeric>

extern uint32_t num_threads;

namespace {
  inline void add_to_acc(vec<uint64_t>& acc_v, uint64_t& u_acc, const vec<uint64_t>& source_v, uint64_t u)
  {
    simde__m512i s = simde_mm512_loadu_si512(acc_v.data());
    s = simde_mm512_add_epi64(s, simde_mm512_loadu_si512(source_v.data()));
    simde_mm512_storeu_si512(acc_v.data(), s);
    u_acc += u;
  }
} // namespace

template<typename T>
QIE<T>::QIE(const params_t<T>& params, const sketch_sptr_t& sketch, const lshf_sptr_t& lshf, const vec<qseq_t>& batch_v)
  : params(params)
  , sketch(sketch)
  , lshf(lshf)
  , batch_v(batch_v)
  , k(lshf->get_k())
  , h(lshf->get_h())
  , llhf(std::make_shared<LLH<T>>(k, h, sketch->get_rho(), params.hdist_th, params.dist_th))
{
  enum_only = params.enum_only;
  skip_test = (params.sample_size == 0);
  keep_hist = (!enum_only || !skip_test);
  coordinates_only = enum_only && skip_test;
  if (keep_hist) {
    // For SIMD alignment, bound is set to 8
    acc_v.assign(hdist_bound + 1, 0);
    scratch_v.assign(hdist_bound + 1, 0);
  }
}

template<typename T>
void QIE<T>::sample_background(DIM<T>& dim, const size_t first_record)
{
  if (skip_test) return;
  vec<uint64_t> lengths;
  lengths.reserve(records_v.size() - first_record);
  for (size_t ri = first_record; ri < records_v.size(); ++ri)
    lengths.push_back(records_v[ri].nbins);
  std::sort(lengths.begin(), lengths.end());
  lengths.erase(std::unique(lengths.begin(), lengths.end()), lengths.end());
  for (const uint64_t nwin_bins : lengths) {
    const auto samples = dim.sample_random_intervals(nwin_bins, bix);
    samples_v.insert(samples_v.end(), samples.begin(), samples.end());
  }
}

template<typename T>
void QIE<T>::map_sequences(std::ostream& sout, const str& rname)
{
  for (bix = 0; bix < batch_v.size(); ++bix) {
    const char* cseq = batch_v[bix].seq.data();
    const uint64_t len = batch_v[bix].seq.size();

    if (len < static_cast<uint64_t>(k)) {
      warn_pmsg(batch_v[bix].qid, "skipped: sequence shorter than k-mer length ", "(len=", len, ", k=", k, ")");
      continue;
    }

    enmers = len - k + 1;
    nbins = (enmers + params.bin_size - 1) >> params.bin_shift;
    if (nbins < 2) {
      warn_pmsg(
        batch_v[bix].qid, "skipped: fewer than two bins after binning ", "(len=", len, ", bin_size=", params.bin_size, ")");
      continue;
    }
    if (params.tau_bin > nbins) {
      warn_pmsg(batch_v[bix].qid, "minimum length is exceeded; using the full query as the effective minimum ");
    }

    const uint64_t tau_eff = std::min(params.tau_bin, nbins) - 1;
    const size_t srprev = records_v.size(); // record index before this query

    if (params.canonical) {
      DIM<T> dim(params, llhf, nbins, enmers);
      auto ctx = make_scan_ctx(*sketch, params.bin_shift, params.hdist_th);
      scan_mers_range<false>(ctx, cseq, 0, enmers, dim_agg_t<T>{dim, nullptr});
      dim.inclusive_scan();

      if (!coordinates_only) dim.compute_prefhistsum();

      double d_q = nanx();
      vec<uint64_t> v_q;
      uint64_t u_q = 0, t_q = 0;
      dim.total_histogram(v_q, u_q, t_q);
      d_q = llhf->mle(v_q.data(), u_q);

      dim.set_query_distance(d_q);
      dim.extrema_scan();

      if (enum_only) {
        extract_simple_intervals(dim, false, tau_eff, d_q);
        if (coordinates_only) continue; // skip MLE or significance
      } else {
        extract_ordered_intervals(dim, false, tau_eff, d_q);
      }
      sample_background(dim, srprev);

      for (size_t ri = srprev; ri < records_v.size(); ++ri) {
        record_t& r = records_v[ri];
        r.d_q = d_q;
        r.d_diff = nanx();
      }

      add_to_acc(acc_v, u_acc, v_q, u_q);
    } else {
      DIM<T> dim_fw(params, llhf, nbins, enmers);
      DIM<T> dim_rc(params, llhf, nbins, enmers);
      scan_mers_range<true>(
        make_scan_ctx(*sketch, params.bin_shift, params.hdist_th), cseq, 0, enmers, dim_agg_t<T>{dim_fw, &dim_rc});

      for (auto* dim : {&dim_fw, &dim_rc}) {
        dim->inclusive_scan();
        if (!coordinates_only) dim->compute_prefhistsum();
      }

      vec<uint64_t> v_q_fw, v_q_rc;
      uint64_t u_q_fw = 0, u_q_rc = 0, t_q_fw = 0, t_q_rc = 0;
      dim_fw.total_histogram(v_q_fw, u_q_fw, t_q_fw);
      dim_rc.total_histogram(v_q_rc, u_q_rc, t_q_rc);
      const double d_q_fw = llhf->mle(v_q_fw.data(), u_q_fw);
      const double d_q_rc = llhf->mle(v_q_rc.data(), u_q_rc);
      const double d_diff = strand_diff(d_q_fw, d_q_rc);

      dim_fw.set_query_distance(d_q_fw);
      dim_rc.set_query_distance(d_q_rc);

      for (auto* dim : {&dim_fw, &dim_rc}) {
        dim->extrema_scan();
      }

      if (enum_only) {
        extract_simple_intervals(dim_fw, false, tau_eff, d_q_fw);
        extract_simple_intervals(dim_rc, true, tau_eff, d_q_rc);
        if (coordinates_only) continue; // skip MLE or significance
      } else {
        extract_ordered_intervals(dim_fw, false, tau_eff, d_q_fw);
        extract_ordered_intervals(dim_rc, true, tau_eff, d_q_rc);
      }

      // The lower-distance strand is the reference: rc when the difference > 0, else fw
      const bool is_rc = d_diff > 0.0;
      sample_background(is_rc ? dim_rc : dim_fw, srprev);

      for (size_t ri = srprev; ri < records_v.size(); ++ri) {
        record_t& r = records_v[ri];
        r.d_q = r.is_rc ? d_q_rc : d_q_fw;
        r.d_diff = d_diff;
      }

      add_to_acc(acc_v, u_acc, is_rc ? v_q_rc : v_q_fw, is_rc ? u_q_rc : u_q_fw);
    }
  }

  if (keep_hist) {
    d_acc = llhf->mle(acc_v.data(), u_acc);
  }
  if (!skip_test) {
    gamma_fit_t fit;
    for (auto& r : records_v) {
      // Intact (full-query) rows have no length-matched background on their own
      // query, so the null would be degenerate; leave percentile/qvalue as NaN.
      if (r.is_intact()) continue;
      test_significance(r, samples_v, params.sample_size, batch_v[r.bix].qid, &fit);
    }
    benjamini_hochberg_correction(records_v);
  }
  report_contiguous(sout, rname);
}

template<typename T>
void QIE<T>::extract_simple_intervals(DIM<T>& dim, bool is_rc, uint64_t tau_eff, double d_q_bg)
{
  if constexpr (std::is_same_v<T, double>) {
    dim.extract_intervals_mx(tau_eff, 1, nbins);
    dim.expand_intervals(params.chisq);
    for (const auto& iv : dim.get_intervals_v(0))
      emit_record(dim, iv.a, iv.b + 1, 0, is_rc, d_q_bg);
  } else {
    for (size_t ix = 0; ix < WIDTH; ++ix) {
      dim.extract_intervals_mx(tau_eff, 1, nbins, ix);
      dim.expand_intervals(params.chisq, ix);
      for (const auto& iv : dim.get_intervals_v(ix))
        emit_record(dim, iv.a, iv.b + 1, ix, is_rc, d_q_bg);
    }
  }
}

template<typename T>
void QIE<T>::extract_ordered_intervals(DIM<T>& dim, bool is_rc, uint64_t tau_eff, double d_q_bg)
{
  const uint64_t nbins = dim.get_nbins();
  bp_v.clear();

  const auto& thrank = dim.get_thrank_v();
  if (thrank.empty()) return;

  size_t sbprev = 0;
  auto merge_from = [&](size_t start) {
    if (bp_v.size() <= start) return;
    auto nit = std::next(bp_v.begin(), start);
    auto cmp = [](const bp_t& a, const bp_t& b) { return a.a_bin < b.a_bin; };
    std::sort(nit, bp_v.end(), cmp);
    std::inplace_merge(bp_v.begin(), nit, bp_v.end(), cmp);
  };

  for (size_t ix : thrank) {
    merge_from(sbprev);

    uint64_t prev = 1;
    // Find gaps between existing intervals and extract new ones.
    for (const auto& s : bp_v) {
      if (s.a_bin >= prev + tau_eff + 1) dim.extract_intervals_mx(tau_eff, prev, s.a_bin - 1, ix);
      prev = std::max(prev, s.b_bin);
    }
    if (nbins >= prev + tau_eff) {
      dim.extract_intervals_mx(tau_eff, prev, nbins, ix);
    }
    dim.expand_intervals(params.chisq, ix);

    const auto& iv_v = dim.get_intervals_v(ix);
    const size_t nprev = bp_v.size();
    for (const auto& iv : iv_v) {
      bp_v.push_back({iv.a, iv.b + 1, ix});
    }
    sbprev = nprev;
  }
  merge_from(sbprev);

  if (bp_v.empty()) {
    if (!dim.get_has_skips()) {
      // No intervals extracted; report the full query.
      emit_record(dim, 1, nbins + 1, size_t(-1), is_rc, d_q_bg);
    } else {
      // Report one background record per maximal skip-free segment (b_bin exclusive).
      // Segments must meet the same minimum length (tau_eff + 1 bins) as extracted intervals.
      uint64_t a = 1;
      for (uint64_t x = 1; x <= nbins; ++x) {
        if (dim.is_skip(x)) {
          if (x > a + tau_eff) emit_record(dim, a, x, size_t(-1), is_rc, d_q_bg);
          a = x + 1;
        }
      }
      if (nbins >= a + tau_eff) emit_record(dim, a, nbins + 1, size_t(-1), is_rc, d_q_bg);
    }
  } else {
    for (const auto& s : bp_v) {
      emit_record(dim, s.a_bin, s.b_bin, s.ix, is_rc, d_q_bg);
    }
  }
}

template<typename T>
void QIE<T>::emit_record(DIM<T>& dim, uint64_t a_bin, uint64_t b_bin, size_t th_ix, bool is_rc, double d_q_bg)
{
  const uint64_t L = enmers + k - 1;

  likelihood_estimate_t est;
  if (!coordinates_only) {
    uint64_t u, t;
    dim.extract_histogram(a_bin - 1, b_bin - 1, scratch_v, u, t);
    est = compute_likelihood_estimate(*llhf, scratch_v.data(), u, t, d_q_bg);
    if (!est.has_hits) ++nunmapped;
  }

  const interval_t bin_iv{a_bin, b_bin};
  const interval_t seq_iv = get_coordinates(bin_iv, params.bin_shift, enmers, k);
  records_v.emplace_back(bix, L, seq_iv, bin_iv, is_rc, est.d, est.I, th_ix, est.lr_bg, est.lr_ub);
}

template<typename T>
xy_t QIE<T>::get_distance_bin(const record_t& r, const vec<double>& th_v) const
{
  xy_t d_range{d_eps, d_ub};
  if (r.th_ix != size_t(-1)) {
    const double t_i = at(llhf->get_extrema(), r.th_ix);
    const bool is_low = std::isnan(r.d_q) || t_i <= r.d_q;
    const size_t pos = static_cast<size_t>(std::lower_bound(th_v.begin(), th_v.end(), t_i) - th_v.begin());
    if (is_low) {
      // Matched low threshold: distance is in [th(i-1), th(i)).
      d_range.first = (pos > 0) ? th_v[pos - 1] : d_eps;
      d_range.second = t_i;
    } else {
      // Matched high threshold: distance is in (th(i), th(i+1)].
      d_range.first = t_i;
      d_range.second = (pos + 1 < th_v.size()) ? th_v[pos + 1] : d_ub;
    }
  } else {
    const double d = is_valid_distance(r.d) ? r.d : r.d_q;
    if (is_valid_distance(d)) d_range = bracket_distance(d, th_v);
  }
  return d_range;
}

template<typename T>
void QIE<T>::report_contiguous(std::ostream& sout, const str& rname) const
{
  vec<double> th_v(WIDTH);
  for (size_t i = 0; i < WIDTH; ++i)
    th_v[i] = at(llhf->get_extrema(), i);
  std::sort(th_v.begin(), th_v.end());

  for (const auto& r : records_v) {
    const uint8_t mask = (r.th_ix != size_t(-1)) ? static_cast<uint8_t>(1u << r.th_ix) : 0;
    const auto d_range = get_distance_bin(r, th_v);
    std::ostringstream d_bin;
    d_bin.flags(sout.flags());
    d_bin.precision(sout.precision());
    d_bin << '(' << d_range.first << ", " << d_range.second << ')';

    if (params.canonical) {
      write_tsv(sout,
                batch_v[r.bix].qid,
                r.L,
                r.seq_iv.a,
                r.seq_iv.b,
                rname,
                r.d,
                static_cast<uint32_t>(mask),
                d_bin.str(),
                r.d_q,
                d_acc,
                r.percentile,
                r.fold,
                r.qvalue,
                r.I,
                r.lr_bg,
                r.lr_ub)
        << '\n';
    } else {
      write_tsv(sout,
                batch_v[r.bix].qid,
                r.L,
                r.seq_iv.a,
                r.seq_iv.b,
                report_strand(r.is_rc, r.d_diff),
                static_cast<uint32_t>(r.is_rc),
                rname,
                r.d,
                static_cast<uint32_t>(mask),
                d_bin.str(),
                r.d_q,
                r.d_diff,
                d_acc,
                r.percentile,
                r.fold,
                r.qvalue,
                r.I,
                r.lr_bg,
                r.lr_ub)
        << '\n';
    }
  }
}

template class QIE<double>;
template class QIE<cm512_t>;

template class LLH<double>;
template class LLH<cm512_t>;

bool MapSC::validate_configuration()
{
  bool is_invalid = false;
  if (thresholds_v.size() != 1 && thresholds_v.size() != 8) {
    is_invalid = true;
    cerr_msg("--dist-th requires exactly 1 or 8 thresholds; got ", thresholds_v.size());
  }
  for (size_t i = 0; i < thresholds_v.size(); ++i) {
    if (thresholds_v[i] <= 0.0) {
      is_invalid = true;
      cerr_msg("--dist-th[", i, "] must be positive: ", thresholds_v[i]);
    }
  }
  {
    auto sorted_v = thresholds_v;
    std::sort(sorted_v.begin(), sorted_v.end());
    if (const auto it = std::adjacent_find(sorted_v.begin(), sorted_v.end()); it != sorted_v.end()) {
      is_invalid = true;
      cerr_msg("--dist-th values must be unique; duplicate: ", *it);
    }
  }
  if (hdist_th > hdist_bound) {
    is_invalid = true;
    cerr_msg("--hdist-th must be in [0, ", hdist_bound, "] with the current SIMD histogram layout; got ", hdist_th);
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

void MapSC::map()
{
  set_precision(*output_stream, 5);

  qseq_sptr_t qs = std::make_shared<QSeq>(target_path);

  // read_next_batch returns false on EOF *after* appending the final partial
  // batch, so count from the accumulated vector rather than the loop iterations.
  while (qs->read_next_batch()) {
  }
  total_qseq = qs->get_batch_v().size();

  const vec<uint64_t> sketch_offsets = read_sketch_offsets(sketch_path);
  const uint32_t nsketches = static_cast<uint32_t>(sketch_offsets.size());
  const uint32_t nthreads = std::max(1u, std::min(num_threads, nsketches));
  cerr_msg("Processing ", nsketches, " sketches w/ ", nthreads, " thread(s)...");

  std::vector<strstream> results(nsketches);
  std::vector<uint64_t> nunmapped_v(nsketches, 0);
  std::atomic<uint32_t> count_p{0};
  std::mutex cerr_mtx;

  ThreadPool pool(nthreads);
  pool.parallel_for(nsketches, 1, [&](const uint64_t i) {
    // Each sketch gets its own RNG stream: results are independent of scheduling.
    init_thread_rng(static_cast<uint32_t>(i) + 1);
    // Each task opens its own file handle so no stream sharing occurs
    std::ifstream sketch_stream(sketch_path, std::ifstream::binary);
    sketch_sptr_t sketch = std::make_shared<Sketch>(sketch_path);
    sketch->load_from_offset(sketch_stream, sketch_offsets[i]);
    sketch_stream.close();
    bool canonical = sketch->is_canonical();

    strstream sout;
    set_precision(sout, 5);
    if (thresholds_v.size() == 1) {
      params_t<double> params(thresholds_v.front(), hdist_th, tau, chisq, bin_shift, sample_size, canonical, enum_only);
      QIE<double> qie(params, sketch, sketch->get_lshf_sptr(), qs->get_batch_v());
      qie.map_sequences(sout, sketch->get_rname());
      nunmapped_v[i] = qie.get_nunmapped();
    } else {
      params_t<cm512_t> params({0}, hdist_th, tau, chisq, bin_shift, sample_size, canonical, enum_only);
      std::copy(thresholds_v.begin(), thresholds_v.end(), params.dist_th.begin());
      QIE<cm512_t> qie(params, sketch, sketch->get_lshf_sptr(), qs->get_batch_v());
      qie.map_sequences(sout, sketch->get_rname());
      nunmapped_v[i] = qie.get_nunmapped();
    }

    results[i] = std::move(sout);

    uint32_t num_p = count_p.fetch_add(1, std::memory_order_relaxed) + 1;
    {
      std::lock_guard<std::mutex> lock(cerr_mtx);
      std::cerr << "\rProcessed sketch " << num_p << "/" << nsketches << "..." << std::flush;
      if (num_p == nsketches) std::cerr << std::endl;
    }
  });

  const uint64_t nunmapped = std::accumulate(nunmapped_v.begin(), nunmapped_v.end(), uint64_t(0));
  if (nunmapped > 0) {
    cerr_msg("Unmapped intervals (no k-mer hits): ", nunmapped, " (distance reported as NA)");
  }

  for (uint32_t i = 0; i < nsketches; ++i) {
    if (results[i].tellp() > 0) *(output_stream) << results[i].rdbuf();
  }
}

MapSC::MapSC(CLI::App& sc)
{
  sc.add_option("target-path", target_path, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible)")
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("sketch-path", sketch_path, "Reference sketch file <path>")->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match [4]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("--chisq", chisq, "Chi-square threshold [33.00051]")->check(CLI::NonNegativeNumber);
  sc.add_option("-d,--dist-th", thresholds_v, "Distance threshold(s) - provide exactly 1 or 8 values")
    ->required()
    ->expected(1, 8);
  sc.add_option("-l", tau, "Minimum interval length in k-mers")->required()->check(CLI::PositiveNumber);
  sc.add_option("-b,--bin-shift", bin_shift, "Group consecutive k-mers into bins of size 2^b [0]")->check(CLI::Range(0, 16));
  sc.add_flag("--enum-only,!--no-enum-only", enum_only, "Enumerate intervals without MLE distance estimation [false]");
  sc.add_option("--sample-size", sample_size, "Samples for significance test (0: skip) [200]")->check(CLI::NonNegativeNumber);
  sc.callback([&]() {
    if (!validate_configuration()) {
      error_exit("Invalid configuration!");
    }
    if (!output_path.empty()) {
      output_file.open(output_path);
      output_stream = &output_file;
    }
  });
}
