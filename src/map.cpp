#include "map.hpp"
#include "random.hpp"
#include <boost/math/tools/minima.hpp>

namespace {
  template<typename T>
  inline double at(T v, const size_t idx)
  {
    if constexpr (std::is_same_v<T, double>) {
      return v;
    } else {
      return v[idx];
    }
  }

  template<typename T>
  std::pair<double, double> mle(const llh_sptr_t<T>& llhf, const vec<uint64_t>& v, uint64_t u, uint64_t t)
  {
    llhf->set_counts(v.data(), u);
    auto f = [&](const double& D) { return (*llhf)(D); };
    const double ub = (t == 0) ? 0.75 + eps : 0.5;
    return boost::math::tools::brent_find_minima(f, eps, ub, 24);
  }
} // namespace

template<typename T>
QIE<T>::QIE(const params_t<T>& params,
            const sketch_sptr_t& sketch,
            const lshf_sptr_t& lshf,
            const vec<str>& seq_batch,
            const vec<str>& qid_batch)
  : params(params)
  , sketch(sketch)
  , lshf(lshf)
  , seq_batch(seq_batch)
  , qid_batch(qid_batch)
  , batch_size(seq_batch.size())
  , k(lshf->get_k())
  , h(lshf->get_h())
  , m(lshf->get_m())
{
  llhf = std::make_shared<LLH<T>>(k, h, sketch->get_rho(), params.hdist_th, params.dist_th);
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  mask_bp = u64m >> ((32 - k) * 2);
  if (!params.enum_only) {
    v_acc.assign(hdist_bound + 1, 0); // For SIMD alignment; set to 8.
    stride_len = std::ceil((params.tau_bin * subsample_factor) + 1);
  }
}

template<typename T>
void QIE<T>::map_sequences(std::ostream& sout, const str& rid)
{
  for (bix = 0; bix < batch_size; ++bix) {
    const char* cseq = seq_batch[bix].data();
    const uint64_t len = seq_batch[bix].size();
    onmers = 0;
    if (len < static_cast<uint64_t>(k)) {
      // TODO: should we have a more verbose mode where this reported?
      continue;
    }
    enmers = len - k + 1;
    nbins = (enmers + params.bin_size - 1) >> params.bin_shift;
    if (nbins < 2) {
      warn_msg("The bin size is too high for the query sequence length: " + std::to_string(len));
      continue;
    }
    if (params.tau_bin > nbins) {
      warn_msg("The minimum length threshold is too high for a query sequence of length: " + std::to_string(len));
      continue;
    }

    DIM<T> dim_fw(params, llhf, nbins, enmers);
    DIM<T> dim_rc(params, llhf, nbins, enmers);
    search_mers(cseq, len, dim_fw, dim_rc);

    // Extract intervals for all thresholds
    const uint64_t tau_eff = std::min(params.tau_bin, nbins) - 1;
    for (auto* dim : {&dim_fw, &dim_rc}) {
      dim->inclusive_scan();
      dim->extrema_scan();
      // dim->release_accumulators();
      if constexpr (std::is_same_v<T, double>) {
        /* dim->extract_intervals_sx(tau_eff); */
        dim->extract_intervals_mx(tau_eff);
        dim->expand_intervals(params.chisq);
      } else {
        for (size_t i = 0; i < WIDTH; ++i) {
          /* dim->extract_intervals_sx(tau_eff, i); */
          dim->extract_intervals_mx(tau_eff, i);
          dim->expand_intervals(params.chisq, i);
        }
      }
    }

    if (params.enum_only) {
      if constexpr (std::is_same_v<T, double>) {
        report_intervals(sout, rid, dim_fw, false);
        report_intervals(sout, rid, dim_rc, true);
      } else {
        for (size_t i = 0; i < WIDTH; ++i) {
          report_intervals(sout, rid, dim_fw, false, i);
          report_intervals(sout, rid, dim_rc, true, i);
        }
      }
    } else {
      const auto [pos_bv, neg_bv] = llhf->get_sign_bv();
      for (auto* dim : {&dim_fw, &dim_rc}) {
        dim->compute_prefhistsum();
        dim->map_contiguous_segments(pos_bv, '<');
        dim->map_contiguous_segments(neg_bv, '>');
      }

      // Extract per-query histograms, compute per-query MLE, accumulate into genome-wide histogram
      vec<uint64_t> v_q_fw, v_q_rc;
      uint64_t u_q_fw = 0, u_q_rc = 0, t_q_fw = 0, t_q_rc = 0;
      dim_fw.extract_histogram(0, nbins, v_q_fw, u_q_fw, t_q_fw);
      dim_rc.extract_histogram(0, nbins, v_q_rc, u_q_rc, t_q_rc);
      const double d_q_fw = compute_mle_dist(v_q_fw, u_q_fw, t_q_fw);
      const double d_q_rc = compute_mle_dist(v_q_rc, u_q_rc, t_q_rc);

      // Determine reference strand as the one with lower distance
      const bool is_r = std::isnan(d_q_rc) || (!std::isnan(d_q_fw) && d_q_fw <= d_q_rc);

      dim_fw.set_rstrand(is_r);
      dim_rc.set_rstrand(!is_r);

      collect_segments(dim_fw, d_q_fw, false);
      collect_segments(dim_rc, d_q_rc, true);
      save_qstride(is_r ? dim_fw : dim_rc);
      const auto& v_q = is_r ? v_q_fw : v_q_rc;
      const uint64_t u_q = is_r ? u_q_fw : u_q_rc;
      simde__m512i v_sacc = simde_mm512_loadu_si512(v_acc.data());
      v_sacc = simde_mm512_add_epi64(v_sacc, simde_mm512_loadu_si512(v_q.data()));
      simde_mm512_storeu_si512(v_acc.data(), v_sacc);
      u_acc += u_q;
    }
  }

  if (!params.enum_only) {
    uint64_t t_acc = 0;
    for (auto c : v_acc)
      t_acc += c;
    d_acc = compute_mle_dist(v_acc, u_acc, t_acc);
    fit_gamma_significance();
    emit_segments(sout, rid);
  }
}

template<typename T>
DIM<T>::DIM(const params_t<T>& params, const llh_sptr_t<T>& llhf, uint64_t nbins, uint64_t nmers)
  : params(params)
  , llhf(llhf)
  , nbins(nbins)
  , nmers(nmers)
{
  fdc_v.resize(nbins); // Alternative?: fdc_v.reserve(nbins);
  sdc_v.resize(nbins); // Alternative?: sdc_v.reserve(nbins);
  if (!params.enum_only) {
    // Row 0 is zeros (sentinel), same layout as fdps_v/sdps_v.
    // Later, aggregate_mer() accumulates into rows 1..nbins; compute_prefhistsum() converts in-place.
    hdisthist_v.assign((nbins + 1) * (params.hdist_th + 1), 0);
  }
}

template<typename T>
void QIE<T>::emit_segments(std::ostream& sout, const str& rid) const
{ // TODO: Revisit?
  for (const auto& r : records_v) {
    write_tsv(sout,
              qid_batch[r.bix],
              r.L,
              r.a,
              r.b,
              r.strand,
              // r.rstrand ? "T" : "F", // TODO: Enable for visualization.
              rid,
              r.d,
              static_cast<uint32_t>(r.mask),
              r.sign,
              r.d_q,
              d_acc,
              r.percentile,
              r.fold)
      << '\n';
  }
}

template<typename T>
void QIE<T>::report_intervals(std::ostream& sout, const str& rid, DIM<T>& dim, bool rc, size_t idx)
{ // TODO: Revisit?
  const str strand = rc ? "-" : "+";
  const double dist_th = at(params.dist_th, idx);
  const uint64_t nbins = dim.get_nbins();
  const uint64_t L = enmers + k - 1;
  for (uint64_t i = 0;; ++i) {
    const interval_t ab = dim.get_interval(i, idx);
    if (ab.first >= nbins) break;
    const uint64_t a = ab.first << params.bin_shift;
    const uint64_t b = std::min(((ab.second + 1) << params.bin_shift), enmers) + k - 1;
    assert(a < b);
    // write_tsv(sout, qid_batch[bix], L, a, b, strand, rid, dist_th) << '\n';
    write_tsv(sout, qid_batch[bix], L, a, b - 1, strand, rid, dist_th) << '\n';
  }
}

template<typename T>
void QIE<T>::fit_gamma_significance()
{ // TODO: Debug.
  const uint64_t G = stride_len;

  // Find the largest bin count (G aligned) across all stored queries.
  uint64_t max_nbins = 0;
  for (const auto& qs : qstrides_v) {
    max_nbins = std::max(max_nbins, (qs.nbins / G) * G);
  }
  if (max_nbins < G) {
    warn_msg("All queries are too short for significance testing with strides of length " + std::to_string(G));
    return;
  }

  // Build the full geometric grid of lengths up to the maximum number of bins.
  vec<uint64_t> L_v;
  for (uint64_t L = G; L <= max_nbins;) {
    L_v.push_back(L);
    uint64_t step = std::ceil((L + G) * (grid_growth - 1.0) / G) * G;
    L += step;
  }

  assert(!L_v.empty());

  // Map each record to its nearest grid index.
  vec<size_t> idx_v(records_v.size());
  vec<bool> seen_v(L_v.size(), false);
  for (size_t ridx = 0; ridx < records_v.size(); ++ridx) {
    auto it = std::lower_bound(L_v.begin(), L_v.end(), records_v[ridx].nbins_s);
    if (it == L_v.end()) {
      --it; // Round to largest grid length
    } else if (it != L_v.begin()) {
      // Round to whichever neighbor is closer
      auto prev = std::prev(it);
      if ((records_v[ridx].nbins_s - *prev) < (*it - records_v[ridx].nbins_s)) {
        it = prev;
      }
    }
    const size_t grid_idx = static_cast<size_t>(it - L_v.begin());
    idx_v[ridx] = grid_idx;
    seen_v[grid_idx] = true;
  }

  // Sample only for seen grid lengths and fit a global Gamma for each.
  vec<vec<snull_t>> samples_vvec(L_v.size());
  vec<GammaModel::params_t> gp_v(L_v.size());
  vec<bool> valid_v(L_v.size(), false);

  for (size_t grid_idx = 0; grid_idx < L_v.size(); ++grid_idx) {
    if (!seen_v[grid_idx]) continue;
    sample_distances(L_v[grid_idx], samples_vvec[grid_idx]);

    if (samples_vvec[grid_idx].size() >= GammaModel::min_nsamples) {
      std::sort(samples_vvec[grid_idx].begin(), samples_vvec[grid_idx].end(), [](const snull_t& a, const snull_t& b) {
        return a.d < b.d;
      });
      if (!params.ecdf_test) {
        vec<double> d_v(samples_vvec[grid_idx].size());
        for (size_t i = 0; i < d_v.size(); ++i)
          d_v[i] = samples_vvec[grid_idx][i].d;
        gp_v[grid_idx] = GammaModel::fit(d_v);
      }
      valid_v[grid_idx] = true;
    } else {
      warn_msg("Failed to sample sufficient number of distances for the length " + std::to_string(L_v[grid_idx]));
    }
  }

  // Score each record, filtering out null samples that overlap with the record itself.
  using boost::math::gamma_distribution;
  for (size_t ridx = 0; ridx < records_v.size(); ++ridx) {
    auto& r = records_v[ridx];
    const size_t grid_idx = idx_v[ridx];
    if (!valid_v[grid_idx]) {
      warn_msg("No valid null distribution for this segment length; skipping p-value");
      r.percentile = std::numeric_limits<double>::quiet_NaN();
      r.fold = std::numeric_limits<double>::quiet_NaN();
      continue;
    }

    // Collect non-overlapping null distances for this record
    const size_t qidx = r.bix;
    const uint64_t bin_a = (r.a - 1) >> params.bin_shift;
    const uint64_t bin_b = bin_a + r.nbins_s + 1;
    vec<double> d_v;
    d_v.reserve(samples_vvec[grid_idx].size());
    for (const auto& s : samples_vvec[grid_idx]) {
      if ((s.qidx == qidx) && (s.start < bin_b) && (bin_a < s.end)) continue;
      d_v.push_back(s.d);
    }

    if (d_v.size() < GammaModel::min_nsamples) {
      warn_msg("Too few null samples for p-value (after overlap filter or sparse grid)");
      r.percentile = std::numeric_limits<double>::quiet_NaN();
      r.fold = std::numeric_limits<double>::quiet_NaN();
      continue;
    }

    double cdf_val = std::numeric_limits<double>::quiet_NaN();
    double median = std::numeric_limits<double>::quiet_NaN();
    if (params.ecdf_test) {
      // ECDF-based test, no parameter estimation
      const size_t N = d_v.size();
      const auto lb = std::lower_bound(d_v.begin(), d_v.end(), r.d);
      const size_t rank = static_cast<size_t>(lb - d_v.begin());
      cdf_val = static_cast<double>(rank) / static_cast<double>(N);
      median = (N % 2 == 0) ? 0.5 * (d_v[N / 2 - 1] + d_v[N / 2]) : d_v[N / 2];
    } else {
      // Gamma scoring with parameter estimation
      const bool has_overlap = (d_v.size() < samples_vvec[grid_idx].size());
      GammaModel::params_t gp;

      if (!has_overlap) {
        // No overlap, hence just use the precomputed global Gamma
        gp = gp_v[grid_idx];
      } else {
        // Refit Gamma from filtered (non-overlapping) samples
        gp = GammaModel::fit(d_v);
      }

      gamma_distribution<double> g(gp.alpha, gp.beta);
      cdf_val = boost::math::cdf(g, r.d);
      median = quantile(g, 0.5);
    }
    if (r.rstrand) {
      // Two-sided for the strand closer to the reference
      r.percentile = 2.0 * std::min(cdf_val, 1.0 - cdf_val);
    } else {
      // One-sided for the opposite strand
      r.percentile = cdf_val;
    }
    r.fold = (median > eps) ? r.d / median : std::numeric_limits<double>::quiet_NaN();
  }
}

template<typename T>
void QIE<T>::sample_distances(uint64_t L, vec<snull_t>& samples_v) const
{ // TODO: Debug.
  const uint32_t W = params.hdist_th + 1;
  const uint64_t G = stride_len;
  const uint64_t N = L / G;
  vec<uint64_t> v(hdist_bound + 1);

  // Calculate weights as the number of valid starting positions for each query stride.
  vec<uint64_t> weights(qstrides_v.size());
  uint64_t total_weight = 0;
  for (size_t qi = 0; qi < qstrides_v.size(); ++qi) {
    const uint64_t nstrides = (qstrides_v[qi].nbins / G) + 1;
    weights[qi] = (nstrides > N) ? (nstrides - N) : 0; // Assign 0 weight if there are not enough strides
    total_weight += weights[qi];
  }
  if (total_weight == 0) return;

  // Sample windows with probability proportional to available positions.
  std::discrete_distribution<size_t> rvidx(weights.begin(), weights.end());

  for (uint64_t si = 0; si < params.nsamples; ++si) {
    const size_t qidx = rvidx(gen);
    const auto& qs = qstrides_v[qidx];
    const uint64_t nstrides = (qs.nbins / G) + 1;

    const uint64_t max_idx = nstrides - 1 - N;
    std::uniform_int_distribution<uint64_t> rvpos(0, max_idx);
    const uint64_t pidx = rvpos(gen);

    // Compute histogram as row[pidx + N] - row[pidx]
    std::fill(v.begin(), v.end(), 0);
    uint64_t t = 0;
    for (uint32_t d = 0; d < W; ++d) {
      v[d] = qs.hdisthist_v[((pidx + N) * W) + d] - qs.hdisthist_v[(pidx * W) + d];
      t += v[d];
    }

    // Compute miss count as the total nmers in the window minus hits
    const uint64_t bin_a = pidx * G;
    const uint64_t bin_b = bin_a + L;
    const uint64_t a = (bin_a << params.bin_shift);
    const uint64_t b = std::min((bin_b << params.bin_shift), qs.nmers);
    const uint64_t u = (b - a) - t;

    // Compute distance via LLH + Brent; widen upper bound when no hits
    double d = mle(llhf, v, u, t).first;
    if (std::isnan(d)) d = 0.75;
    samples_v.push_back({qidx, bin_a, bin_b, d});
  }
}

template<typename T>
void DIM<T>::aggregate_mer(uint32_t hdist_min, uint64_t i)
{
  // i is the bin index; multiple k-mers in the same bin accumulate here.
  if (hdist_min <= params.hdist_th) {
    t_q++;
    if (!params.enum_only) hdisthist_v[((i + 1) * (params.hdist_th + 1)) + hdist_min]++;
    add_to(sdc_v[i], llhf->get_sdc(hdist_min));
    add_to(fdc_v[i], llhf->get_fdc(hdist_min));
  } else {
    u_q++;
    add_to(sdc_v[i], llhf->get_sdc());
    add_to(fdc_v[i], llhf->get_fdc());
  }
}

template<typename T>
void DIM<T>::release_accumulators() noexcept
{
  fdc_v.clear();
  fdc_v.shrink_to_fit();
  sdc_v.clear();
  sdc_v.shrink_to_fit();
}

template<typename T>
interval_t DIM<T>::get_interval(uint64_t i, size_t idx) const
{
  if ((idx < eintervals_v.size()) && (i < eintervals_v[idx].size())) {
    return eintervals_v[idx][i];
  } else {
    return {nbins, nbins};
  }
}

template<typename T>
void QIE<T>::search_mers(const char* cseq, uint64_t len, DIM<T>& dim_fw, DIM<T>& dim_rc)
{
  uint64_t i = 0, j = 0, l = 0;
  uint32_t orrix, rcrix;
  uint64_t orenc64_bp, orenc64_lr, rcenc64_bp;
  for (; i < len; ++i) {
    if (__builtin_expect(SEQ_NT4_TABLE[cseq[i]] >= 4, 0)) {
      l = 0; // TODO: What to do for missing ones?
      continue;
    }
    ++l;
    if (l < k) {
      continue; // TODO: How to propagate the missing ones?
    }
    j = i - k + 1;
    if (l == k) {
      compute_encoding(cseq + j, cseq + i + 1, orenc64_lr, orenc64_bp);
    } else {
      update_encoding(cseq + i, orenc64_lr, orenc64_bp);
    }
    orenc64_bp &= mask_bp;
    orenc64_lr &= mask_lr;
    rcenc64_bp = revcomp_bp64(orenc64_bp, k);
    onmers++;
    const uint64_t bin_j = j >> params.bin_shift; // The bin index for this k-mer position
#ifdef CANONICAL
    if (rcenc64_bp < orenc64_bp) {
      orrix = lshf->compute_hash(orenc64_bp);
      const uint32_t off_fw = sketch->partial_offset(orrix);
      sketch->prefetch_offset_inc(off_fw);
      const enc_t enc_lr_fw = lshf->drop_ppos_lr(orenc64_lr);
      sketch->prefetch_offset_enc(off_fw);
      uint32_t hdist_fw;
      if (sketch->scan_bucket(off_fw, enc_lr_fw, hdist_fw)) {
        dim_fw.aggregate_mer(hdist_fw, bin_j);
      }
    } else {
      rcrix = lshf->compute_hash(rcenc64_bp);
      const uint32_t off_rc = sketch->partial_offset(rcrix);
      sketch->prefetch_offset_inc(off_rc);                                  // Phase 1
      const enc_t enc_lr_rc = lshf->drop_ppos_lr(bp64_to_lr64(rcenc64_bp)); // Phase 2
      sketch->prefetch_offset_enc(off_rc);                                  // Phase 3
      uint32_t hdist_rc;
      if (sketch->scan_bucket(off_rc, enc_lr_rc, hdist_rc)) {
        dim_rc.aggregate_mer(hdist_rc, bin_j);
      }
    }
#else
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
      dim_fw.aggregate_mer(hdist_fw, bin_j);
    }
    uint32_t hdist_rc;
    if (sketch->scan_bucket(off_rc, enc_lr_rc, hdist_rc)) {
      dim_rc.aggregate_mer(hdist_rc, bin_j);
    }
#endif /* CANONICAL */
  }
}

template<typename T>
void DIM<T>::inclusive_scan()
{
  assert(nbins > 0);
  const uint64_t s = nbins + 1;

  fdps_v.resize(s);
  sdps_v.resize(s);

  if constexpr (std::is_same_v<T, double>) {
    fdps_v[0] = 0.0;
    sdps_v[0] = 0.0;
    for (uint64_t i = 1; i < s; ++i) {
      fdps_v[i] = fdps_v[i - 1] + fdc_v[i - 1];
      sdps_v[i] = sdps_v[i - 1] + sdc_v[i - 1];
    }
  } else {
    fdps_v[0].fill(0.0);
    sdps_v[0].fill(0.0);
    simde__m512d fdps_acc = simde_mm512_setzero_pd();
    simde__m512d sdps_acc = simde_mm512_setzero_pd();
    for (uint64_t i = 1; i < s; ++i) {
      const simde__m512d fdc = simde_mm512_loadu_pd(fdc_v[i - 1].data());
      const simde__m512d sdc = simde_mm512_loadu_pd(sdc_v[i - 1].data());
      fdps_acc = simde_mm512_add_pd(fdps_acc, fdc);
      sdps_acc = simde_mm512_add_pd(sdps_acc, sdc);
      simde_mm512_storeu_pd(fdps_v[i].data(), fdps_acc);
      simde_mm512_storeu_pd(sdps_v[i].data(), sdps_acc);
    }
  }
}

template<typename T>
void DIM<T>::extrema_scan()
{
  const uint64_t s = nbins + 1;
  fdpmax_v.resize(s + 1);
  fdsmin_v.resize(s + 1);

  if constexpr (std::is_same_v<T, double>) {
    fdpmax_v[0] = -std::numeric_limits<double>::max();
    fdsmin_v[0] = std::numeric_limits<double>::max();
    std::inclusive_scan(
      fdps_v.begin() + 1, fdps_v.end(), fdpmax_v.begin() + 1, [](double a, double b) { return std::max(a, b); });
    std::inclusive_scan(
      fdps_v.rbegin(), fdps_v.rend() - 1, fdsmin_v.rbegin() + 1, [](double a, double b) { return std::min(a, b); });
    fdpmax_v[s] = std::numeric_limits<double>::max();
    fdsmin_v[s] = -std::numeric_limits<double>::max();
  } else {
    fdpmax_v[0].fill(-std::numeric_limits<double>::max());
    fdsmin_v[0].fill(std::numeric_limits<double>::max());
    simde__m512d fdpmax_acc = simde_mm512_loadu_pd(fdpmax_v[0].data());
    simde__m512d fdsmin_acc = simde_mm512_loadu_pd(fdsmin_v[0].data());
    for (uint64_t i = 1; i < s; ++i) {
      const simde__m512d fdps_front = simde_mm512_loadu_pd(fdps_v[i].data());
      const simde__m512d fdps_back = simde_mm512_loadu_pd(fdps_v[s - i].data());
      fdpmax_acc = simde_mm512_max_pd(fdpmax_acc, fdps_front);
      fdsmin_acc = simde_mm512_min_pd(fdsmin_acc, fdps_back);
      simde_mm512_storeu_pd(fdpmax_v[i].data(), fdpmax_acc);
      simde_mm512_storeu_pd(fdsmin_v[s - i].data(), fdsmin_acc);
    }
    fdpmax_v[s].fill(std::numeric_limits<double>::max());
    fdsmin_v[s].fill(-std::numeric_limits<double>::max());
  }
}

template<typename T>
void DIM<T>::extract_intervals_mx(const uint64_t tau, const size_t idx)
{
  uint64_t b_curr = 1;
  uint64_t b_prev = std::numeric_limits<uint64_t>::max(); // no interval yet

  if (at(fdps_v[nbins], idx) < at(fdps_v[1], idx)) {
    rintervals_v[idx].emplace_back(1, nbins);
    return;
  } // TODO: Is this correct???

  for (uint64_t a = 1; a <= nbins; ++a) {
    const double fdpmax_a = at(fdpmax_v[a - 1], idx);
    const double fdps_a = at(fdps_v[a], idx);

    if (fdpmax_a >= fdps_a) {
      continue; // Skip if a is not a prefix maxima of the prefix sum
    }

    while ((b_curr + 1) <= nbins && (at(fdsmin_v[b_curr + 1], idx) < fdps_a)) {
      ++b_curr; // Right maximal
    } // We increment b_curr to the last b with fdsmin_v[b] < fdps_a
    // There is no other a*>a where b*<b, so no valid (a, b) is missed

    const uint64_t b_star = b_curr;
    if (b_star < (a + tau)) {
      continue; // Skip if no valid right endpoint in [a+tau, nbins]
    }
    if (at(fdps_v[b_star], idx) >= fdps_a) {
      continue; // Skip if negative-sum check fails
    }
    // if (__builtin_expect(b_star == b_prev, 0)) {
    if (b_star == b_prev) {
      continue; // An earlier (leftmost) a already claimed this b*
    }

    if (at(fdps_v[b_star], idx) >= fdpmax_a) { // Left maximal
      rintervals_v[idx].emplace_back(a, b_star);
      b_prev = b_star;
    }
  }
}

template<typename T>
void DIM<T>::extract_intervals_sx(const uint64_t tau, const size_t idx)
{
  // O(n + k) where k = number of suffix minimum positions of fdps_v (k << n).
  // Every valid right endpoint b* is a suffix minimum of fdps_v.
  // Suffix minimum values are strictly increasing left-to-right,
  // Hence, the pointer into the list is monotone across record highs which is O(k) total.
  struct xy_t
  {
    uint64_t pos;
    double val;
  };
  vec<xy_t> xy_v;
  {
    double y_min = std::numeric_limits<double>::max();
    for (uint64_t j = nbins; j >= 1; --j) {
      const double v = at(fdps_v[j], idx);
      if (v < y_min) {
        y_min = v;
        xy_v.push_back({j, v});
      }
    }
    std::reverse(xy_v.begin(), xy_v.end());
  }
  if (xy_v.empty()) return;

  size_t yix_min = 0;
  double running_max = -std::numeric_limits<double>::max();
  uint64_t b_prev = std::numeric_limits<uint64_t>::max();

  for (uint64_t a = 1; a <= nbins; ++a) {
    const double fdps_a = at(fdps_v[a], idx);
    const double fdpmax_a = running_max;
    if (fdps_a > running_max) running_max = fdps_a;

    if (fdpmax_a >= fdps_a) continue; // Skip if a is not a record high

    // Advance yix_min to the last suffix minimum with val < fdps_a.
    while (yix_min + 1 < xy_v.size() && xy_v[yix_min + 1].val < fdps_a) {
      ++yix_min;
    }
    if (xy_v[yix_min].val >= fdps_a) continue; // Skip if no valid right endpoint with val < fdps_a

    const uint64_t b_star = xy_v[yix_min].pos;
    const double fdps_bstar = xy_v[yix_min].val;

    if (b_star < a + tau) continue; // Skip if no valid right endpoint in [a+tau, nbins]
    if (b_star == b_prev) continue; // Skip if b* was already claimed

    if (fdps_bstar >= fdpmax_a) { // Left maximal
      rintervals_v[idx].emplace_back(a, b_star);
      b_prev = b_star;
    }
  }
}

template<typename T>
void DIM<T>::expand_intervals(const double chisq_th, const size_t idx)
{
  if (rintervals_v[idx].empty()) {
    return;
  }

  double fdiff, sdiff, chisq_val;
  uint64_t a, ap, b, bp;
  ap = rintervals_v[idx][0].first;
  bp = rintervals_v[idx][0].second;

  for (uint64_t i = 1; i < rintervals_v[idx].size(); ++i) {
    a = rintervals_v[idx][i].first;
    b = rintervals_v[idx][i].second;
    fdiff = at(fdps_v[b], idx) - at(fdps_v[ap], idx);
    sdiff = at(sdps_v[ap], idx) - at(sdps_v[b], idx);
    // chisq_val = (sdiff > 0.0) ? (fdiff * fdiff) / sdiff : std::numeric_limits<double>::infinity();
    chisq_val = ((fdiff * fdiff) + eps) / (sdiff + eps);

    if ((chisq_val < chisq_th) && (a < bp)) {
      a = ap;
      // b = std::max(bp, b); // This is not necessary due to maximality and monotonicity of a's
    } else {
      eintervals_v[idx].emplace_back(ap, bp);
    }

    ap = a;
    bp = b;
  }

  fdiff = at(fdps_v[bp], idx) - at(fdps_v[ap], idx);
  sdiff = at(sdps_v[ap], idx) - at(sdps_v[bp], idx);
  // chisq_val = (sdiff > 0.0) ? (fdiff * fdiff) / sdiff : std::numeric_limits<double>::infinity();
  chisq_val = ((fdiff * fdiff) + eps) / (sdiff + eps);
  eintervals_v[idx].emplace_back(ap, bp);
  rintervals_v[idx].clear();
  rintervals_v[idx].shrink_to_fit();
}

template<typename T>
void DIM<T>::compute_prefhistsum()
{
  if (params.enum_only) return;
  const uint32_t W = params.hdist_th + 1;
  // Add each row to the previous in-place to get prefix sums
  for (uint64_t i = 0; i < nbins; ++i) {
    for (uint32_t d = 0; d < W; ++d) {
      hdisthist_v[((i + 1) * W) + d] += hdisthist_v[(i * W) + d];
    }
  }
}

template<typename T>
void DIM<T>::map_contiguous_segments(uint8_t th_bv, char sign)
{
  if (th_bv == 0) return;

  // Collect breakpoints via merge of the sorted (a, b) sequences per threshold.
  // Within each threshold, a's are strictly increasing and b's are strictly increasing.
  // Boundary breakpoints 1 and nbins ensure the entire sequence [1, nbins] is segmented.
  vec<uint64_t> pts = {1, nbins};
  for (size_t ti = 0; ti < WIDTH; ++ti) {
    if (!(th_bv & (1u << ti))) continue;
    const auto& e_v = eintervals_v[ti];
    if (e_v.empty()) continue;
    const size_t prev_size = pts.size();
    // Two-pointer merge of the sorted a-stream and sorted b-stream.
    size_t ai = 0, bi = 0;
    const size_t ne = e_v.size();
    while (ai < ne && bi < ne) {
      const uint64_t bin_a = e_v[ai].first;
      const uint64_t bin_b = e_v[bi].second;
      if (bin_a <= bin_b) {
        pts.push_back(bin_a);
        ++ai;
      } else {
        pts.push_back(bin_b);
        ++bi;
      }
    }
    while (ai < ne)
      pts.push_back(e_v[ai++].first);
    while (bi < ne)
      pts.push_back(e_v[bi++].second);
    std::inplace_merge(pts.begin(), pts.begin() + prev_size, pts.end());
  }
  pts.erase(std::unique(pts.begin(), pts.end()), pts.end());
  if (pts.size() < 2) return;

  // Each consecutive pair of breakpoints defines a contiguous segment [a, b].
  // The set of thresholds satisfied is constant within each segment.
  arr<size_t, WIDTH> ti_ix = {};
  for (size_t pi = 0; pi + 1 < pts.size(); ++pi) {
    const uint64_t a = pts[pi];
    const uint64_t b = pts[pi + 1];
    uint8_t mask = 0;
    for (size_t ti = 0; ti < WIDTH; ++ti) {
      if (!(th_bv & (1u << ti))) continue;
      // Advance past any interval whose right endpoint is before a.
      while (ti_ix[ti] < eintervals_v[ti].size() && eintervals_v[ti][ti_ix[ti]].second < a) {
        ++ti_ix[ti];
      }
      // By breakpoint construction [a,b] is either fully inside or fully outside the current interval.
      // Thus, first <= a is a sufficient check.
      if (ti_ix[ti] < eintervals_v[ti].size() && eintervals_v[ti][ti_ix[ti]].first <= a) {
        mask |= static_cast<uint8_t>(1u << ti);
      }
    }
    const double d = estimate_interval_distance(a - 1, b);
    segments_v.push_back({a, b - 1, d, mask, sign});
  }
  segments_v.back().end += 1;
}

template<typename T>
void DIM<T>::extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const
{
  const uint32_t W = params.hdist_th + 1;
  v.assign(hdist_bound + 1, 0);
  const simde__mmask8 mask = static_cast<simde__mmask8>((1u << W) - 1);
  const simde__m512i vb = simde_mm512_maskz_loadu_epi64(mask, &hdisthist_v[b * W]);
  const simde__m512i va = simde_mm512_maskz_loadu_epi64(mask, &hdisthist_v[a * W]);
  const simde__m512i vd = simde_mm512_sub_epi64(vb, va);
  simde_mm512_storeu_si512(v.data(), vd);
  const simde__m256i lo = simde_mm512_castsi512_si256(vd);
  const simde__m256i hi = simde_mm512_extracti64x4_epi64(vd, 1);
  const simde__m256i s4 = simde_mm256_add_epi64(lo, hi);
  const simde__m128i s4_lo = simde_mm256_castsi256_si128(s4);
  const simde__m128i s4_hi = simde_mm256_extracti128_si256(s4, 1);
  const simde__m128i s2 = simde_mm_add_epi64(s4_lo, s4_hi);
  t = simde_mm_extract_epi64(s2, 0) + simde_mm_extract_epi64(s2, 1);
  const uint64_t mers_b = std::min(b << params.bin_shift, nmers);
  const uint64_t mers_a = std::min(a << params.bin_shift, nmers);
  u = (mers_b - mers_a) - t;
}

template<typename T>
void QIE<T>::collect_segments(const DIM<T>& dim, double d_q, bool rc)
{
  const char strand = rc ? '-' : '+';
  const uint64_t L = enmers + k - 1;
  const auto& segments_v = dim.get_segments();

  const size_t prev_size = records_v.size();
  for (const auto& ab : segments_v) {
    // if (std::isnan(ab.d)) continue; // This might help for filtering
    const uint64_t a = ((ab.start - 1) << params.bin_shift) + 1;
    const uint64_t b = std::min((ab.end << params.bin_shift), enmers);
    assert(a < b);
    records_v.emplace_back(bix, L, a, b, strand, d_q, ab, dim.get_rstrand());
  }
  if (records_v.size() > prev_size) {
    records_v.back().b += k - 1;
  }
}

template<typename T>
void QIE<T>::save_qstride(const DIM<T>& dim)
{
  const uint32_t W = params.hdist_th + 1;
  const uint64_t G = stride_len;
  // Sample rows at 0, G, 2G, ..., floor(nbins/G)*G.
  const uint64_t nbins_q = dim.get_nbins();
  const uint64_t nstrides = (nbins_q / G) + 1;

  qstride_t qs;
  qs.nbins = nbins_q;
  qs.nmers = dim.get_nmers();
  qs.hdisthist_v.resize(nstrides * W);

  const auto& src = dim.get_hdisthist();
  for (uint64_t i = 0; i < nstrides; ++i) {
    const uint64_t rix = i * G;
    std::copy_n(src.begin() + (rix * W), W, qs.hdisthist_v.begin() + i * W);
  }

  qstrides_v.push_back(std::move(qs));
}

template<typename T>
double QIE<T>::compute_mle_dist(const vec<uint64_t>& v, uint64_t u, uint64_t t)
{
  return mle(llhf, v, u, t).first;
}

template<typename T>
double DIM<T>::estimate_interval_distance(uint64_t a, uint64_t b)
{
  vec<uint64_t> v;
  uint64_t u, t;
  extract_histogram(a, b, v, u, t);
  return mle(llhf, v, u, t).first;
}

template class QIE<double>;
template class QIE<cm512_t>;

template class DIM<double>;
template class DIM<cm512_t>;

template class LLH<double>;
template class LLH<cm512_t>;
