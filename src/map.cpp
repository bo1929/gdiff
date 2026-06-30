#include "map.hpp"
#include <boost/math/tools/minima.hpp>

namespace {
  struct pv_t
  {
    uint64_t pos;
    double val;
  };

  template<typename T>
  inline double at(T v, const size_t ix)
  {
    if constexpr (std::is_same_v<T, double>) {
      return v;
    } else {
      return v[ix];
    }
  }

  inline void add_to_acc(vec<uint64_t>& v_acc, uint64_t& u_acc, const vec<uint64_t>& v, uint64_t u)
  {
    simde__m512i s = simde_mm512_loadu_si512(v_acc.data());
    s = simde_mm512_add_epi64(s, simde_mm512_loadu_si512(v.data()));
    simde_mm512_storeu_si512(v_acc.data(), s);
    u_acc += u;
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
  diststat = std::make_shared<DistanceStat<T>>(params, llhf);
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  mask_bp = u64m >> ((32 - k) * 2);
  enum_only = params.enum_only;
  skip_test = (params.sample_size == 0);
  keep_hist = (!enum_only || !skip_test);
  coordinates_only = enum_only && skip_test;
  if (keep_hist) {
    // For SIMD alignment, bound is set to 8
    v_acc.assign(hdist_bound + 1, 0);
    v_scratch.assign(hdist_bound + 1, 0);
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
      warn_pmsg(qid_batch[bix], "skipped: sequence shorter than k-mer length ", "(len=", len, ", k=", k, ")");
      continue;
    }

    enmers = len - k + 1;
    nbins = (enmers + params.bin_size - 1) >> params.bin_shift;
    if (nbins < 2) {
      warn_pmsg(
        qid_batch[bix], "skipped: fewer than two bins after binning ", "(len=", len, ", bin_size=", params.bin_size, ")");
      continue;
    }
    if (params.tau_bin > nbins) {
      warn_pmsg(qid_batch[bix], "minimum length is exceeded; using the full query as the effective minimum ");
    }

    const uint64_t tau_eff = std::min(params.tau_bin, nbins) - 1;
    const size_t srprev = records_v.size(); // record index before this query

    if (params.canonical) {
      DIM<T> dim(params, llhf, nbins, enmers);
      search_mers(cseq, len, dim);
      dim.inclusive_scan();

      if (!coordinates_only) dim.compute_prefhistsum();

      double d_q = nanx();
      vec<uint64_t> v_q;
      uint64_t u_q = 0, t_q = 0;
      dim.total_histogram(v_q, u_q, t_q);
      d_q = llhf->mle(v_q.data(), u_q);
      d_q = validate_distance(d_q);

      dim.set_query_distance(d_q);
      dim.extrema_scan();

      if (enum_only) {
        extract_simple_intervals(dim, false, tau_eff);
        if (coordinates_only) continue; // skip MLE or significance
      }

      if (!skip_test) diststat->sample_null_pool(dim, tau_eff, bix);
      if (!enum_only) extract_ordered_intervals(dim, false, tau_eff);

      for (size_t ri = srprev; ri < records_v.size(); ++ri) {
        record_t& r = records_v[ri];
        r.d_q = d_q;
        r.d_diff = nanx();
      }

      add_to_acc(v_acc, u_acc, v_q, u_q);
    } else {
      DIM<T> dim_fw(params, llhf, nbins, enmers);
      DIM<T> dim_rc(params, llhf, nbins, enmers);
      search_mers(cseq, len, dim_fw, dim_rc);

      // Build prefix sums/histograms; d_q chooses per-threshold signs before extrema_scan.
      for (auto* dim : {&dim_fw, &dim_rc}) {
        dim->inclusive_scan();
        if (!coordinates_only) dim->compute_prefhistsum();
      }

      // Strand-wide MLE distances.
      double d_q_fw = nanx(), d_q_rc = nanx(), d_diff = nanx();
      vec<uint64_t> v_q_fw, v_q_rc;
      uint64_t u_q_fw = 0, u_q_rc = 0, t_q = 0;
      dim_fw.total_histogram(v_q_fw, u_q_fw, t_q);
      dim_rc.total_histogram(v_q_rc, u_q_rc, t_q);
      d_q_fw = llhf->mle(v_q_fw.data(), u_q_fw);
      d_q_rc = llhf->mle(v_q_rc.data(), u_q_rc);
      d_diff = strand_diff(d_q_fw, d_q_rc);
      d_q_fw = validate_distance(d_q_fw);
      d_q_rc = validate_distance(d_q_rc);

      // Determine direction per threshold, then signed extrema for extraction.
      dim_fw.set_query_distance(d_q_fw);
      dim_rc.set_query_distance(d_q_rc);

      for (auto* dim : {&dim_fw, &dim_rc}) {
        dim->extrema_scan();
      }

      if (enum_only) {
        extract_simple_intervals(dim_fw, false, tau_eff);
        extract_simple_intervals(dim_rc, true, tau_eff);
        if (coordinates_only) continue; // skip MLE or significance
      }

      // The lower-distance strand is the reference: rc when the difference > 0, else fw
      const bool is_rc = (!std::isnan(d_diff)) && d_diff > 0.0;
      if (!skip_test) diststat->sample_null_pool(is_rc ? dim_rc : dim_fw, tau_eff, bix);

      if (!enum_only) {
        extract_ordered_intervals(dim_fw, false, tau_eff);
        extract_ordered_intervals(dim_rc, true, tau_eff);
      }

      for (size_t ri = srprev; ri < records_v.size(); ++ri) {
        record_t& r = records_v[ri];
        r.d_q = r.is_rc ? d_q_rc : d_q_fw;
        r.d_diff = d_diff;
      }

      add_to_acc(v_acc, u_acc, is_rc ? v_q_rc : v_q_fw, is_rc ? u_q_rc : u_q_fw);
    }
  }

  if (keep_hist) d_acc = llhf->mle(v_acc.data(), u_acc);
  if (!skip_test) {
    for (auto& r : records_v)
      diststat->test_significance(r, params.sample_size, qid_batch[r.bix]);
    diststat->benjamini_hochberg_correction(records_v);
  }
  report_contiguous(sout, rid);
}

template<typename T>
DIM<T>::DIM(const params_t<T>& params, const llh_sptr_t<T>& llhf, uint64_t nbins, uint64_t nmers)
  : params(params)
  , llhf(llhf)
  , nbins(nbins)
  , nmers(nmers)
  , keep_hist(!params.enum_only || params.sample_size > 0)
{
  fdc_v.resize(nbins); // Alternative?: fdc_v.reserve(nbins);
  sdc_v.resize(nbins); // Alternative?: sdc_v.reserve(nbins);
  thneg_v.fill(false);
  thrank_v.resize(WIDTH);
  for (size_t i = 0; i < WIDTH; ++i)
    thrank_v[i] = i;
  if (keep_hist) {
    // Note that aggregate_mer() accumulates into rows 1 to nbins.
    // At the end, compute_prefhistsum() converts in-place.
    hdisthist_v.assign((nbins + 1) * (params.hdist_th + 1), 0);
    // First delta (HD threshold) + 1 values are zeros, same layout as fdps_v and sdps_v.
  } else {
    // Otherwise, just keep a global histogram for the query.
    hdisthist_v.assign(hdist_bound + 1, 0);
  }
}

template<typename T>
void DIM<T>::set_query_distance(const double d_q)
{ // {{{ OK
  const bool is_valid = std::isfinite(d_q);

  if constexpr (std::is_same_v<T, double>) {
    const double t = params.dist_th;
    thneg_v.front() = is_valid && t > d_q;
    thrank_v = {0};
  } else {
    arr<vi_t, WIDTH> tp;
    for (size_t i = 0; i < WIDTH; ++i) {
      const double t = at(params.dist_th, i);
      thneg_v[i] = is_valid && t > d_q;
      tp[i] = {t, i};
    }

    if (!is_valid) {
      std::sort(tp.begin(), tp.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    } else {
      std::sort(tp.begin(), tp.end(), [&](const auto& a, const auto& b) {
        const double sa = a.first <= d_q ? a.first : (2.0 * d_ub - a.first);
        const double sb = b.first <= d_q ? b.first : (2.0 * d_ub - b.first);
        return sa < sb;
      });
    }

    thrank_v.resize(WIDTH);
    for (size_t i = 0; i < WIDTH; ++i)
      thrank_v[i] = tp[i].second;
  }

  apply_threshold_signs();
} // }}}

template<typename T>
void DIM<T>::apply_threshold_signs()
{ // {{{ OK
  // Flip the sign of fdps_v lanes whose threshold is on the opposite side of d_q.
  if constexpr (std::is_same_v<T, double>) {
    if (!thneg_v.front()) return;
    for (uint64_t i = 1; i < fdps_v.size(); ++i) {
      fdps_v[i] = -fdps_v[i];
    }
  } else {
    alignas(64) double xor_mask[WIDTH];
    bool flip_sign = false;
    for (size_t i = 0; i < WIDTH; ++i) {
      xor_mask[i] = thneg_v[i] ? -0.0 : 0.0;
      if (thneg_v[i]) {
        flip_sign = true;
      }
    }
    if (!flip_sign) return;
    const simde__m512d s = simde_mm512_loadu_pd(xor_mask);
    for (uint64_t i = 1; i < fdps_v.size(); ++i) {
      simde__m512d v = simde_mm512_loadu_pd(fdps_v[i].data());
      v = simde_mm512_xor_pd(v, s);
      simde_mm512_storeu_pd(fdps_v[i].data(), v);
    }
  }
} // }}}

template<typename T>
void DIM<T>::aggregate_mer(uint32_t hdist_min, uint64_t i)
{
  // The bin index is i, multiple k-mers in the same bin accumulate here
  if (hdist_min <= params.hdist_th) {
    t_q++;
    if (keep_hist) {
      hdisthist_v[((i + 1) * (params.hdist_th + 1)) + hdist_min]++;
    } else {
      hdisthist_v[hdist_min]++;
    }
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
{ // {{{ OK: (probably not quite needed)
  fdc_v.clear();
  fdc_v.shrink_to_fit();
  sdc_v.clear();
  sdc_v.shrink_to_fit();
} // }}}

template<typename T>
void QIE<T>::search_mers(const char* cseq, uint64_t len, DIM<T>& dim)
{
  uint64_t i = 0, j = 0, l = 0;
  uint32_t orrix, rcrix;
  uint64_t orenc64_bp, orenc64_lr, rcenc64_bp;
  for (; i < len; ++i) {
    if (__builtin_expect(SEQ_NT4_TABLE[cseq[i]] >= 4, 0)) {
      // TODO: What to do for missing ones? Masked repeats should be handled here, too
      l = 0;
      continue;
    }
    ++l;
    if (l < k) {
      // TODO: What to do for missing ones? Masked repeats should be handled here, too
      continue;
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
    if (rcenc64_bp < orenc64_bp) {
      orrix = lshf->compute_hash(orenc64_bp);
      const uint32_t off_fw = sketch->partial_offset(orrix);
      sketch->prefetch_offset_inc(off_fw);
      const enc_t enc_lr_fw = lshf->drop_ppos_lr(orenc64_lr);
      sketch->prefetch_offset_enc(off_fw);
      uint32_t hdist_fw;
      if (sketch->scan_bucket(off_fw, enc_lr_fw, hdist_fw)) {
        dim.aggregate_mer(hdist_fw, bin_j);
      }
    } else {
      rcrix = lshf->compute_hash(rcenc64_bp);
      const uint32_t off_rc = sketch->partial_offset(rcrix);
      sketch->prefetch_offset_inc(off_rc);                                  // Phase 1
      const enc_t enc_lr_rc = lshf->drop_ppos_lr(bp64_to_lr64(rcenc64_bp)); // Phase 2
      sketch->prefetch_offset_enc(off_rc);                                  // Phase 3
      uint32_t hdist_rc;
      if (sketch->scan_bucket(off_rc, enc_lr_rc, hdist_rc)) {
        dim.aggregate_mer(hdist_rc, bin_j);
      }
    }
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
      // TODO: What to do for missing ones? Masked repeats should be handled here, too
      l = 0;
      continue;
    }
    ++l;
    if (l < k) {
      // TODO: What to do for missing ones? Masked repeats should be handled here, too
      continue;
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
  }
}

template<typename T>
void DIM<T>::inclusive_scan()
{ // {{{ OK
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
} // }}}

template<typename T>
void DIM<T>::extrema_scan()
{ // {{{ OK
  const uint64_t s = nbins + 1;
  fdpmax_v.resize(s + 1);
  fdsmin_v.resize(s + 1);

  if constexpr (std::is_same_v<T, double>) {
    fdpmax_v[0] = ninf();
    fdsmin_v[0] = pinf();
    std::inclusive_scan(
      fdps_v.begin() + 1, fdps_v.end(), fdpmax_v.begin() + 1, [](double a, double b) { return std::max(a, b); });
    std::inclusive_scan(
      fdps_v.rbegin(), fdps_v.rend() - 1, fdsmin_v.rbegin() + 1, [](double a, double b) { return std::min(a, b); });
    fdpmax_v[s] = pinf();
    fdsmin_v[s] = ninf();
  } else {
    fdpmax_v[0].fill(ninf());
    fdsmin_v[0].fill(pinf());
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
    fdpmax_v[s].fill(pinf());
    fdsmin_v[s].fill(ninf());
  }
} // }}}

// Find maximal intervals [a, b] within [lix, rix] where the prefix sum drops below a prior maximum.
//
// Conditions for a valid interval (a, b):
//   1. a is a strict prefix maximum:  fdps[a] > fdpmax[a-1]
//   2. b is right-maximal:            fdsmin[b] < fdps[a] but fdsmin[b+1] >= fdps[a]
//   3. Minimum length:                b >= a + tau
//   4. Negative sum:                  fdps[b] < fdps[a]
//   5. Left-maximal (non-redundant):  fdps[b] >= fdpmax[a-1]
//   6. b not claimed by earlier a:    b != b_prev
//
// Early return: if prefix sum ends below where it started (fdps[rix] < fdps[lix]), the entire [lix, rix] is one interval.
template<typename T>
void DIM<T>::extract_intervals_mx(const uint64_t tau, const uint64_t lix, const uint64_t rix, const size_t ix)
{
  uint64_t b_curr = lix;
  uint64_t b_prev = std::numeric_limits<uint64_t>::max();

  if (rix >= lix + tau && at(fdps_v[rix], ix) < at(fdps_v[lix], ix)) {
    intervals_v[ix].emplace_back(lix, rix);
    return;
  }

  for (uint64_t a = lix; a <= rix; ++a) {
    const double fdpmax_a = at(fdpmax_v[a - 1], ix);
    const double fdps_a = at(fdps_v[a], ix);

    if (fdpmax_a >= fdps_a) continue; // Condition 1

    while ((b_curr + 1) <= rix && (at(fdsmin_v[b_curr + 1], ix) < fdps_a))
      ++b_curr; // Condition 2

    const uint64_t b_star = b_curr;
    if (b_star < (a + tau)) continue;               // Condition 3
    if (at(fdps_v[b_star], ix) >= fdps_a) continue; // Condition 4
    if (b_star == b_prev) continue;                 // Condition 6

    if (at(fdps_v[b_star], ix) >= fdpmax_a) { // Condition 5
      intervals_v[ix].emplace_back(a, b_star);
      b_prev = b_star;
    }
  }
}

template<typename T>
void DIM<T>::extract_intervals_sx(const uint64_t tau, const uint64_t lix, const uint64_t rix, const size_t ix)
{ // {{{ ~OK: more or less identical to mx
  // Every valid right endpoint b* is a suffix minimum of fdps_v
  // Suffix minimum values are strictly increasing left-to-right,
  // Hence, the pointer into the list is monotone across record highs which is O(k) total
  const uint64_t gap_len = rix - lix + 1;
  if (gap_len >= 1 + tau && at(fdps_v[rix], ix) < at(fdps_v[lix], ix)) {
    intervals_v[ix].emplace_back(lix, rix);
    return;
  }

  vec<pv_t> pv_v;
  {
    double y_min = pinf();
    for (uint64_t j = rix; j >= lix; --j) {
      const double v = at(fdps_v[j], ix);
      if (v < y_min) {
        y_min = v;
        pv_v.push_back({j, v});
      }
    }
    std::reverse(pv_v.begin(), pv_v.end());
  }
  if (pv_v.empty()) return;

  size_t yix_min = 0;
  double running_max = ninf();
  uint64_t b_prev = std::numeric_limits<uint64_t>::max();

  for (uint64_t a = lix; a <= rix; ++a) {
    const double fdps_a = at(fdps_v[a], ix);
    const double fdpmax_a = running_max;
    if (fdps_a > running_max) running_max = fdps_a;

    if (fdpmax_a >= fdps_a) continue; // Skip if a is not a record high

    // Advance yix_min to the last suffix minimum with val < fdps_a
    while (yix_min + 1 < pv_v.size() && pv_v[yix_min + 1].val < fdps_a) {
      ++yix_min;
    }
    // Skip if no valid right endpoint with val < fdps_a
    if (pv_v[yix_min].val >= fdps_a) continue;

    const uint64_t b_star = pv_v[yix_min].pos;
    const double fdps_bstar = pv_v[yix_min].val;

    if (b_star < a + tau) continue; // Skip if no valid right endpoint in [a+tau, nbins]
    if (b_star == b_prev) continue; // Skip if b* was already claimed

    if (fdps_bstar >= fdpmax_a) {              // Left maximal
      intervals_v[ix].emplace_back(a, b_star); // 1-based inclusive coordinates
      b_prev = b_star;
    }
  }
} // }}}

template<typename T>
void DIM<T>::expand_intervals(const double chisq_th, const size_t ix)
{ // {{{ OK
  auto& iv_ix = intervals_v[ix];
  if (iv_ix.empty()) return;

  double fdiff, sdiff, chisq_val;
  uint64_t a, ap, b, bp;
  size_t w = 0;
  ap = iv_ix[0].a;
  bp = iv_ix[0].b;

  for (size_t i = 1; i < iv_ix.size(); ++i) {
    a = iv_ix[i].a;
    b = iv_ix[i].b;
    fdiff = at(fdps_v[b], ix) - at(fdps_v[ap], ix);
    sdiff = at(sdps_v[ap], ix) - at(sdps_v[b], ix);
    // chisq_val = (sdiff > 0.0) ? (fdiff * fdiff) / sdiff : std::numeric_limits<double>::infinity();
    chisq_val = ((fdiff * fdiff) + eps) / (sdiff + eps);

    if ((chisq_val < chisq_th) && (a < bp)) {
      a = ap;
      // b = std::max(bp, b); // This is not necessary due to maximality and monotonicity of a's
    } else {
      iv_ix[w++] = {ap, bp}; // 1-based inclusive coordinates
    }

    ap = a;
    bp = b;
  }

  iv_ix[w++] = {ap, bp}; // 1-based inclusive coordinates
  iv_ix.resize(w);
} // }}}

template<typename T>
void DIM<T>::compute_prefhistsum()
{ // {{{ OK
  if (!keep_hist) return;
  const uint32_t W = params.hdist_th + 1;
  // Add each row to the previous in-place to get prefix sums
  for (uint64_t i = 0; i < nbins; ++i) {
    for (uint32_t d = 0; d < W; ++d) {
      hdisthist_v[((i + 1) * W) + d] += hdisthist_v[(i * W) + d];
    }
  }
} // }}}

template<typename T>
void DIM<T>::total_histogram(vec<uint64_t>& v, uint64_t& u, uint64_t& t) const
{
  if (keep_hist) {
    extract_histogram(0, nbins, v, u, t);
  } else {
    v = hdisthist_v;
    u = u_q; // explicit HD misses only; omits unscanned positions (e.g. Ns)
    t = t_q;
  }
}

template<typename T>
void DIM<T>::extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const
{ // {{{ ~OK (TODO: counting misses or Ns)
  const uint32_t W = params.hdist_th + 1;
  v.resize(hdist_bound + 1);
  // v.assign(hdist_bound + 1, 0);
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
  t = simde_mm_extract_epi64(s2, 0) + simde_mm_extract_epi64(s2, 1);
  const uint64_t mers_b = std::min(b << params.bin_shift, nmers);
  const uint64_t mers_a = std::min(a << params.bin_shift, nmers);
  // The miss count (u) is based on all positions in [a,b), but search_mers() skips Ns.
  u = (mers_b - mers_a) - t;
  // A better solution is needed for Ns in this case, skipping does not work.
} // }}}

template<typename T>
void QIE<T>::extract_simple_intervals(DIM<T>& dim, bool is_rc, uint64_t tau_eff)
{ // {{{ ~OK
  if constexpr (std::is_same_v<T, double>) {
    dim.extract_intervals_mx(tau_eff, 1, nbins);
    dim.expand_intervals(params.chisq);
    for (const auto& iv : dim.get_intervals(0))
      emit_record(dim, iv.a, iv.b + 1, 0, is_rc);
  } else {
    for (size_t ix = 0; ix < WIDTH; ++ix) {
      dim.extract_intervals_mx(tau_eff, 1, nbins, ix);
      dim.expand_intervals(params.chisq, ix);
      for (const auto& iv : dim.get_intervals(ix))
        emit_record(dim, iv.a, iv.b + 1, ix, is_rc);
    }
  }
} // }}}

template<typename T>
void QIE<T>::extract_ordered_intervals(DIM<T>& dim, bool is_rc, uint64_t tau_eff)
{ // {{{ ~OK
  const uint64_t nbins = dim.get_nbins();
  bp_v.clear();

  const auto& thrank = dim.get_thrank();
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

    const auto& iv_v = dim.get_intervals(ix);
    for (const auto& iv : iv_v) {
      bp_v.push_back({iv.a, iv.b + 1, ix});
    }

    sbprev = bp_v.size();
  }
  merge_from(sbprev);

  if (bp_v.empty()) {
    // No intervals extracted; report the full query.
    emit_record(dim, 1, nbins + 1, size_t(-1), is_rc);
  } else {
    // uint64_t prev = 1;
    for (const auto& s : bp_v) {
      // TODO: Is full segmentation desired? Even when they are too short?
      // if (s.a_bin > prev) emit_record(dim, prev, s.a_bin, size_t(-1), is_rc);
      emit_record(dim, s.a_bin, s.b_bin, s.ix, is_rc);
      // prev = s.b_bin;
    }
    //if (prev <= nbins) emit_record(dim, prev, nbins + 1, size_t(-1), is_rc);
  }
} // }}}

template<typename T>
void QIE<T>::emit_record(DIM<T>& dim, uint64_t a_bin, uint64_t b_bin, size_t th_ix, bool is_rc)
{ // {{{ OK
  const uint64_t L = enmers + k - 1;
  const bool is_last = (b_bin == dim.get_nbins() + 1);

  double d = nanx();
  double I = nanx();
  if (!coordinates_only) {
    uint64_t u, t;
    dim.extract_histogram(a_bin - 1, b_bin - 1, v_scratch, u, t);
    d = llhf->mle(v_scratch.data(), u);
    I = llhf->compute_fisher_info(d);
    // MLE at the upper boundary is unreliable - mark NA.
    d = validate_distance(d);
  }

  const interval_t bin_iv{a_bin, b_bin};
  const interval_t seq_iv = get_coordinates(bin_iv, params.bin_shift, enmers, k, is_last);
  // assert(seq_iv.a < seq_iv.b);

  records_v.emplace_back(bix, L, seq_iv, bin_iv, is_rc, d, I, th_ix);
} // }}}

template<typename T>
xy_t QIE<T>::get_distance_bin(const record_t& r, const arr<double, WIDTH>& th_v) const
{ // {{{ OK
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
      d_range.second = (pos + 1 < WIDTH) ? th_v[pos + 1] : d_ub;
    }
  } else if (std::isfinite(r.d) || std::isfinite(r.d_q)) {
    // No threshold triggered: bracket the interval's own distance
    const double d_anchor = std::isfinite(r.d) ? r.d : r.d_q;
    const auto it = std::lower_bound(th_v.begin(), th_v.end(), d_anchor);
    if (it != th_v.begin()) d_range.first = *(it - 1);
    if (it != th_v.end()) d_range.second = *it;
  }
  return d_range;
} // }}}

template<typename T>
void QIE<T>::report_contiguous(std::ostream& sout, const str& rid) const
{ // {{{ ???
  // TODO: Revisit this and design a better format!
  arr<double, WIDTH> th_v{};
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
                qid_batch[r.bix],
                r.L,
                r.seq_iv.a,
                r.seq_iv.b,
                rid,
                r.d,
                static_cast<uint32_t>(mask),
                d_bin.str(),
                r.d_q,
                d_acc,
                r.percentile,
                r.fold,
                r.qvalue)
        << '\n';
    } else {
      write_tsv(sout,
                qid_batch[r.bix],
                r.L,
                r.seq_iv.a,
                r.seq_iv.b,
                report_strand(r.is_rc, r.d_diff),
                static_cast<uint32_t>(r.is_rc),
                rid,
                r.d,
                static_cast<uint32_t>(mask),
                d_bin.str(),
                r.d_q,
                r.d_diff,
                d_acc,
                r.percentile,
                r.fold,
                r.qvalue)
        << '\n';
    }
  }
} // }}}

template class QIE<double>;
template class QIE<cm512_t>;

template class DIM<double>;
template class DIM<cm512_t>;

template class LLH<double>;
template class LLH<cm512_t>;
