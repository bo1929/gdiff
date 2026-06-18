#include "map.hpp"
#include <boost/math/tools/minima.hpp>

namespace {
  template<typename T>
  inline double at(T v, const size_t ix)
  {
    if constexpr (std::is_same_v<T, double>) {
      return v;
    } else {
      return v[ix];
    }
  }

  // Generates two independent N(0,1) samples from two U(0,1) draws
  double sample_box_muller(std::mt19937& rng)
  {
    std::uniform_real_distribution<double> U(0.0, 1.0);
    const double u1 = U(rng);
    const double u2 = U(rng);
    const double r = std::sqrt(-2.0 * std::log(u1));
    const double phi = 2.0 * M_PI * u2;
    // return {r * std::cos(phi)
    return r * std::cos(phi); // r * std::sin(phi) is the other sample
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
  // Precompute sorted threshold order based on the absolute value, positive before negative if tied.
  if constexpr (std::is_same_v<T, double>) {
    sthix_v = {0};
  } else {
    for (size_t ti = 0; ti < WIDTH; ++ti)
      sthix_v.push_back(ti);
    std::sort(sthix_v.begin(), sthix_v.end(), [&](size_t a, size_t b) {
      const double va = std::abs(at(llhf->get_extrema(), a));
      const double vb = std::abs(at(llhf->get_extrema(), b));
      if (va != vb) return va < vb;
      return at(llhf->get_sign(), a) > at(llhf->get_sign(), b); // positive first
    });
  }
  for (size_t ix : sthix_v) {
    if (at(llhf->get_sign(), ix) > 0)
      ltix_v.push_back(ix);
    else
      gtix_v.push_back(ix);
  }
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

    DIM<T> dim_fw(params, llhf, nbins, enmers);
    DIM<T> dim_rc(params, llhf, nbins, enmers);
    search_mers(cseq, len, dim_fw, dim_rc);

    const uint64_t tau_eff = std::min(params.tau_bin, nbins) - 1;
    const size_t srprev = records_v.size(); // record index before this query

    // Extract intervals for all thresholds.
    for (auto* dim : {&dim_fw, &dim_rc}) {
      dim->inclusive_scan();
      dim->extrema_scan();
      if (!coordinates_only) {
        dim->compute_prefhistsum();
      }
    }

    if (enum_only) {
      extract_simple_intervals(dim_fw, false, tau_eff);
      extract_simple_intervals(dim_rc, true, tau_eff);
      if (coordinates_only) continue; // skip MLE or significance
    }

    // Compute the per-query MLE on each strand.
    vec<uint64_t> v_q_fw, v_q_rc;
    uint64_t u_q_fw = 0, u_q_rc = 0, t_q;
    dim_fw.extract_histogram(0, nbins, v_q_fw, u_q_fw, t_q);
    dim_rc.extract_histogram(0, nbins, v_q_rc, u_q_rc, t_q);
    const double d_q_fw = llhf->mle(v_q_fw.data(), u_q_fw);
    const double d_q_rc = llhf->mle(v_q_rc.data(), u_q_rc);
    const double d_diff = strand_diff(d_q_fw, d_q_rc);

    // The lower-distance strand is the reference: rc when the difference > 0, else fw.
    const bool is_rc = !std::isnan(d_q_rc) && (std::isnan(d_q_fw) || d_q_rc < d_q_fw);
    if (!skip_test) sample_null_pool(is_rc ? dim_rc : dim_fw, tau_eff);

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

  if (keep_hist) d_acc = llhf->mle(v_acc.data(), u_acc);
  if (!skip_test) test_significance(params.sample_size);
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
  if (keep_hist) {
    // Note that aggregate_mer() accumulates into rows 1 to nbins
    // At the end, compute_prefhistsum() converts in-place
    hdisthist_v.assign((nbins + 1) * (params.hdist_th + 1), 0);
    // First delta (HD threshold) + 1 values are zeros, same layout as fdps_v/sdps_v
  }
}

template<typename T>
void DIM<T>::aggregate_mer(uint32_t hdist_min, uint64_t i)
{
  // The bin index is i, multiple k-mers in the same bin accumulate here
  if (hdist_min <= params.hdist_th) {
    t_q++;
    if (keep_hist) hdisthist_v[((i + 1) * (params.hdist_th + 1)) + hdist_min]++;
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
interval_t DIM<T>::get_interval(uint64_t i, size_t ix) const
{
  if ((ix < intervals_v.size()) && (i < intervals_v[ix].size())) {
    return intervals_v[ix][i];
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

// Find maximal intervals [a, b] within [lix, rix] where the prefix sum drops below a prior maximum
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
{
  // Every valid right endpoint b* is a suffix minimum of fdps_v
  // Suffix minimum values are strictly increasing left-to-right,
  // Hence, the pointer into the list is monotone across record highs which is O(k) total
  const uint64_t gap_len = rix - lix + 1;
  if (gap_len >= 1 + tau && at(fdps_v[rix], ix) < at(fdps_v[lix], ix)) {
    intervals_v[ix].emplace_back(lix, rix);
    return;
  }

  struct xy_t
  {
    uint64_t pos;
    double val;
  };
  vec<xy_t> xy_v;
  {
    double y_min = std::numeric_limits<double>::max();
    for (uint64_t j = rix; j >= lix; --j) {
      const double v = at(fdps_v[j], ix);
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

  for (uint64_t a = lix; a <= rix; ++a) {
    const double fdps_a = at(fdps_v[a], ix);
    const double fdpmax_a = running_max;
    if (fdps_a > running_max) running_max = fdps_a;

    if (fdpmax_a >= fdps_a) continue; // Skip if a is not a record high

    // Advance yix_min to the last suffix minimum with val < fdps_a
    while (yix_min + 1 < xy_v.size() && xy_v[yix_min + 1].val < fdps_a) {
      ++yix_min;
    }
    // Skip if no valid right endpoint with val < fdps_a
    if (xy_v[yix_min].val >= fdps_a) continue;

    const uint64_t b_star = xy_v[yix_min].pos;
    const double fdps_bstar = xy_v[yix_min].val;

    if (b_star < a + tau) continue; // Skip if no valid right endpoint in [a+tau, nbins]
    if (b_star == b_prev) continue; // Skip if b* was already claimed

    if (fdps_bstar >= fdpmax_a) {              // Left maximal
      intervals_v[ix].emplace_back(a, b_star); // 1-based inclusive coordinates
      b_prev = b_star;
    }
  }
}

template<typename T>
void DIM<T>::expand_intervals(const double chisq_th, const size_t ix)
{
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
}

template<typename T>
void DIM<T>::compute_prefhistsum()
{
  if (!keep_hist) return;
  const uint32_t W = params.hdist_th + 1;
  // Add each row to the previous in-place to get prefix sums
  for (uint64_t i = 0; i < nbins; ++i) {
    for (uint32_t d = 0; d < W; ++d) {
      hdisthist_v[((i + 1) * W) + d] += hdisthist_v[(i * W) + d];
    }
  }
}

template<typename T>
void DIM<T>::extract_histogram(uint64_t a, uint64_t b, vec<uint64_t>& v, uint64_t& u, uint64_t& t) const
{
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
  // TODO: The miss count (u) is based on all positions in [a,b), but search_mers() skips Ns.
  u = (mers_b - mers_a) - t;
  // TODO: A better solution is needed for Ns in this case, skipping does not work.
}

template<typename T>
void QIE<T>::extract_simple_intervals(DIM<T>& dim, bool is_rc, uint64_t tau_eff)
{
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
}

template<typename T>
void QIE<T>::extract_ordered_intervals(DIM<T>& dim, bool is_rc, uint64_t tau_eff)
{
  const uint64_t nbins = dim.get_nbins();
  bp_v.clear();

  for (const auto* ix_v : {&ltix_v, &gtix_v}) {
    if (ix_v->empty()) continue;

    const size_t sbprev = bp_v.size();

    auto merge_bp = [&]() { // Merge newly appended entries in a sorted manner.
      if (bp_v.size() <= sbprev) return;
      auto nit = std::next(bp_v.begin(), sbprev);
      auto cmp = [](const bp_t& a, const bp_t& b) { return a.a_bin < b.a_bin; };
      std::sort(nit, bp_v.end(), cmp);
      std::inplace_merge(bp_v.begin(), nit, bp_v.end(), cmp);
    };

    for (size_t ix : *ix_v) {
      merge_bp();

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

      for (const auto& iv : dim.get_intervals(ix))
        bp_v.push_back({iv.a, iv.b + 1, ix});
    }

    merge_bp();
  }

  if (bp_v.empty()) {
    // No intervals extracted; report the full query.
    // emit_record(dim, 1, nbins + 1, size_t(-1), is_rc);
  } else {
    // uint64_t prev = 1;
    for (const auto& s : bp_v) {
      // if (s.a_bin > prev) emit_record(dim, prev, s.a_bin, size_t(-1), is_rc);
      emit_record(dim, s.a_bin, s.b_bin, s.ix, is_rc);
      // prev = s.b_bin;
    }
    // if (prev <= nbins) emit_record(dim, prev, nbins + 1, size_t(-1), is_rc);
  }
}

template<typename T>
void QIE<T>::emit_record(DIM<T>& dim, uint64_t a_bin, uint64_t b_bin, size_t th_ix, bool is_rc)
{
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
  assert(seq_iv.a < seq_iv.b);

  records_v.emplace_back(bix, L, seq_iv, bin_iv, is_rc, d, I, th_ix);
}

template<typename T>
void QIE<T>::sample_null_pool(DIM<T>& dim, uint64_t tau_eff)
{
  const uint64_t nbins_q = dim.get_nbins();
  const uint64_t win_len = tau_eff + 1; // fixed window length in bins
  if (nbins_q < win_len) return;

  const uint64_t max_nwindows = nbins_q / win_len;
  const uint64_t nsamples = std::min<uint64_t>(params.sample_size, max_nwindows);

  std::uniform_int_distribution<uint64_t> rvstart(0, nbins_q - win_len);
  vec<uint64_t> v(hdist_bound + 1);

  for (uint64_t s = 0; s < nsamples; ++s) {
    const uint64_t x = rvstart(gen);
    const uint64_t a_bin = x + 1;
    const uint64_t b_bin = x + win_len + 1;

    uint64_t u, t;
    dim.extract_histogram(a_bin - 1, b_bin - 1, v, u, t);
    double d = llhf->mle(v.data(), u);
    double I = llhf->compute_fisher_info(d);
    d = validate_distance(d);
    if (std::isfinite(d) && std::isfinite(I)) null_v.push_back({d, I, bix, {a_bin, b_bin}});
  }
}

template<typename T>
bool QIE<T>::sample_metropolis_hastings(const xy_t& init)
{
  // Metropolis-Hastings posterior draws kept (S) and burn-in iterations (B) per null window.
  bool is_valid = true;
  // Calculates the negative log-likelihood of the current distance
  auto f = [&](const double& D) { return (*llhf)(D); };
  double d = init.first;
  // The target posterior is proportional to exp(-nll(D)).
  const double step = init.second; // Proposal scale; sqrt(1/I) approximates the posterior sigma
  double nll = f(d);

  std::uniform_real_distribution<double> ruv(0, 1);

  samples_v.reserve(samples_v.size() + S);

  for (uint64_t iter = 0; iter < B + S; ++iter) {
    // Reflect proposal off boundaries to stay within [d_eps, d_ub].
    double d_i = d + step * sample_box_muller(gen);
    // Reflect instead of reject - better mixing near boundaries
    if (d_i <= 0.0) d_i = -d_i;
    if (d_i >= d_ub) d_i = 2.0 * d_ub - d_i;
    d_i = std::clamp(d_i, eps, d_ub - eps);

    const double nll_i = f(d_i);
    // log acceptance ratio = logL(prop) - logL(cur) = f(cur) - f(prop)
    const double log_alpha = nll - nll_i;
    if (std::log(ruv(gen)) < log_alpha) {
      d = d_i;
      nll = nll_i;
    }
    if (iter >= B) {
      const double I = llhf->compute_fisher_info(d);
      if ((std::isfinite(I) && std::isfinite(d)) && (d > 0.0 && I > 0.0)) {
        samples_v.push_back({d, I});
      } else {
        is_valid = false;
      }
    }
  }
  return is_valid;
}

template<typename T>
void QIE<T>::filter_sample(const record_t& r, vec<xy_t>& rsample_v, uint64_t sample_size) const
{
  // Drop same-query null windows that overlap the candidate interval; other overlaps are kept.
  rsample_v.clear();
  rsample_v.reserve(null_v.size());
  for (const auto& s : null_v) {
    if (s.bix == r.bix && overlaps_half_open(s.bin_iv, r.bin_iv)) continue;
    //TODO: Do we really need Fisher information, I, at this point (in general)?
    if (!std::isfinite(s.I)) continue;
    const double d = validate_distance(s.d);
    if (!std::isfinite(d)) continue;
    rsample_v.push_back({d, s.I});
  }

  if (rsample_v.size() <= sample_size) return;

  // Reservoir sample down to sample_size (algorithm R).
  size_t w = 0;
  for (size_t i = 0; i < rsample_v.size(); ++i) {
    if (w < sample_size) {
      rsample_v[w++] = rsample_v[i];
    } else {
      const size_t j = std::uniform_int_distribution<size_t>(0, i)(gen);
      if (j < sample_size) rsample_v[j] = rsample_v[i];
    }
  }
  rsample_v.resize(static_cast<size_t>(sample_size));
}

template<typename T>
void QIE<T>::test_significance(const uint64_t sample_size)
{
  vec<xy_t> rsample_v;
  for (auto& r : records_v) {
    if (r.is_intact()) {
      warn_pmsg(qid_batch[r.bix], "has a full length interval; skipping significance test");
      continue;
    }

    filter_sample(r, rsample_v, sample_size);

    if (rsample_v.size() < GammaModel::min_nsamples) {
      warn_pmsg(qid_batch[r.bix], "not enough null samples; skipping significance test");
      continue;
    }
    const auto [prob, median] = score_gamma(r, rsample_v);
    if (!std::isfinite(prob) || !std::isfinite(median)) {
      warn_pmsg(qid_batch[r.bix], "gamma fit failed; skipping significance test");
      continue;
    }
    const bool two_sided = std::isfinite(r.d_diff) && (r.is_rc == (r.d_diff > 0.0));
    r.percentile = two_sided ? (2.0 * std::min(prob, 1.0 - prob)) : prob;
    if (std::isfinite(r.d) && median > eps) r.fold = r.d / median;
  }

  // TODO: Check if NaN p-values are handled properly
  benjamini_hochberg_correction();
}

template<typename T>
xy_t QIE<T>::score_gamma(const record_t& r, const vec<xy_t>& samples_v) const
{ // TODO: NaNs and out of bound distances are not handled well
  const double d_obs = validate_distance(r.d);
  vec<double> d_v;
  d_v.reserve(samples_v.size());
  for (const auto& s : samples_v) {
    if (const double d = validate_distance(s.first); std::isfinite(d)) {
      d_v.push_back(d);
    }
  }
  auto result = GammaModel::score_from_samples(d_obs, d_v, d_eps, d_ub - d_eps);
  if (std::isfinite(result.second)) result.second = validate_distance(result.second);
  return result;
}

// Benjamini-Hochberg correction per strand (fw and rc separately).
template<typename T>
void QIE<T>::benjamini_hochberg_correction()
{
  for (bool is_rc : {false, true}) {
    vec<size_t> idx;
    for (size_t i = 0; i < records_v.size(); ++i)
      if (records_v[i].is_rc == is_rc && std::isfinite(records_v[i].percentile)) idx.push_back(i);
    if (idx.empty()) continue;

    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return records_v[a].percentile < records_v[b].percentile; });

    const double m = static_cast<double>(idx.size());
    double q_min = 1.0;
    for (size_t rank = idx.size(); rank >= 1; --rank) {
      record_t& r = records_v[idx[rank - 1]];
      q_min = std::min(q_min, std::min(1.0, r.percentile * m / static_cast<double>(rank)));
      r.qvalue = q_min;
    }
  }
}

template<typename T>
void QIE<T>::report_contiguous(std::ostream& sout, const str& rid) const
{ // Revisit this
  for (const auto& r : records_v) {
    const uint8_t mask = (r.th_ix != size_t(-1)) ? static_cast<uint8_t>(1u << r.th_ix) : 0;
    xy_t d_range{d_eps, d_ub};
    if (r.th_ix != size_t(-1)) {
      const double ext = at(llhf->get_extrema(), r.th_ix);
      if (ext > 0.0)
        d_range.second = std::min(d_range.second, ext);
      else
        d_range.first = std::max(d_range.first, -ext);
    }

    std::ostringstream d_bin;
    d_bin.flags(sout.flags());
    d_bin.precision(sout.precision());
    d_bin << '(' << d_range.first << ", " << d_range.second << ')';

    write_tsv(sout,
              qid_batch[r.bix],
              r.L,
              r.seq_iv.a,
              r.seq_iv.b,
              report_strand(r.is_rc, r.d_diff),
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

template class QIE<double>;
template class QIE<cm512_t>;

template class DIM<double>;
template class DIM<cm512_t>;

template class LLH<double>;
template class LLH<cm512_t>;
