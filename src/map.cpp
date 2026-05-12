#include "map.hpp"
#include "gamma.hpp"
#include "random.hpp"
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
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  mask_bp = u64m >> ((32 - k) * 2);
  if (!params.enum_only) {
    v_acc.assign(hdist_bound + 1, 0); // For SIMD alignment, set to 8
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
      warn_pmsg(qid_batch[bix], "minimum length is exceeded; using the full contig as the effective minimum ");
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
        enumerate_intervals(sout, rid, dim_fw, false);
        enumerate_intervals(sout, rid, dim_rc, true);
      } else {
        for (size_t i = 0; i < WIDTH; ++i) {
          enumerate_intervals(sout, rid, dim_fw, false, i);
          enumerate_intervals(sout, rid, dim_rc, true, i);
        }
      }
    } else {
      const auto [pos_bv, neg_bv] = llhf->get_sign_bv();
      const uint8_t th_bv = pos_bv | neg_bv; // This is always 11111111
      for (auto* dim : {&dim_fw, &dim_rc}) {
        dim->compute_prefhistsum();
      }

      // Extract per-query histograms, compute per-query MLE, accumulate into genome-wide histogram
      vec<uint64_t> v_q_fw, v_q_rc;
      uint64_t u_q_fw = 0, u_q_rc = 0, t_q_fw = 0, t_q_rc = 0;
      dim_fw.extract_histogram(0, nbins, v_q_fw, u_q_fw, t_q_fw);
      dim_rc.extract_histogram(0, nbins, v_q_rc, u_q_rc, t_q_rc);
      const double d_q_fw = compute_mle_dist(v_q_fw, u_q_fw, t_q_fw);
      const double d_q_rc = compute_mle_dist(v_q_rc, u_q_rc, t_q_rc);

      // Reference is the one with the lower distance among fw/rc strands.
      const bool is_ref = std::isnan(d_q_rc) || (!std::isnan(d_q_fw) && d_q_fw <= d_q_rc);
      save_mesh(is_ref ? dim_fw : dim_rc);

      collect_contiguous(dim_fw, th_bv, d_q_fw, false, is_ref);
      collect_contiguous(dim_rc, th_bv, d_q_rc, true, !is_ref);
      add_to_acc(v_acc, u_acc, is_ref ? v_q_fw : v_q_rc, is_ref ? u_q_fw : u_q_rc);
    }
  }

  if (!params.enum_only) {
    uint64_t t_acc = 0;
    for (auto c : v_acc) {
      t_acc += c;
    }
    d_acc = compute_mle_dist(v_acc, u_acc, t_acc);
    test_significance(params.nsamples, 250);
    report_contiguous(sout, rid);
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
    // Note that aggregate_mer() accumulates into rows 1 to nbins
    // At the end, compute_prefhistsum() converts in-place
    hdisthist_v.assign((nbins + 1) * (params.hdist_th + 1), 0);
    // First params.hdist_th values are zeros, same layout as fdps_v/sdps_v
  }
}

template<typename T>
void DIM<T>::aggregate_mer(uint32_t hdist_min, uint64_t i)
{
  // The bin index is i, multiple k-mers in the same bin accumulate here
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
interval_t DIM<T>::get_interval(uint64_t i, size_t ix) const
{
  if ((ix < eintervals_v.size()) && (i < eintervals_v[ix].size())) {
    return eintervals_v[ix][i];
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

template<typename T>
void DIM<T>::extract_intervals_mx(const uint64_t tau, const size_t ix)
{
  uint64_t b_curr = 1;
  uint64_t b_prev = std::numeric_limits<uint64_t>::max(); // no interval yet

  if (nbins >= 1 + tau && at(fdps_v[nbins], ix) < at(fdps_v[1], ix)) {
    rintervals_v[ix].emplace_back(1, nbins);
    return;
  }

  for (uint64_t a = 1; a <= nbins; ++a) {
    const double fdpmax_a = at(fdpmax_v[a - 1], ix);
    const double fdps_a = at(fdps_v[a], ix);

    if (fdpmax_a >= fdps_a) {
      continue; // Skip if a is not a prefix maxima of the prefix sum
    }

    while ((b_curr + 1) <= nbins && (at(fdsmin_v[b_curr + 1], ix) < fdps_a)) {
      ++b_curr; // Right maximal
    } // We increment b_curr to the last b with fdsmin_v[b] < fdps_a
    // There is no other a*>a where b*<b, so no valid (a, b) is missed

    const uint64_t b_star = b_curr;
    if (b_star < (a + tau)) {
      continue; // Skip if no valid right endpoint in [a+tau, nbins]
    }
    if (at(fdps_v[b_star], ix) >= fdps_a) {
      continue; // Skip if negative-sum check fails
    }
    // if (__builtin_expect(b_star == b_prev, 0)) {
    if (b_star == b_prev) {
      continue; // An earlier (leftmost) a already claimed this b*
    }

    if (at(fdps_v[b_star], ix) >= fdpmax_a) {   // Left maximal
      rintervals_v[ix].emplace_back(a, b_star); // 1-based inclusive coordinates
      b_prev = b_star;
    }
  }
}

template<typename T>
void DIM<T>::extract_intervals_sx(const uint64_t tau, const size_t ix)
{
  // Every valid right endpoint b* is a suffix minimum of fdps_v
  // Suffix minimum values are strictly increasing left-to-right,
  // Hence, the pointer into the list is monotone across record highs which is O(k) total
  if (nbins >= 1 + tau && at(fdps_v[nbins], ix) < at(fdps_v[1], ix)) {
    rintervals_v[ix].emplace_back(1, nbins);
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
    for (uint64_t j = nbins; j >= 1; --j) {
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

  for (uint64_t a = 1; a <= nbins; ++a) {
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

    if (fdps_bstar >= fdpmax_a) {               // Left maximal
      rintervals_v[ix].emplace_back(a, b_star); // 1-based inclusive coordinates
      b_prev = b_star;
    }
  }
}

template<typename T>
void DIM<T>::expand_intervals(const double chisq_th, const size_t ix)
{
  if (rintervals_v[ix].empty()) {
    return;
  }

  double fdiff, sdiff, chisq_val;
  uint64_t a, ap, b, bp;
  ap = rintervals_v[ix][0].first;
  bp = rintervals_v[ix][0].second;

  for (uint64_t i = 1; i < rintervals_v[ix].size(); ++i) {
    a = rintervals_v[ix][i].first;
    b = rintervals_v[ix][i].second;
    fdiff = at(fdps_v[b], ix) - at(fdps_v[ap], ix);
    sdiff = at(sdps_v[ap], ix) - at(sdps_v[b], ix);
    // chisq_val = (sdiff > 0.0) ? (fdiff * fdiff) / sdiff : std::numeric_limits<double>::infinity();
    chisq_val = ((fdiff * fdiff) + eps) / (sdiff + eps);

    if ((chisq_val < chisq_th) && (a < bp)) {
      a = ap;
      // b = std::max(bp, b); // This is not necessary due to maximality and monotonicity of a's
    } else {
      eintervals_v[ix].emplace_back(ap, bp); // 1-based inclusive coordinates
    }

    ap = a;
    bp = b;
  }

  fdiff = at(fdps_v[bp], ix) - at(fdps_v[ap], ix);
  sdiff = at(sdps_v[ap], ix) - at(sdps_v[bp], ix);
  // chisq_val = (sdiff > 0.0) ? (fdiff * fdiff) / sdiff : std::numeric_limits<double>::infinity();
  chisq_val = ((fdiff * fdiff) + eps) / (sdiff + eps);
  eintervals_v[ix].emplace_back(ap, bp); // 1-based inclusive coordinates
  rintervals_v[ix].clear();
  rintervals_v[ix].shrink_to_fit();
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
  u = (mers_b - mers_a) - t;
}

template<typename T>
void QIE<T>::save_mesh(const DIM<T>& dim)
{
  static constexpr size_t npoints_min = 8;
  static constexpr size_t npoints_max = 128;

  const uint32_t W = params.hdist_th + 1;
  const uint64_t nbins_q = dim.get_nbins();
  const auto& hist = dim.get_hdisthist();
  const size_t s = static_cast<size_t>(std::sqrt(static_cast<double>(nbins_q)));
  const size_t npoints = std::clamp(s, npoints_min, npoints_max);

  mesh_t m(bix, nbins_q, dim.get_nmers());

  if (nbins_q <= npoints + 1) {
    // Short contig: every boundary 0, 1, ..., nbins_q
    m.points_v.resize(nbins_q + 1);
    std::iota(m.points_v.begin(), m.points_v.end(), uint64_t(0));
  } else {
    // Long contig: sample npoints interior boundaries from {1, ..., nbins_q-1}.
    m.points_v.reserve(npoints + 2);
    m.points_v.push_back(0);
    for (uint64_t j = nbins_q - npoints; j < nbins_q; ++j) {
      const uint64_t r = std::uniform_int_distribution<uint64_t>(1, j)(gen);
      const auto start = m.points_v.begin() + 1;
      const auto it = std::lower_bound(start, m.points_v.end(), r);
      if (it == m.points_v.end() || *it != r)
        m.points_v.insert(it, r);
      else
        m.points_v.push_back(j);
    }
    // Interior entries are sorted distinct values in {1, ..., nbins_q-1}
    m.points_v.push_back(nbins_q);
  }

  const size_t P = m.points_v.size();
  m.hists_v.resize(P * W);
  for (size_t i = 0; i < P; ++i) {
    const uint64_t p = m.points_v[i];
    std::copy(hist.begin() + p * W, hist.begin() + (p + 1) * W, m.hists_v.begin() + i * W);
  }

  meshes_v.push_back(std::move(m));
}

template<typename T>
bool QIE<T>::sample_mesh_distance(const uint64_t mix, const uint64_t bix, const interval_t& ab)
{
  const auto& m = meshes_v[mix];
  const uint32_t W = params.hdist_th + 1;
  vec<uint64_t> v(hdist_bound + 1);

  // Uniform pair of distinct point indices
  std::uniform_int_distribution<uint64_t> rvpt(0, m.points_v.size() - 1);
  uint64_t i = rvpt(gen), j = rvpt(gen);
  if (i == j) {
    return false;
  }
  if (i > j) {
    std::swap(i, j);
  }

  // Reject any window that overlaps with the current interval
  // Samples from other contigs are independent by construction
  if (m.bix == bix && m.points_v[i] < ab.second && ab.first < m.points_v[j]) {
    return false;
  }

  // std::fill(v.begin(), v.end(), 0);
  uint64_t t = 0;
  for (uint32_t x = 0; x < W; ++x) {
    v[x] = m.hists_v[j * W + x] - m.hists_v[i * W + x];
    t += v[x];
  }
  const uint64_t mers_a = std::min(m.points_v[i] << params.bin_shift, m.nmers);
  const uint64_t mers_b = std::min(m.points_v[j] << params.bin_shift, m.nmers);
  const uint64_t u = (mers_b - mers_a) - t;

  const double d = compute_mle_dist(v, u, t);
  const double I = llhf->compute_fisher_info(d);
  if ((std::isfinite(I) && std::isfinite(d)) && (d > 0.0 && I > 0.0)) {
    samples_v.push_back({d, I});
    return true;
  } else {
    return false;
  }
}

template<typename T>
bool QIE<T>::sample_metropolis_hastings(const xy_t& init, uint64_t S, uint64_t B)
{
  bool is_valid = true;
  auto f = [&](const double& D) { return (*llhf)(D); };
  double d = init.first;
  const double step = init.second; // Optimal MH step is sigma
  double llh = f(d);

  std::uniform_real_distribution<double> ruv(0, 1);

  samples_v.reserve(samples_v.size() + S);

  for (uint64_t iter = 0; iter < B + S; ++iter) {
    // Reflect proposal off boundaries to stay in (0,1)
    double d_i = d + step * sample_box_muller(gen);
    // Reflect instead of reject - better mixing near boundaries
    if (d_i <= 0.0) d_i = -d_i;
    if (d_i >= 0.5) d_i = 1.0 - d_i;
    d_i = std::clamp(d_i, eps, 0.5 - eps);

    const double llh_i = f(d_i);
    const double log_alpha = llh_i - llh;
    if (std::log(ruv(gen)) < log_alpha) {
      d = d_i;
      llh = llh_i;
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
void QIE<T>::collect_contiguous(DIM<T>& dim, uint8_t th_bv, double d_q, bool is_rc, bool is_ref)
{
  if (th_bv == 0) return;

  const uint64_t nbins = dim.get_nbins();

  // Collect 1-based half-open bin breakpoints from all active thresholds.
  // Note that, eintervals_v stores 1-based inclusive bin intervals.
  vec<uint64_t> pts = {1, nbins + 1};
  for (size_t ti = 0; ti < WIDTH; ++ti) {
    if (!(th_bv & (1u << ti)) || dim.get_eintervals(ti).empty()) continue;
    for (const auto& iv : dim.get_eintervals(ti)) {
      // Convert to 1-based [first, second+1) bin-boundary coordinates.
      pts.push_back(iv.first);
      pts.push_back(iv.second + 1);
    }
  }

  std::sort(pts.begin(), pts.end());
  pts.erase(std::unique(pts.begin(), pts.end()), pts.end());

  // Scan pointer per threshold into its sorted interval list
  arr<size_t, WIDTH> ti_ix = {};

  const uint64_t L = enmers + k - 1;

  for (size_t pi = 0; pi + 1 < pts.size(); ++pi) {
    const uint64_t a_bin = pts[pi];
    const uint64_t b_bin = pts[pi + 1];
    const bool is_lpi = (pi + 2 == pts.size());

    // Determine which thresholds cover [a_bin, b_bin) and tighten the distance bin bounds
    // Positive extrema implies upper bound (<)
    // Negative extrema implies lower bound (>)
    uint8_t mask = 0;
    xy_t d_range{0.0, 1.0};
    for (size_t ti = 0; ti < WIDTH; ++ti) {
      if (!(th_bv & (1u << ti))) continue;
      const auto& e_v = dim.get_eintervals(ti);
      while (ti_ix[ti] < e_v.size() && e_v[ti_ix[ti]].second < a_bin)
        ++ti_ix[ti];
      if (ti_ix[ti] >= e_v.size() || e_v[ti_ix[ti]].first > a_bin) continue;
      mask |= static_cast<uint8_t>(1u << ti);
      const double ext = at(llhf->get_extrema(), ti);
      if (ext > 0.0)
        d_range.second = std::min(d_range.second, ext);
      else
        d_range.first = std::max(d_range.first, -ext);
    }

    // Convert to 0-based half-open bin-boundary coordinates for histogram extraction.
    vec<uint64_t> v;
    uint64_t u, t;
    dim.extract_histogram(a_bin - 1, b_bin - 1, v, u, t);
    const double d = compute_mle_dist(v, u, t);
    const double I = llhf->compute_fisher_info(d);

    const auto [a, b] = get_rinterval({a_bin, b_bin}, params.bin_shift, enmers, k, is_lpi);
    assert(a < b);
    records_v.emplace_back(bix, L, a, b, a_bin, b_bin, is_rc, is_ref, d, d_q, I, mask, d_range);
  }
}

template<typename T>
double QIE<T>::compute_mle_dist(const vec<uint64_t>& v, uint64_t u, uint64_t t)
{
  llhf->set_counts(v.data(), u);
  auto f = [&](const double& D) { return (*llhf)(D); };
  // const double ub = (t == 0) ? 0.75 + eps : 0.5;
  const double ub = 0.5 - eps;
  xy_t result = boost::math::tools::brent_find_minima(f, eps, ub, 24);
  if (std::isnan(result.first)) result.first = ub; // TODO: How to handle these?
  return result.first;
}

template<typename T>
void QIE<T>::test_significance(const uint64_t nsamples, const uint64_t ntries)
{
  // Weight each query by its bin count to sample proportionally.
  // Null summaries with fewer than 2 points cannot yield any interval.
  vec<uint64_t> weights(meshes_v.size());
  uint64_t total_weight = 0;
  for (size_t qi = 0; qi < meshes_v.size(); ++qi) {
    if (meshes_v[qi].points_v.size() < 2) continue;
    weights[qi] = meshes_v[qi].nbins;
    total_weight += weights[qi];
  }
  if (total_weight == 0) {
    warn_msg("No valid meshes found; skipping significance test");
    return;
  }
  std::discrete_distribution<size_t> rvix(weights.begin(), weights.end());
  samples_v.reserve(nsamples);

  for (auto& r : records_v) {
    samples_v.clear();
    uint64_t t = 0;
    uint64_t n = 0;
    const interval_t ab = r.get_interval();
    if (r.is_intact()) {
      warn_pmsg(qid_batch[r.bix], "has a full length interval; skipping significance test");
      continue;
    }
    while ((n < nsamples) && (t < ntries)) {
      ++t;
      const uint64_t mix = rvix(gen);
      if (sample_mesh_distance(mix, r.bix, ab)) {
        const auto [d, I] = samples_v.back();
        const size_t nsprev = samples_v.size() - 1;
        if (sample_metropolis_hastings({d, std::sqrt(1.0 / I)})) {
          ++n;
        } else {
          samples_v.resize(nsprev);
          warn_pmsg(qid_batch[r.bix], "metropolis-hastings returned an invalid sample; skipping significance test");
        }
      } else {
        warn_pmsg(qid_batch[r.bix], "mesh sampling returned an invalid sample; skipping significance test");
      }
    }

    if (samples_v.size() < (nsamples / 2)) {
      warn_pmsg(qid_batch[r.bix], "not enough non-overlapping meshes; skipping significance test");
      continue;
    }
    const auto [prob, median] = score_gamma(r, samples_v);
    if (std::isnan(prob)) {
      warn_pmsg(qid_batch[r.bix], "gamma fit failed; skipping significance test");
      continue;
    }
    // The null-sampled strand gets a two-sided test; the opposite strand keeps the low-distance tail.
    r.percentile = r.is_ref ? (2.0 * std::min(prob, 1.0 - prob)) : prob;
    r.fold = (median > eps) ? (r.d / median) : std::numeric_limits<double>::quiet_NaN();
  }
}

template<typename T>
xy_t QIE<T>::score_gamma(const record_t& r, const vec<xy_t>& samples_v) const
{
  // Fit Gamma directly to MCMC posterior draws from reference likelihoods.
  // The draws already represent the latent distribution -- no noise deconvolution needed.
  // Test record's noise sigma_r is accounted for in the marginal_cdf evaluation below.
  vec<double> d_v;
  for (const auto& s : samples_v) {
    d_v.push_back(s.first);
  }
  const GammaModel::params_t gp = GammaModel::fit_from_samples(d_v);
  if (!(gp.shape > 0.0) || !(gp.scale > 0.0)) {
    return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
  }
  const double sigma_r = (r.I > 0.0 && std::isfinite(r.I)) ? (1.0 / std::sqrt(r.I)) : 0.0;
  const double prob = GammaModel::gamma_cdf(r.d, gp.shape, gp.scale);

  // Median of the latent Gamma distribution by binary search on [eps, 1-eps].
  // Uses plain gamma_cdf since the fit was on latent draws -- no convolution with sigma_r.
  auto F = [&](double x) { return GammaModel::gamma_cdf(x, gp.shape, gp.scale); };

  double d_low = GammaModel::eps, d_high = 1.0 - GammaModel::eps;
  for (int it = 0; it < 40; ++it) {
    const double d_m = 0.5 * (d_low + d_high);
    if (F(d_m) < 0.5) {
      d_low = d_m;
    } else {
      d_high = d_m;
    }
  }
  const double median = 0.5 * (d_low + d_high);
  return {prob, median};
}

template<typename T>
void QIE<T>::report_contiguous(std::ostream& sout, const str& rid) const
{ // TODO: Revisit?
  for (const auto& r : records_v) {
    const char strand = r.is_ref ? '-' : '+';
    std::ostringstream d_bin;
    d_bin.flags(sout.flags());
    d_bin.precision(sout.precision());
    d_bin << '(' << r.d_low << ", " << r.d_high << ')';
    write_tsv(sout,
              qid_batch[r.bix],
              r.L,
              r.a,
              r.b,
              strand,
              rid,
              r.d,
              static_cast<uint32_t>(r.mask),
              d_bin.str(),
              r.d_q,
              d_acc,
              r.percentile,
              r.fold)
      << '\n';
  }
}

template<typename T>
void QIE<T>::enumerate_intervals(std::ostream& sout, const str& rid, DIM<T>& dim, bool is_rc, size_t ix)
{ // TODO: Revisit?
  const char strand = is_rc ? '-' : '+';
  const double dist_th = at(params.dist_th, ix);
  const uint64_t nbins = dim.get_nbins();
  const uint64_t L = enmers + k - 1;
  for (uint64_t i = 0;; ++i) {
    const interval_t ab = dim.get_interval(i, ix);
    if (ab.first >= nbins) break;
    const auto [a, b] = get_einterval(ab, params.bin_shift, enmers, k);
    assert(a < b);
    write_tsv(sout, qid_batch[bix], L, a, b, strand, rid, dist_th) << '\n';
  }
}

template class QIE<double>;
template class QIE<cm512_t>;

template class DIM<double>;
template class DIM<cm512_t>;

template class LLH<double>;
template class LLH<cm512_t>;
