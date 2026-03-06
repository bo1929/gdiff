#include "map.hpp"

#define EPS 1e-10

template<typename T>
DIM<T>::DIM(llh_sptr_t<T> llhf, uint32_t hdist_th, uint64_t n)
  : llhf(llhf)
  , hdist_th(hdist_th)
  , n(n)
{
  hdisthist_v.resize(hdist_th + 1, 0);

  // fdc_v.reserve(n);
  // sdc_v.reserve(n);
  fdc_v.resize(n);
  sdc_v.resize(n);
  fdps_v.reserve(n + 1);
  sdps_v.reserve(n + 1);
  fdpmax_v.reserve(n + 2);
  fdsmin_v.reserve(n + 2);

  // This is for deciding between fw and rc.
  // if constexpr (std::is_same_v<T, double>) {
  //   fdt = 0.0;
  //   sdt = 0.0;
  // } else {
  //   fdt.fill(0.0);
  //   sdt.fill(0.0);
  // }
}

template<typename T>
inline double DIM<T>::at(T v, const size_t idx)
{
  if constexpr (std::is_same_v<T, double>) {
    return v;
  } else {
    return v[idx];
  }
}

template<typename T>
void DIM<T>::aggregate_mer(sketch_sptr_t sketch, uint32_t rix, enc_t enc_lr, uint64_t i)
{
  // This is binned, and i is the bin index.
  // Multiple k-mers that fall in the same bin all accumulate here.
  const uint32_t hdist_min = sketch->search_mer(rix, enc_lr);

  if (hdist_min <= hdist_th) {
    merhit_count++;
    hdisthist_v[hdist_min]++;
    add_to(sdc_v[i], llhf->get_sdc(hdist_min));
    add_to(fdc_v[i], llhf->get_fdc(hdist_min));
  } else {
    mermiss_count++;
    add_to(sdc_v[i], llhf->get_sdc());
    add_to(fdc_v[i], llhf->get_fdc());
  }

  // This is for deciding between fw and rc.
  // if constexpr (std::is_same_v<T, double>) {
  //   // fdt += fdc_v.back();
  //   // sdt += sdc_v.back();
  //   fdt += fdc_v[i];
  //   sdt += sdc_v[i];
  // } else {
  //   simde__m512d vfdt = simde_mm512_load_pd(fdt.data());
  //   simde__m512d vsdt = simde_mm512_load_pd(sdt.data());
  //   simde__m512d fdc = simde_mm512_load_pd(fdc_v[i].data());
  //   simde__m512d sdc = simde_mm512_load_pd(sdc_v[i].data());
  //   vfdt = simde_mm512_add_pd(vfdt, fdc);
  //   vsdt = simde_mm512_add_pd(vsdt, sdc);
  //   simde_mm512_store_pd(fdt.data(), vfdt);
  //   simde_mm512_store_pd(sdt.data(), vsdt);
  // }
}

template<typename T>
void DIM<T>::inclusive_scan()
{
  assert(n > 0);
  const uint64_t s = n + 1;

  fdps_v.resize(s);
  sdps_v.resize(s);
  fdpmax_v.resize(s + 1);
  fdsmin_v.resize(s + 1);

  if constexpr (std::is_same_v<T, double>) {
    fdps_v[0] = 0.0;
    sdps_v[0] = 0.0;
    for (uint64_t i = 1; i < s; ++i) {
      fdps_v[i] = fdps_v[i - 1] + fdc_v[i - 1];
      sdps_v[i] = sdps_v[i - 1] + sdc_v[i - 1];
    }

    fdpmax_v.front() = -std::numeric_limits<double>::max();
    std::inclusive_scan(
      fdps_v.begin(), fdps_v.end(), fdpmax_v.begin() + 1, [](double a, double b) { return std::max(a, b); });

    fdsmin_v.back() = std::numeric_limits<double>::max();
    std::inclusive_scan(
      fdps_v.rbegin(), fdps_v.rend(), fdsmin_v.rbegin() + 1, [](double a, double b) { return std::min(a, b); });

  } else {
    fdps_v[0].fill(0.0);
    sdps_v[0].fill(0.0);
    simde__m512d fdps_acc = simde_mm512_setzero_pd();
    simde__m512d sdps_acc = simde_mm512_setzero_pd();
    for (uint64_t i = 1; i < s; ++i) {
      const simde__m512d fdc = simde_mm512_load_pd(fdc_v[i - 1].data());
      const simde__m512d sdc = simde_mm512_load_pd(sdc_v[i - 1].data());
      fdps_acc = simde_mm512_add_pd(fdps_acc, fdc);
      sdps_acc = simde_mm512_add_pd(sdps_acc, sdc);
      simde_mm512_store_pd(fdps_v[i].data(), fdps_acc);
      simde_mm512_store_pd(sdps_v[i].data(), sdps_acc);
    }

    fdpmax_v[0].fill(-std::numeric_limits<double>::max());
    simde__m512d fdpmax_acc = simde_mm512_load_pd(fdpmax_v[0].data());
    for (uint64_t i = 0; i < s; ++i) {
      const simde__m512d fdps = simde_mm512_load_pd(fdps_v[i].data());
      fdpmax_acc = simde_mm512_max_pd(fdpmax_acc, fdps);
      simde_mm512_store_pd(fdpmax_v[i + 1].data(), fdpmax_acc);
    }

    uint64_t j = s;
    fdsmin_v[j].fill(std::numeric_limits<double>::max());
    simde__m512d fdsmin_acc = simde_mm512_load_pd(fdsmin_v[j].data());
    for (uint64_t i = 0; i < s; ++i) {
      j = s - i - 1;
      const simde__m512d fdps = simde_mm512_load_pd(fdps_v[j].data());
      fdsmin_acc = simde_mm512_min_pd(fdsmin_acc, fdps);
      simde_mm512_store_pd(fdsmin_v[j].data(), fdsmin_acc);
    }
  }
}

template<typename T>
void DIM<T>::extract_intervals(const uint64_t tau, const size_t idx)
{
  for (uint64_t a = 1, b = 1; a <= n; ++a) {
    const double fdpmax_a = at(fdpmax_v[a - 1], idx);
    const double fdps_a = at(fdps_v[a], idx);

    if (fdpmax_a >= fdps_a) {
      continue;
    }

    if (b < (a + tau)) {
      b = a + tau - 1;
    }
    if (__builtin_expect(b > n, 0)) {
      break;
    }

    if (at(fdsmin_v[b + 1], idx) >= at(fdps_v[a], idx)) {
      continue;
    }

    b++;
    while (b <= n) {
      const double fdps_b = at(fdps_v[b], idx);
      const bool negative_sum = fdps_b < fdps_a;
      const bool left_maximal = fdpmax_a <= fdps_b;
      const bool right_maximal = fdps_a <= at(fdsmin_v[b + 1], idx);
      if (negative_sum && left_maximal && right_maximal) {
        rintervals_v[idx].emplace_back(a, b);
        break;
      }
      b++;
    }
  }
}

template<typename T>
uint64_t DIM<T>::expand_intervals(const double chisq_th, const size_t idx)
{
  if (rintervals_v[idx].empty()) {
    return 0;
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
    chisq_val = (fdiff * fdiff) / (sdiff + EPS);

    if ((chisq_val < chisq_th) && (a < bp)) {
      a = ap;
    } else {
      eintervals_v[idx].emplace_back(ap, bp);
      chisq_v[idx].push_back(chisq_val);
    }

    ap = a;
    bp = b;
  }

  fdiff = at(fdps_v[bp], idx) - at(fdps_v[ap], idx);
  sdiff = at(sdps_v[ap], idx) - at(sdps_v[bp], idx);
  // chisq_val = (sdiff > 0.0) ? (fdiff * fdiff) / sdiff : std::numeric_limits<double>::infinity();
  chisq_val = (fdiff * fdiff) / (sdiff + EPS);
  eintervals_v[idx].emplace_back(ap, bp);
  chisq_v[idx].push_back(chisq_val);

  return eintervals_v[idx].size();
}

// template<typename T>
// void DIM<T>::optimize_loglikelihood()
// {
//   llhf->set_counts(hdisthist_v.data(), mermiss_count);
//   /* std::pair<double, double> sol_r = boost::math::tools::brent_find_minima((*llhf), 1e-10, 0.5, 16); */
//   /* d_llh = sol_r.first; */
//   /* v_llh = sol_r.second; */
// }

template<typename T>
interval_t DIM<T>::get_interval(uint64_t i, size_t idx)
{
  if ((idx < eintervals_v.size()) && (i < eintervals_v[idx].size())) {
    return eintervals_v[idx][i];
  } else {
    return {n, n};
  }
}

template<typename T>
QIE<T>::QIE(sketch_sptr_t sketch, lshf_sptr_t lshf, const vec<str>& seq_batch, const vec<str>& qid_batch, params_t<T> params)
  : sketch(sketch)
  , lshf(lshf)
  , k(lshf->get_k())
  , h(lshf->get_h())
  , m(lshf->get_m())
  , batch_size(seq_batch.size())
  , params(params)
  , min_length(params.min_length)
  , bin_shift(params.bin_shift)
  , bin_size(1 << params.bin_shift)
  , chisq(params.chisq)
  , seq_batch(seq_batch)
  , qid_batch(qid_batch)
{
  llhf = std::make_shared<LLH<T>>(k, h, sketch->get_rho(), params.hdist_th, params.dist_th);
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  mask_bp = u64m >> ((32 - k) * 2);
}

template<typename T>
void QIE<T>::map_sequences(std::ostream& sout, const str& rid)
{
  for (bix = 0; bix < batch_size; ++bix) {
    const char* cseq = seq_batch[bix].data();
    const uint64_t len = seq_batch[bix].size();
    onmers = 0;
    if (len < static_cast<uint64_t>(k)) {
      // const uint64_t n = len;
      // if constexpr (std::is_same_v<T, double>) {
      //   sout << WRITE_CINTERVAL(qid_batch[bix], n, n, n, "+", rid, at(params.dist_th, 0)) << '\n';
      // } else {
      //   for (size_t idx = 0; idx < WIDTH; ++idx)
      //     sout << WRITE_CINTERVAL(qid_batch[bix], n, n, n, "+", rid, at(params.dist_th, idx)) << '\n';
      // }
      continue;
    }
    enmers = len - k + 1;
    nbins = (enmers + bin_size - 1) >> bin_shift;

    DIM<T> dim_fw(llhf, params.hdist_th, nbins);
    DIM<T> dim_rc(llhf, params.hdist_th, nbins);
    search_mers(cseq, len, dim_fw, dim_rc);

    // Convert min_length (in base-pairs / k-mer units) to bin units for tau
    const uint64_t tau_bins = (min_length + bin_size - 1) >> bin_shift;

    if constexpr (std::is_same_v<T, double>) {
      dim_fw.inclusive_scan();
      dim_rc.inclusive_scan();
      dim_fw.extract_intervals(std::min(tau_bins, nbins) - 1);
      dim_rc.extract_intervals(std::min(tau_bins, nbins) - 1);
      dim_fw.expand_intervals(chisq);
      dim_rc.expand_intervals(chisq);
      report_intervals(sout, rid, dim_fw, false);
      report_intervals(sout, rid, dim_rc, true);
    } else {
      // TODO: Do not do this twice?
      dim_fw.inclusive_scan();
      dim_rc.inclusive_scan();
      for (size_t i = 0; i < WIDTH; ++i) {
        dim_fw.extract_intervals(std::min(tau_bins, nbins) - 1, i);
        dim_rc.extract_intervals(std::min(tau_bins, nbins) - 1, i);
        dim_fw.expand_intervals(chisq, i);
        dim_rc.expand_intervals(chisq, i);
        report_intervals(sout, rid, dim_fw, false, i);
        report_intervals(sout, rid, dim_rc, true, i);
      }
    }
  }
}

template<typename T>
void QIE<T>::search_mers(const char* cseq, uint64_t len, DIM<T>& dim_fw, DIM<T>& dim_rc)
{
  uint32_t i = 0, j = 0, l = 0;
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
    const uint64_t bix_j = j >> bin_shift; // The bin index for this k-mer position
#ifdef CANONICAL
    if (rcenc64_bp < orenc64_bp) {
      orrix = lshf->compute_hash(orenc64_bp);
      if (__builtin_expect(sketch->check_partial(orrix), 1)) {
        dim_fw.aggregate_mer(sketch, orrix, lshf->drop_ppos_lr(orenc64_lr), bix_j);
      }
    } else {
      rcrix = lshf->compute_hash(rcenc64_bp);
      if (__builtin_expect(sketch->check_partial(rcrix), 1)) {
        dim_rc.aggregate_mer(sketch, rcrix, lshf->drop_ppos_lr(bp64_to_lr64(rcenc64_bp)), bix_j);
      }
    }
#else
    orrix = lshf->compute_hash(orenc64_bp);
    if (__builtin_expect(sketch->check_partial(orrix), 1)) {
      dim_fw.aggregate_mer(sketch, orrix, lshf->drop_ppos_lr(orenc64_lr), bix_j);
    }
    rcrix = lshf->compute_hash(rcenc64_bp);
    if (__builtin_expect(sketch->check_partial(rcrix), 1)) {
      dim_rc.aggregate_mer(sketch, rcrix, lshf->drop_ppos_lr(bp64_to_lr64(rcenc64_bp)), bix_j);
    }
#endif /* CANONICAL */
  }
}

template<typename T>
void QIE<T>::report_intervals(std::ostream& sout, const str& rid, DIM<T>& dim, bool rc, size_t idx)
{ // TODO: Revisit?
  const str strand = rc ? "-" : "+";
  double dist_th = at(params.dist_th, idx);
  uint64_t n = enmers + k - 1;
  interval_t x;
  uint64_t i = 0;
  x = dim.get_interval(i, idx);
  while (x.first < nbins) {
    const uint64_t start = x.first << bin_shift;
    const uint64_t end = std::min(((x.second + 1) << bin_shift) - 1, enmers - 1) + k - 1;
    sout << WRITE_CINTERVAL(qid_batch[bix], n, start, end, strand, rid, dist_th) << '\n';
    x = dim.get_interval(++i, idx);
  }
}

template<typename T>
inline double QIE<T>::at(T v, const size_t idx)
{
  if constexpr (std::is_same_v<T, double>) {
    return v;
  } else {
    return v[idx];
  }
}

template class QIE<double>;
template class QIE<cm512_t>;

template class DIM<double>;
template class DIM<cm512_t>;

template class LLH<double>;
template class LLH<cm512_t>;
