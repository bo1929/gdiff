#include "map.hpp"

template<typename T>
DIM<T>::DIM(llh_sptr_t<T> llhf, uint32_t hdist_th, uint64_t en_mers)
  : llhf(llhf)
  , hdist_th(hdist_th)
  , en_mers(en_mers)
{
  hdisthist_v.resize(hdist_th + 1, 0);

  // fdc_v.reserve(en_mers);
  // sdc_v.reserve(en_mers);
  fdc_v.resize(en_mers);
  sdc_v.resize(en_mers);
  fdps_v.reserve(en_mers + 1);
  sdps_v.reserve(en_mers + 1);
  fdpmax_v.reserve(en_mers + 2);
  fdsmin_v.reserve(en_mers + 2);

  if constexpr (std::is_same_v<T, double>) {
    fdt = 0.0;
    sdt = 0.0;
  } else {
    fdt.fill(0.0);
    sdt.fill(0.0);
  }
  // for (size_t i = 0; i < WIDTH; ++i) {
  //   rintervals_v[i].clear();
  //   eintervals_v[i].clear();
  //   chisq_v[i].clear();
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
  const uint32_t hdist_min = sketch->search_mer(rix, enc_lr);

  if (hdist_min <= hdist_th) {
    merhit_count++;
    hdisthist_v[hdist_min]++;
    // sdc_v.push_back(llhf->get_sdc(hdist_min));
    // fdc_v.push_back(llhf->get_fdc(hdist_min));
    sdc_v[i] = llhf->get_sdc(hdist_min);
    fdc_v[i] = llhf->get_fdc(hdist_min);
  } else {
    mermiss_count++;
    // sdc_v.push_back(llhf->get_sdc());
    // fdc_v.push_back(llhf->get_fdc());
    sdc_v[i] = llhf->get_sdc();
    fdc_v[i] = llhf->get_fdc();
  }

  if constexpr (std::is_same_v<T, double>) {
    // fdt += fdc_v.back();
    // sdt += sdc_v.back();
    fdt += fdc_v[i];
    sdt += sdc_v[i];
  } else {
    simde__m512d vfdt = simde_mm512_load_pd(fdt.data());
    simde__m512d vsdt = simde_mm512_load_pd(sdt.data());
    simde__m512d fdc = simde_mm512_load_pd(fdc_v[i].data());
    simde__m512d sdc = simde_mm512_load_pd(sdc_v[i].data());
    vfdt = simde_mm512_add_pd(vfdt, fdc);
    vsdt = simde_mm512_add_pd(vsdt, sdc);
    simde_mm512_store_pd(fdt.data(), vfdt);
    simde_mm512_store_pd(sdt.data(), vsdt);
  }
}

template<typename T>
void DIM<T>::inclusive_scan()
{
  assert(en_mers > 0);
  const uint64_t s = en_mers + 1;

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
  for (uint64_t a = 1, b = 1; a <= en_mers; ++a) {
    const double fdpmax_a = at(fdpmax_v[a - 1], idx);
    const double fdps_a = at(fdps_v[a], idx);

    if (fdpmax_a >= fdps_a) {
      continue;
    }

    if (b < (a + tau)) {
      b = a + tau - 1;
    }
    if (__builtin_expect(b > en_mers, 0)) {
      break;
    }

    if (at(fdsmin_v[b + 1], idx) >= at(fdps_v[a], idx)) {
      continue;
    }

    b++;
    while (b <= en_mers) {
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
    chisq_val = (fdiff * fdiff) / sdiff;

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
  chisq_val = (fdiff * fdiff) / sdiff;
  eintervals_v[idx].emplace_back(ap, bp);
  chisq_v[idx].push_back(chisq_val);

  return eintervals_v[idx].size();
}

template<typename T>
void DIM<T>::report_intervals(std::ostream& output_stream, const std::string& identifier, const size_t idx)
{ // TODO: Revisit.
  if (eintervals_v[idx].empty()) {
    output_stream << identifier << ',' << en_mers << ',' << en_mers << ',' << en_mers << ",0,0\n";
  } else {
    for (uint64_t i = 0; i < eintervals_v[idx].size(); ++i) {
      const uint64_t a = eintervals_v[idx][i].first;
      const uint64_t b = eintervals_v[idx][i].second;
      output_stream << identifier << ',' << a << ',' << b - 1 << ',' << en_mers << ','
                    << (at(fdps_v[b], idx) - at(fdps_v[a], idx)) << ',' << chisq_v[idx][i] << '\n';
    }
  }
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
QIE<T>::QIE(sketch_sptr_t sketch, lshf_sptr_t lshf, qseq_sptr_t qs, params_t<T> params)
  : sketch(sketch)
  , lshf(lshf)
  , k(lshf->get_k())
  , h(lshf->get_h())
  , m(lshf->get_m())
  , batch_size(qs->cbatch_size)
  , params(params)
  , min_length(params.min_length)
  , chisq(params.chisq)
{
  std::swap(qs->seq_batch, seq_batch);
  std::swap(qs->identifer_batch, identifer_batch);
  llhf = std::make_shared<LLH<T>>(k, h, sketch->get_rho(), params.hdist_th, params.dist_th);
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  mask_bp = u64m >> ((32 - k) * 2);
}

template<typename T>
void QIE<T>::map_sequences(std::ostream& output_stream)
{
  strstream batch_stream;
  // batch_stream.precision(4);
  for (bix = 0; bix < batch_size; ++bix) {
    const char* cseq = seq_batch[bix].data();
    const uint64_t len = seq_batch[bix].size();
    onmers = 0;
    en_mers = len - k + 1;

    DIM<T> dim_or(llhf, params.hdist_th, en_mers);
    DIM<T> dim_rc(llhf, params.hdist_th, en_mers);
    search_mers(cseq, len, dim_or, dim_rc);

    if constexpr (std::is_same_v<T, double>) {
      const bool use_rc = (llhf->get_sign() * dim_or.get_fdt()) >= (llhf->get_sign() * dim_rc.get_fdt());
      DIM<T>& dim = use_rc ? dim_rc : dim_or;
      dim.inclusive_scan();
      dim.extract_intervals(std::min(min_length, en_mers) - 1);
      dim.expand_intervals(chisq);
      dim.report_intervals(batch_stream, identifer_batch[bix]);
    } else {
      // TODO: Do not do this twice?
      dim_or.inclusive_scan();
      dim_rc.inclusive_scan();
      for (size_t i = 0; i < WIDTH; ++i) {
        const bool use_rc = (llhf->get_sign()[i] * dim_or.get_fdt()[i]) >= (llhf->get_sign()[i] * dim_rc.get_fdt()[i]);
        DIM<T>& dim = use_rc ? dim_rc : dim_or;
        dim.extract_intervals(std::min(min_length, en_mers) - 1, i);
        dim.expand_intervals(chisq, i);
        dim.report_intervals(batch_stream, identifer_batch[bix], i);
      }
    }
  }
  output_stream << batch_stream.rdbuf();
}

template<typename T>
void QIE<T>::search_mers(const char* cseq, uint64_t len, DIM<T>& dim_or, DIM<T>& dim_rc)
{
  uint32_t i = 0, j = 0, l = 0;
  uint32_t orrix, rcrix;
  uint64_t orenc64_bp, orenc64_lr, rcenc64_bp;
  for (; i < len; ++i) {
    if (__builtin_expect(SEQ_NT4_TABLE[cseq[i]] >= 4, 0)) {
      l = 0;
      continue;
    }
    ++l;
    if (l < k) {
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
#ifdef CANONICAL
    if (rcenc64_bp < orenc64_bp) {
      orrix = lshf->compute_hash(orenc64_bp);
      if (__builtin_expect(sketch->check_partial(orrix), 1)) {
        dim_or.aggregate_mer(sketch, orrix, lshf->drop_ppos_lr(orenc64_lr), j);
      }
    } else {
      rcrix = lshf->compute_hash(rcenc64_bp);
      if (__builtin_expect(sketch->check_partial(rcrix), 1)) {
        dim_rc.aggregate_mer(sketch, rcrix, lshf->drop_ppos_lr(bp64_to_lr64(rcenc64_bp)), j);
      }
    }
#else
    orrix = lshf->compute_hash(orenc64_bp);
    if (__builtin_expect(sketch->check_partial(orrix), 1)) {
      dim_or.aggregate_mer(sketch, orrix, lshf->drop_ppos_lr(orenc64_lr), j);
    }
    rcrix = lshf->compute_hash(rcenc64_bp);
    if (__builtin_expect(sketch->check_partial(rcrix), 1)) {
      dim_rc.aggregate_mer(sketch, rcrix, lshf->drop_ppos_lr(bp64_to_lr64(rcenc64_bp)), j);
    }
#endif /* CANONICAL */
  }
}

template class QIE<double>;
template class QIE<cm512_t>;

template class DIM<double>;
template class DIM<cm512_t>;

template class LLH<double>;
template class LLH<cm512_t>;
