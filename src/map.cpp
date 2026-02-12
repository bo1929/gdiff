#include "map.hpp"

using std::make_pair;
using std::make_tuple;
using std::max;
using std::min;
using std::pair;
using std::tuple;

template<typename T>
std::ostream& operator<<(std::ostream& os, const vec<T>& v)
{
  for (const auto& element : v) {
    os << element << " ";
  }
  return os;
}

vec<double> compute_prefix_sum(const vec<double>& X)
{
  int n = X.size();
  double alpha = 1.0 / 1000;
  vec<double> P(n + 1, 0.0);
  P[1] = P[0] + X[0];
  for (int i = 1; i < n; ++i) {
    P[i + 1] = (P[i] + X[i]);
  }
  // for (int i = 2; i < n; ++i) {
  //   P[i] = P[i] * alpha + P[i - 1] * (1 - alpha);
  // }
  // for (int i = 0; i < n; ++i) {
  //   P[i] = std::lround(P[i] / 2) - std::lround(P[i] / 2) % 2;
  //   // std::cout << P[i + 1] << std::endl;
  // }
  return P;
}

vec<double> prefix_maxima(const vec<double>& P)
{
  int n = P.size();
  vec<double> M(n);
  double running_max = -std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i) {
    running_max = max(running_max, P[i]);
    M[i] = running_max;
  }
  return M;
}

vec<double> suffix_min(const vec<double>& P)
{
  int n = P.size();
  vec<double> S_min(n + 1, std::numeric_limits<double>::infinity());
  for (int i = n - 1; i >= 0; --i)
    S_min[i] = min(P[i], S_min[i + 1]);
  return S_min;
}

SBatch::SBatch(sketch_sptr_t sketch,
               qseq_sptr_t qs,
               uint32_t hdist_th,
               double dist_th,
               uint64_t min_length,
               double chi_sq,
               bool divergent)
  : sketch(sketch)
  , hdist_th(hdist_th)
  , dist_th(dist_th)
  , min_length(min_length)
  , chi_sq(chi_sq)
  , divergent(divergent)
{
  lshf = sketch->get_lshf();
  k = lshf->get_k();
  h = lshf->get_h();
  m = lshf->get_m();
  batch_size = qs->cbatch_size;
  std::swap(qs->seq_batch, seq_batch);
  std::swap(qs->identifer_batch, identifer_batch);
  llhfunc = optimize::HDistHistLLH(h, k, hdist_th);
  uint64_t u64m = std::numeric_limits<uint64_t>::max();
  mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  mask_bp = u64m >> ((32 - k) * 2);
  rho = sketch->get_rho();
  sd.resize(hdist_th + 2);
  sdd.resize(hdist_th + 2);
  for (uint32_t x = 0; x <= hdist_th; ++x) {
    sd[x] = llhfunc.smd(dist_th, x);
    sdd[x] = llhfunc.smdd(dist_th, x);
  }
  sd[hdist_th + 1] = llhfunc.mpd(dist_th, rho);
  sdd[hdist_th + 1] = llhfunc.mpdd(dist_th, rho);
}

void SBatch::map_sequences(std::ostream& output_stream)
{
  strstream batch_stream;
  for (bix = 0; bix < batch_size; ++bix) {
    const char* seq = seq_batch[bix].data();
    uint64_t len = seq_batch[bix].size();
    enmers = len - k + 1;
    onmers = 0;

    ssum1 = 0;
    csum1 = 0;
    ssum2 = 0;
    csum2 = 0;
    s1.resize(enmers);
    c1.resize(enmers);
    s2.resize(enmers);
    c2.resize(enmers);
    prefix_sum_s.resize(enmers + 1);
    prefix_sum_c.resize(enmers + 1);

    std::fill(s1.begin(), s1.end(), 0);
    std::fill(c1.begin(), c1.end(), 0);
    std::fill(s2.begin(), s2.end(), 0);
    std::fill(c2.begin(), c2.end(), 0);
    std::fill(prefix_sum_s.begin(), prefix_sum_s.end(), 0);
    std::fill(prefix_sum_c.begin(), prefix_sum_c.end(), 0);

    prefix_sum_s[0] = 0;
    prefix_sum_c[0] = 0;

    SSummary or_summary(enmers, hdist_th);
    SSummary rc_summary(enmers, hdist_th);
    search_mers(seq, len, or_summary, rc_summary);

    /* std::cout << ssum1  << "/" << ssum2 << std::endl; */
    if (ssum1 < ssum2) {
      for (int j = 0; j < prefix_sum_c.size() - 1; ++j) {
        prefix_sum_c[j + 1] = prefix_sum_c[j] + c1[j];
        prefix_sum_s[j + 1] = prefix_sum_s[j] + s1[j];
      }
    } else {
      for (int j = 0; j < prefix_sum_c.size() - 1; ++j) {
        prefix_sum_c[j + 1] = prefix_sum_c[j] + c2[j];
        prefix_sum_s[j + 1] = prefix_sum_s[j] + s2[j];
      }
    }

    vec<double> prefmax(prefix_sum_s.size() + 1);
    std::fill(prefmax.begin(), prefmax.end(), 0);
    prefmax[0] = -std::numeric_limits<double>::max();
    std::inclusive_scan(
      prefix_sum_s.begin(), prefix_sum_s.end(), prefmax.begin() + 1, [](double a, double b) { return std::max(a, b); });

    vec<double> suffmin(prefix_sum_s.size());
    std::inclusive_scan(
      prefix_sum_s.rbegin(), prefix_sum_s.rend(), suffmin.rbegin(), [](double a, double b) { return std::min(a, b); });
    suffmin.push_back(std::numeric_limits<double>::max());

    int n = enmers;
    uint64_t tau = std::min(min_length, len) - 100;
    int b = 1;
    int ap = 0;
    int bp = 0;
    vec<int> start_idx;
    vec<int> end_idx;

    for (int a = 1; a <= n; ++a) {
      /* batch_stream << "X:" <<  a << std::endl; */
      if (prefmax[a - 1] >= prefix_sum_s[a]) {
        continue;
      }
      if (b < (a + tau)) {
        b = a + tau - 1;
      }
      if (b > n) {
        break;
      }
      if (suffmin[b + 1] >= prefix_sum_s[a]) {
        /* batch_stream << suffmin[b - 1] << std::endl; */
        /* batch_stream << suffmin[b] << std::endl; */
        /* batch_stream << suffmin[b + 1] << std::endl; */
        /* batch_stream << prefix_sum_s[a - 1] << std::endl; */
        /* batch_stream << prefix_sum_s[a] << std::endl; */
        /* batch_stream << prefix_sum_s[a + 1] << std::endl; */
        continue;
      }
      b++;
      while (b <= n) {
        if ((prefmax[a - 1] <= prefix_sum_s[b]) && (prefix_sum_s[b] < prefix_sum_s[a]) &&
            (prefix_sum_s[a] <= suffmin[b + 1])) {
          /* batch_stream << bix << ixx << identifer_batch[bix] << "," << a << "," << b - 1 << ',' << n << "," */
          /*              << prefix_sum_s[b] - prefix_sum_s[a] << "," << 0 << "\n"; */
          start_idx.push_back(a);
          end_idx.push_back(b);
          break;
        }
        b++;
      }
    }
    if (!start_idx.empty()) {
      ap = start_idx[0];
      bp = end_idx[0];

      for (int i = 1; i < start_idx.size(); ++i) {
        int a = start_idx[i];
        int b = end_idx[i];
        double chi_sq_b =
          (prefix_sum_s[b] - prefix_sum_s[ap]) * (prefix_sum_s[b] - prefix_sum_s[ap]) / (prefix_sum_c[ap] - prefix_sum_c[b]);
        if ((chi_sq_b < chi_sq) && (chi_sq_b > 0) && (a < bp)) {
          a = ap;
        } else {
          batch_stream << identifer_batch[bix] << "," << ap << "," << bp - 1 << ',' << n << ","
                       << prefix_sum_s[bp] - prefix_sum_s[ap] << "," << chi_sq_b << "\n";
        }
        ap = a;
        bp = b;
      }
      double chi_sq_b = (prefix_sum_s[bp] - prefix_sum_s[ap]) * (prefix_sum_s[bp] - prefix_sum_s[ap]) /
                        (prefix_sum_c[ap] - prefix_sum_c[bp]);
      batch_stream << identifer_batch[bix] << "," << ap << "," << bp - 1 << ',' << n << ","
                   << prefix_sum_s[bp] - prefix_sum_s[ap] << "," << chi_sq_b << "\n";
    } else {
      batch_stream << identifer_batch[bix] << "," << n << "," << n << ',' << n << "," << 0 << "," << 0 << "\n";
    }
  }
#pragma omp critical
  output_stream << batch_stream.rdbuf();
}

void SBatch::search_mers(const char* seq, uint64_t len, SSummary& or_summary, SSummary& rc_summary)
{
  uint32_t i, l;
  uint32_t orrix, rcrix;
  uint64_t orenc64_bp, orenc64_lr, rcenc64_bp;
  for (i = l = 0; i < len;) {
    if (SEQ_NT4_TABLE[seq[i]] >= 4) {
      l = 0, i++;
      continue;
    }
    l++, i++;
    if (l < k) {
      continue;
    }
    if (l == k) {
      compute_encoding(seq + i - k, seq + i, orenc64_lr, orenc64_bp);
    } else {
      update_encoding(seq + i - 1, orenc64_lr, orenc64_bp);
    }
    orenc64_bp = orenc64_bp & mask_bp;
    orenc64_lr = orenc64_lr & mask_lr;
    rcenc64_bp = revcomp_bp64(orenc64_bp, k);
    onmers++; // TODO: Incorporate missing partial/fraction?
    uint32_t hdist_curr = std::numeric_limits<uint32_t>::max();
    uint32_t hdist_or = hdist_curr, hdist_rc = hdist_curr;
#ifdef CANONICAL
    if (rcenc64_bp < orenc64_bp) {
      orrix = lshf->compute_hash(orenc64_bp);
      if (sketch->check_partial(orrix)) {
        hdist_or = or_summary.add_matching_mer(sketch, orrix, lshf->drop_ppos_lr(orenc64_lr));
      }
    } else {
      rcrix = lshf->compute_hash(rcenc64_bp);
      if (sketch->check_partial(rcrix)) {
        hdist_rc = rc_summary.add_matching_mer(sketch, rcrix, lshf->drop_ppos_lr(bp64_to_lr64(rcenc64_bp)));
      }
    }
#else
    orrix = lshf->compute_hash(orenc64_bp);
    if (sketch->check_partial(orrix)) {
      hdist_or = or_summary.add_matching_mer(sketch, orrix, lshf->drop_ppos_lr(orenc64_lr));
    }
    rcrix = lshf->compute_hash(rcenc64_bp);
    if (sketch->check_partial(rcrix)) {
      hdist_rc = rc_summary.add_matching_mer(sketch, rcrix, lshf->drop_ppos_lr(bp64_to_lr64(rcenc64_bp)));
    }
#endif /* CANONICAL */

    // if (i % 2 == 0) {
    //   continue;
    // }
    hdist_or = std::min(hdist_or, hdist_th + 1);
    hdist_rc = std::min(hdist_rc, hdist_th + 1);

    int j = i - k;
    c1[j] = sdd[hdist_or];
    c2[j] = sdd[hdist_rc];
    if (divergent) {
      s1[j] = -sd[hdist_or];
      s2[j] = -sd[hdist_rc];
      ssum1 += sd[hdist_or];
      ssum2 += sd[hdist_rc];
    } else {
      s1[j] = sd[hdist_or];
      s2[j] = sd[hdist_rc];
      ssum1 += sd[hdist_or];
      ssum2 += sd[hdist_rc];
    }
  }
}

uint32_t SSummary::add_matching_mer(sketch_sptr_t sketch, uint32_t rix, enc_t enc_lr)
{
  uint32_t hdist_min = sketch->search_mer(rix, enc_lr);
  if (hdist_min <= hdist_th) {
    mismatch_count--;
    match_count++;
    // hdisthist_v[hdist_min]++;
  }
  return hdist_min;
}
void SSummary::optimize_likelihood(optimize::HDistHistLLH llhfunc, double rho)
{
  llhfunc.set_parameters(hdisthist_v.data(), mismatch_count, rho);
  /* std::pair<double, double> sol_r = boost::math::tools::brent_find_minima(llhfunc, 1e-10, 0.5, 16); */
  /* d_llh = sol_r.first; */
  /* v_llh = sol_r.second; */
}
