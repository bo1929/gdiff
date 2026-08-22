#include "dist2.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <sys/mman.h>
#include <unistd.h>

#include "common.hpp"
#include "dist.hpp"
#include "enc.hpp"
#include "msg.hpp"
#include "tpool.hpp"

extern uint32_t num_threads;

namespace {

  using steady_clk = std::chrono::steady_clock;

  inline double ms_since(steady_clk::time_point t0)
  {
    return std::chrono::duration<double, std::milli>(steady_clk::now() - t0).count();
  }

  inline bool profile_enabled()
  {
    static const bool on = [] {
      const char* e = std::getenv("GDIFF2_PROFILE");
      return e != nullptr && e[0] != '\0' && e[0] != '0';
    }();
    return on;
  }

  inline uint32_t bucket_hdist_min_fast(const enc_t* ix1, const enc_t* ix2, const enc_t enc) noexcept
  {
    uint32_t hdist_min = std::numeric_limits<uint32_t>::max();
    if (ix1 < ix2) __builtin_prefetch(ix1, 0, 0);

    const enc_t* p = ix1;
    for (; p + 4 <= ix2; p += 4) {
      if (p + 8 <= ix2) __builtin_prefetch(p + 8, 0, 0);
      const uint32_t h0 = popcount_lr32(p[0] ^ enc);
      const uint32_t h1 = popcount_lr32(p[1] ^ enc);
      const uint32_t h2 = popcount_lr32(p[2] ^ enc);
      const uint32_t h3 = popcount_lr32(p[3] ^ enc);
      hdist_min = std::min(hdist_min, std::min(std::min(h0, h1), std::min(h2, h3)));
    }
    for (; p < ix2; ++p) {
      const uint32_t h = popcount_lr32((*p) ^ enc);
      hdist_min = h < hdist_min ? h : hdist_min;
    }
    return hdist_min;
  }

  struct pool_stats_t
  {
    uint64_t n_hash = 0;
    uint64_t n_bix_runs = 0;
    uint64_t n_empty_runs = 0;
    uint64_t n_empty_hashes = 0; // hashes skipped as miss without bucket scan
    uint64_t n_unobserved = 0;
    uint64_t n_enc_cmp = 0;
    uint64_t n_hit = 0;
    uint64_t n_miss = 0;
    double ms_total = 0;
  };

  // Accumulate a query's hash pool against a reference's LSH buckets into
  // per-window hit histograms. Kept exactly as-is from the single-pair dist2:
  // the u = nvalid - matched identity is preserved (no per-window miss counter
  // stored); empty-bucket/mismatched k-mers all land in u downstream.
  pool_stats_t accumulate_pool(const Sketch2& ref,
                               const win_hash_pool_t& pool,
                               uint32_t hdist_th,
                               uint64_t* hist_flat,
                               uint32_t nwins,
                               bool do_profile) noexcept
  {
    pool_stats_t st;
    const SFHM* sfhm = ref.get_sfhm_sptr().get();
    const uint32_t nrows = ref.get_nrows();
    const uint64_t* hashes = pool.hashes_ptr();
    const uint32_t* win_ix = pool.win_ix_ptr();
    const uint64_t n = pool.size();
    if (n == 0) return st;

    const auto t_all0 = do_profile ? steady_clk::now() : steady_clk::time_point{};

    // LSH-sorted pool in bix runs: one bucket resolve per distinct bix. The
    // pool only holds in-range (frac-observed) k-mers; frac-unobserved ones are
    // inferred per window as u = nvalid - t downstream.
    for (uint64_t i = 0; i < n;) {
      const uint32_t bix = static_cast<uint32_t>(hashes[i] >> 32);
      uint64_t j = i + 1;
      while (j < n && (hashes[j] >> 32) == bix)
        ++j;
      if (do_profile) {
        ++st.n_bix_runs;
        st.n_hash += (j - i);
        st.n_miss += (j - i);
      }

      if (bix >= nrows || !ref.bucket_nonempty(bix)) {
        if (do_profile) {
          ++st.n_empty_runs;
          st.n_empty_hashes += (j - i);
        }
        i = j;
        continue;
      }

      const enc_t* beg = sfhm->bucket_ptr_start(bix);
      const enc_t* end = sfhm->bucket_ptr_next(bix);
      const uint64_t blen = static_cast<uint64_t>(end - beg);
      if (beg < end) __builtin_prefetch(beg, 0, 0);

      for (uint64_t k = i; k < j; ++k) {
        const uint32_t wix = win_ix[k];
        if (wix >= nwins) continue;
        const uint32_t hd = bucket_hdist_min_fast(beg, end, static_cast<enc_t>(hashes[k] & 0xffffffffu));
        if (do_profile) st.n_enc_cmp += blen;
        if (hd <= hdist_th) {
          ++hist_flat[static_cast<size_t>(wix) * (hdist_bound + 1) + hd];
          if (do_profile) {
            ++st.n_hit;
            --st.n_miss;
          }
        }
      }
      i = j;
    }

    if (do_profile) st.ms_total = ms_since(t_all0);
    return st;
  }

  void print_pool_stats(const char* tag, const pool_stats_t& st)
  {
    const double avg_run = st.n_bix_runs ? static_cast<double>(st.n_hash) / static_cast<double>(st.n_bix_runs) : 0.0;
    const double empty_frac =
      st.n_bix_runs ? static_cast<double>(st.n_empty_runs) / static_cast<double>(st.n_bix_runs) : 0.0;
    const double avg_blen = (st.n_hash > st.n_empty_hashes)
                              ? static_cast<double>(st.n_enc_cmp) / static_cast<double>(st.n_hash - st.n_empty_hashes)
                              : 0.0;
    cerr_msg("[profile ",
             tag,
             "] hashes=",
             st.n_hash,
             " runs=",
             st.n_bix_runs,
             " avg_hashes/run=",
             avg_run,
             " empty_runs=",
             st.n_empty_runs,
             " empty_run_frac=",
             empty_frac,
             " empty_hashes_skipped=",
             st.n_empty_hashes,
             " n_unobserved=",
             st.n_unobserved,
             " enc_cmp=",
             st.n_enc_cmp,
             " avg_bucket_len=",
             avg_blen,
             " hits=",
             st.n_hit,
             " misses=",
             st.n_miss,
             " accumulate_ms=",
             st.ms_total);
  }

  // "\r"-based progress for the pair jobs, one bar redraw per reference batch
  // (each batch is a barrier: the coordinator prints after it drains). Pairs are
  // complete only when both directions ran, hence the /2.
  void print_pair_progress(uint64_t done_jobs, uint64_t total_jobs)
  {
    if (!isatty(STDERR_FILENO)) return; // keep redirected logs clean
    const int kBarWidth = 40;
    const double frac = total_jobs ? static_cast<double>(done_jobs) / static_cast<double>(total_jobs) : 1.0;
    const int filled = static_cast<int>(frac * kBarWidth + 0.5);
    std::ostringstream os;
    os << "\rpairs [" << std::string(static_cast<size_t>(filled), '#')
       << std::string(static_cast<size_t>(kBarWidth - filled), '.') << "] " << std::setw(3)
       << static_cast<int>(frac * 100.0 + 0.5) << "% (" << std::fixed << std::setprecision(1)
       << 0.5 * static_cast<double>(done_jobs) << "/" << 0.5 * static_cast<double>(total_jobs) << ")" << std::flush;
    std::cerr << os.str();
  }

  // True when `path` begins with the GSK5 magic (a real bundle); otherwise it is
  // assumed to be a list file of bundle paths.
  bool is_bundle_file(const std::filesystem::path& path)
  {
    std::ifstream in(path, std::ios::binary);
    if (!in) return false;
    uint32_t magic = 0;
    in.read(reinterpret_cast<char*>(&magic), sizeof(uint32_t));
    return in.good() && magic == SKETCH2_MAGIC;
  }

} // namespace

Dist2SC::dir_out_t Dist2SC::run_direction(const Sketch2& query,
                                          const Sketch2& reference,
                                          const str& dir_label,
                                          uint32_t hdist_th,
                                          bool output_samples)
{
  if (!query.compatible_with(reference)) {
    error_exit(concat_msg("Incompatible sketch2 pair: ",
                          query.get_rname(),
                          " vs ",
                          reference.get_rname(),
                          " (k/h/w/nrows/LSH/window params must match; use the same --seed)"));
  }
  if (!reference.has_buckets()) {
    error_exit(concat_msg("Reference sketch2 has no buckets: ",
                          reference.get_rname(),
                          " (dist2 needs full bundles in the reference role, not --windows-only)"));
  }
  if (!query.has_windows()) {
    error_exit(concat_msg("Query sketch2 has no windows loaded: ", query.get_rname()));
  }

  const bool do_profile = profile_enabled();
  const auto t_dir0 = do_profile ? steady_clk::now() : steady_clk::time_point{};

  LLH<double> llhf(reference.get_k(), reference.get_h(), reference.get_rho(), hdist_th, 0.0, false);
  const uint32_t k = query.get_k();
  const bool canonical = query.is_canonical();
  const uint32_t nwins = static_cast<uint32_t>(query.get_wins().size());

  const vec<win2_t>& wins_v = query.get_wins();
  vec<uint64_t> hist_fw(static_cast<size_t>(nwins) * (hdist_bound + 1), 0);
  vec<uint64_t> u_fw(nwins, 0);
  for (uint32_t wi = 0; wi < nwins; ++wi)
    u_fw[wi] = wins_v[wi].nvalid_fw;
  const auto t_acc0 = do_profile ? steady_clk::now() : steady_clk::time_point{};
  const pool_stats_t st_fw = accumulate_pool(reference, query.get_pool_fw(), hdist_th, hist_fw.data(), nwins, do_profile);
  double ms_acc = do_profile ? ms_since(t_acc0) : 0.0;

  vec<uint64_t> hist_rc;
  vec<uint64_t> u_rc;
  pool_stats_t st_rc;
  if (!canonical) {
    hist_rc.assign(static_cast<size_t>(nwins) * (hdist_bound + 1), 0);
    u_rc.assign(nwins, 0);
    for (uint32_t wi = 0; wi < nwins; ++wi)
      u_rc[wi] = static_cast<uint64_t>(wins_v[wi].nvalid_rc);
    const auto t_rc0 = do_profile ? steady_clk::now() : steady_clk::time_point{};
    st_rc = accumulate_pool(reference, query.get_pool_rc(), hdist_th, hist_rc.data(), nwins, do_profile);
    if (do_profile) ms_acc += ms_since(t_rc0);
  }

  vec<double> d_v;
  d_v.reserve(nwins);

  struct row_t
  {
    uint32_t wix;
    double d;
    char strand;
    const uint64_t* hist; // selected-strand histogram (fw or rc)
    uint64_t u;           // selected-strand miss count
  };
  vec<row_t> rows_v;
  if (output_samples) rows_v.reserve(nwins);

  for (uint32_t wix = 0; wix < nwins; ++wix) {
    const uint64_t* hfw = hist_fw.data() + static_cast<size_t>(wix) * (hdist_bound + 1);
    uint64_t t_fw = 0;
    for (uint32_t d = 0; d <= hdist_th; ++d)
      t_fw += hfw[d];
    // u = nvalid - matched: frac-unobserved/empty-bucket/mismatched k-mers all
    // count as misses, exactly matching scan.hpp's hdist_th+1 aggregation.
    u_fw[wix] = u_fw[wix] >= t_fw ? u_fw[wix] - t_fw : 0;
    const double d_fw = (t_fw == 0) ? nanx() : llhf.mle(hfw, u_fw[wix]);

    double d = d_fw;
    char strand = canonical ? '.' : '+';
    const uint64_t* hist_sel = hfw;
    uint64_t u_sel = u_fw[wix];
    if (!canonical) {
      const uint64_t* hrc = hist_rc.data() + static_cast<size_t>(wix) * (hdist_bound + 1);
      uint64_t t_rc = 0;
      for (uint32_t di = 0; di <= hdist_th; ++di)
        t_rc += hrc[di];
      u_rc[wix] = u_rc[wix] >= t_rc ? u_rc[wix] - t_rc : 0;
      const double d_rc = (t_rc == 0) ? nanx() : llhf.mle(hrc, u_rc[wix]);
      const auto picked = select_strand_distance(d_fw, d_rc);
      d = picked.first;
      strand = picked.second;
      if (strand == '-') {
        hist_sel = hrc;
        u_sel = u_rc[wix];
      }
    }

    if (is_valid_distance(d)) d_v.push_back(d);

    if (output_samples) rows_v.push_back({wix, d, strand, hist_sel, u_sel});
  }

  std::sort(d_v.begin(), d_v.end());
  const double d_median = linear_quantile(d_v, 0.5);
  const double ms_mle = do_profile ? ms_since(t_dir0) : 0.0;

  dir_out_t out;
  out.d_median = d_median;
  out.n_valid = static_cast<uint64_t>(d_v.size());
  out.n_unmapped = static_cast<uint64_t>(nwins) - out.n_valid;

  if (do_profile) {
    print_pool_stats("fw", st_fw);
    if (!canonical) print_pool_stats("rc", st_rc);
    cerr_msg("[profile dir=",
             dir_label,
             "] acc_ms=",
             ms_acc,
             " mle_ms=",
             ms_mle,
             " total_ms=",
             ms_since(t_dir0),
             " nwins=",
             nwins,
             " ref_nkmers=",
             reference.get_sfhm_sptr()->get_nkmers(),
             " nrows=",
             reference.get_nrows());
  }

  if (output_samples) {
    strstream ss;
    set_precision(ss, 5);
    for (const row_t& row : rows_v) {
      double lr_bg = nanx();
      double lr_ub = nanx();
      if (is_valid_distance(row.d)) {
        uint64_t n_total = row.u;
        for (uint32_t di = 0; di <= hdist_th; ++di)
          n_total += row.hist[di];
        lr_ub = compute_lr_ub(llhf, row.d, n_total);
        if (is_valid_distance(d_median))
          lr_bg = likelihood_ratio_statistic(llhf.nll(d_median, row.hist, row.u), llhf.nll(row.d, row.hist, row.u));
      }
      const win2_t& win = query.get_wins()[row.wix];
      write_tsv(
        ss, dir_label, win.qid, win.start + 1, win.end + k - 1, row.strand, reference.get_rname(), row.d, lr_bg, lr_ub)
        << "\n";
    }
    out.samples = ss.str();
  }

  return out;
}

void Dist2SC::resolve_source(const std::filesystem::path& path, std::vector<uint32_t>& out)
{
  if (is_bundle_file(path)) {
    const bundle_t* bundle = bundle_for(path);
    for (uint32_t rec = 0; rec < bundle->idx.records.size(); ++rec) {
      bool found = false;
      for (uint32_t mi = 0; mi < members.size(); ++mi) {
        if (members[mi].bundle == bundle && members[mi].rec == rec) {
          out.push_back(mi);
          found = true;
          break;
        }
      }
      if (!found) {
        members.push_back({bundle, rec, str{}, nullptr, nullptr});
        out.push_back(static_cast<uint32_t>(members.size() - 1));
      }
    }
    return;
  }

  // Not a bundle: `path` is a list file of bundle paths (one per line), with
  // relative entries resolved against the list file's directory.
  std::ifstream in(path);
  check_fstream(in, "Cannot open input list", path.string());
  const std::filesystem::path base_dir = path.parent_path();
  str line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    line.erase(line.find_last_not_of(" \t\r\n") + 1);
    if (line.empty()) continue;
    std::filesystem::path entry(line);
    if (entry.is_relative()) entry = base_dir / entry;
    resolve_source(entry, out);
  }
}

const Dist2SC::bundle_t* Dist2SC::bundle_for(const std::filesystem::path& path)
{
  for (const auto& b : bundles) {
    if (b->path == path) return b.get();
  }
  auto b = std::make_unique<bundle_t>();
  b->path = path;
  b->idx = read_sketch2_index(path);
  bundles.push_back(std::move(b));
  return bundles.back().get();
}

void Dist2SC::dist()
{
  const bool do_profile = profile_enabled();
  const auto t_all0 = std::chrono::steady_clock::now();

  const bool have_pos_q = !set_a_path.empty();
  const bool have_pos_r = !set_b_path.empty();
  const bool have_list_q = !query_list_path.empty();
  const bool have_list_r = !reference_list_path.empty();

  if (have_pos_q && have_list_q) error_exit("Query set given both positionally and via --query-list; pick one");
  if (have_pos_r && have_list_r) error_exit("Reference set given both positionally and via --reference-list; pick one");

  if (have_list_q) {
    resolve_source(query_list_path, set_a);
  } else if (have_pos_q) {
    resolve_source(set_a_path, set_a);
  } else {
    error_exit("No set given: pass `dist2 <set>` (all-pairs within) or `dist2 <q> <r>` (cross)");
  }

  if (have_list_r)
    resolve_source(reference_list_path, set_b);
  else if (have_pos_r)
    resolve_source(set_b_path, set_b);

  if (set_a.empty()) error_exit("Query set is empty; nothing to compare");

  // Deduplicate (same bundle record listed in several sets/file lines).
  auto dedupe = [](std::vector<uint32_t>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  };
  dedupe(set_a);
  dedupe(set_b);

  const bool within = set_b.empty();
  cerr_msg(
    "members=", members.size(), " query=", set_a.size(), within ? " within-mode" : " cross-mode", " ref=", set_b.size());

  // --- Load all member windows up front (query role). Reference buckets are
  //     loaded per-batch in the comparison loop below.
  const auto t_load0 = do_profile ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
  ThreadPool pool(std::max(1u, num_threads));
  pool.parallel_for(members.size(), 1, [&](uint64_t mi) {
    member_t& m = members[static_cast<size_t>(mi)];
    m.windows = std::make_unique<Sketch2>(m.bundle->path);
    m.windows->load_from_offset(m.bundle->idx.records[m.rec], m.bundle->idx, Sketch2Part::Windows);
  });
  if (do_profile) cerr_msg("[profile] windows_load_ms=", ms_since(t_load0));

  // --- Build the ordered pair list -------------------------------------------
  vec<pair_t> pairs;
  if (within) {
    pairs.reserve(set_a.size() * (set_a.size() - 1) / 2);
    for (size_t i = 0; i < set_a.size(); ++i)
      for (size_t j = i + 1; j < set_a.size(); ++j)
        pairs.push_back({set_a[i], set_a[j]});
  } else {
    pairs.reserve(set_a.size() * set_b.size());
    for (uint32_t q : set_a)
      for (uint32_t r : set_b)
        pairs.push_back({q, r});
  }
  const uint64_t npairs = pairs.size();
  if (npairs == 0) {
    cerr_msg("No pairs to compare (set size < 2 in within mode)");
    return;
  }

  // --- Directional jobs; reference-major so each reference's buckets are
  //     loaded once and evicted before the next reference is loaded.
  struct job_t
  {
    uint32_t pair_ix;
    uint32_t query_ix; // member providing windows
    uint32_t ref_ix;   // member providing buckets
    bool ba;           // false → ab (query = pair.a), true → ba (query = pair.b)
  };
  std::vector<job_t> jobs;
  jobs.reserve(npairs * 2);
  for (uint32_t p = 0; p < npairs; ++p) {
    jobs.push_back({p, pairs[p].a, pairs[p].b, false});
    jobs.push_back({p, pairs[p].b, pairs[p].a, true});
  }
  std::stable_sort(jobs.begin(), jobs.end(), [](const job_t& x, const job_t& y) { return x.ref_ix < y.ref_ix; });

  // --- Reference-major batches ------------------------------------------------
  const auto t_cmp0 = do_profile ? std::chrono::steady_clock::now() : std::chrono::steady_clock::time_point{};
  const uint64_t total_jobs = static_cast<uint64_t>(jobs.size());
  uint64_t done_jobs = 0;
  size_t batch_start = 0;
  while (batch_start < jobs.size()) {
    const uint32_t ref_ix = jobs[batch_start].ref_ix;
    size_t batch_end = batch_start;
    while (batch_end < jobs.size() && jobs[batch_end].ref_ix == ref_ix)
      ++batch_end;

    member_t& ref = members[ref_ix];
    ref.buckets = std::make_unique<Sketch2>(ref.bundle->path);
    ref.buckets->load_from_offset(ref.bundle->idx.records[ref.rec], ref.bundle->idx, Sketch2Part::Buckets);
    const Sketch2* refp = ref.buckets.get();

    pool.parallel_for(batch_end - batch_start, 1, [&](uint64_t i) {
      const job_t& job = jobs[batch_start + static_cast<size_t>(i)];
      const member_t& q = members[job.query_ix];
      const str dir = job.ba ? "ba" : "ab";
      dir_out_t r = run_direction(*q.windows, *refp, dir, hdist_th, output_samples);
      pair_t& pr = pairs[job.pair_ix];
      if (job.ba)
        pr.ba = std::move(r);
      else
        pr.ab = std::move(r);
    });

    // Evict the reference's bucket pages now that the batch is done.
    const sketch2_entry_t& e = ref.bundle->idx.records[ref.rec];
    refp->advise_range(e.buckets_off, e.buckets_len, MADV_DONTNEED);
    ref.buckets.reset();

    done_jobs += static_cast<uint64_t>(batch_end - batch_start);
    print_pair_progress(done_jobs, total_jobs);

    batch_start = batch_end;
  }
  if (isatty(STDERR_FILENO)) std::cerr << '\n';
  if (do_profile) cerr_msg("[profile] compare_ms=", ms_since(t_cmp0));

  // --- Emit in deterministic pair order ----------------------------------------
  std::ostream& os = *output_stream;
  if (output_samples) {
    for (const pair_t& pr : pairs) {
      os << pr.ab.samples;
      os << pr.ba.samples;
    }
  } else {
    set_precision(os, 8);
    for (const pair_t& pr : pairs) {
      const double d_avg =
        (std::isfinite(pr.ab.d_median) && std::isfinite(pr.ba.d_median)) ? (pr.ab.d_median + pr.ba.d_median) / 2.0 : nanx();
      write_tsv(os,
                members[pr.a].windows->get_rname(),
                members[pr.b].windows->get_rname(),
                pr.ab.n_valid,
                pr.ab.d_median,
                pr.ba.n_valid,
                pr.ba.d_median,
                d_avg)
        << '\n';
    }
  }
  os.flush();

  if (do_profile)
    cerr_msg("[profile total] load=",
             ms_since(t_load0),
             " cmp=",
             ms_since(t_cmp0),
             " pairs=",
             npairs,
             " total_ms=",
             ms_since(t_all0));
}

Dist2SC::Dist2SC(CLI::App& sc)
{
  sc.add_option(
      "set", set_a_path, "Sketch2 bundle or list file (the single set for within-pairs; also the query set in cross mode)")
    ->expected(0, 1)
    ->check(CLI::ExistingFile);
  sc.add_option("refset", set_b_path, "Sketch2 bundle or list file (reference set for cross mode)")
    ->expected(0, 1)
    ->check(CLI::ExistingFile);
  sc.add_option("--query-list", query_list_path, "Query set as a list file (alternative to the first positional)")
    ->excludes("set")
    ->check(CLI::ExistingFile);
  sc.add_option(
      "--reference-list", reference_list_path, "Reference set as a list file (alternative to the second positional)")
    ->excludes("refset")
    ->check(CLI::ExistingFile);
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match [4]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.add_flag("--output-samples", output_samples, "Write per-window sample output instead of per-pair summary");
  sc.callback([&]() {
    if (!validate_configuration()) error_exit("Invalid configuration!");
    if (!output_path.empty()) {
      output_file.open(output_path);
      check_fstream(output_file, "Cannot open output file", output_path.string());
      output_stream = &output_file;
    }
  });
}

bool Dist2SC::validate_configuration()
{
  bool is_invalid = false;
  if (hdist_th > hdist_bound) {
    is_invalid = true;
    cerr_msg("--hdist-th must be in [0, ", hdist_bound, "]; got ", hdist_th);
  }
  return !is_invalid;
}