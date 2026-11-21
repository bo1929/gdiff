#include "dist.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <memory>
#include <numeric>
#include <sstream>
#include <sys/mman.h>
#include <utility>
#include <unistd.h>

#include "utils.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "records.hpp"
#include "scan.hpp"
#include "tpool.hpp"
#include "windows.hpp"

uint32_t hdist_th_for(uint64_t nvalid_bp, std::optional<uint32_t> hdist_th)
{
  if (hdist_th) return *hdist_th;
  return nvalid_bp >= (20ull * 1000 * 1000) ? 2u : 3u;
}

namespace {

  // One strand's evidence for every window of a profile: match histograms and misses, flattened.
  struct evidence_t
  {
    vec<uint64_t> hist_v; // nwins * (hdist_bound + 1)
    vec<uint64_t> u_v;    // nwins

    void assign(size_t nwins)
    {
      hist_v.assign(nwins * (hdist_bound + 1), 0);
      u_v.assign(nwins, 0);
    }
    const uint64_t* window_hist(size_t wix) const noexcept { return hist_v.data() + wix * (hdist_bound + 1); }
    uint64_t* window_hist(size_t wix) noexcept { return hist_v.data() + wix * (hdist_bound + 1); }

    // The valid k-mers the histogram did not match are the misses.
    void set_u(size_t wix, uint64_t nvalid, uint32_t hdist_th) noexcept
    {
      const uint64_t matched = hist_total(window_hist(wix), hdist_th);
      u_v[wix] = nvalid > matched ? nvalid - matched : 0;
    }

    uint64_t get_nmers(size_t wix, uint32_t hdist_th) const noexcept
    {
      return u_v[wix] + hist_total(window_hist(wix), hdist_th);
    }

    // NaN when no k-mer matched within the threshold.
    double mle(const LLH<double>& llhf, size_t wix, uint32_t hdist_th) const
    {
      const uint64_t* hist = window_hist(wix);
      return hist_total(hist, hdist_th) == 0 ? nanx() : llhf.mle(hist, u_v[wix]);
    }
  };

  // The pool is grouped by bucket, so each bucket is resolved once per sample.
  void accumulate_pool(const Buckets& buckets,
                       const hash_pool_t& pool,
                       uint32_t hdist_th,
                       size_t nwins,
                       evidence_t& out) noexcept
  {
    const uint64_t* hashes = pool.hashes_ptr();
    const uint16_t* win_ix = pool.win_ix_ptr();
    const uint64_t n = pool.size();
    uint64_t j = 0;
    for (uint64_t i = 0; i < n; i = j) {
      const uint32_t bix = key_bix(hashes[i]);
      for (j = i + 1; j < n && key_bix(hashes[j]) == bix;)
        ++j;
      const enc_t* beg = nullptr;
      const enc_t* end = nullptr;
      if (!buckets.range(bix, beg, end)) continue;
      __builtin_prefetch(beg, 0, 0);
      for (uint64_t e = i; e < j; ++e) {
        const size_t wix = win_ix[e];
        if (wix >= nwins) continue;
        const uint32_t hd = bucket_hdist_min(beg, end, key_enc(hashes[e]));
        if (hd <= hdist_th) ++out.window_hist(wix)[hd];
      }
    }
  }

  // Both strands of one profile; rc stays empty for canonical sketches.
  struct strands_t
  {
    evidence_t fw;
    evidence_t rc;
    bool has_rc = false;
  };

  // Counts the k-mers of one window straight into its slot of each strand.
  struct window_agg_t
  {
    strands_t& strands;
    size_t wix;
    uint32_t hdist_th;

    void operator()(uint64_t, uint32_t hd, bool is_rc) const noexcept
    {
      evidence_t& ev = is_rc ? strands.rc : strands.fw;
      if (hd <= hdist_th) {
        ++ev.window_hist(wix)[hd];
      } else {
        ++ev.u_v[wix];
      }
    }
    void skip_mer(uint64_t) const noexcept {}
  };

  // Per-window counts of one profile.
  strands_t query_windows(const Sketch& query, const Sketch& reference, uint32_t hdist_th)
  {
    const vec<window_t>& wins_v = query.get_windows().wins_v;
    const size_t nwins = wins_v.size();
    const bool canonical = query.is_canonical();
    strands_t strands;
    strands.has_rc = !canonical;
    strands.fw.assign(nwins);
    if (strands.has_rc) strands.rc.assign(nwins);

    if (!query.get_config().keep_seq) {
      const Buckets& buckets = reference.get_buckets();
      accumulate_pool(buckets, query.get_windows().pool_fw, hdist_th, nwins, strands.fw);
      if (strands.has_rc) accumulate_pool(buckets, query.get_windows().pool_rc, hdist_th, nwins, strands.rc);
      // A pool holds only retained k-mers; the rest are recovered as misses.
      for (size_t wix = 0; wix < nwins; ++wix) {
        strands.fw.set_u(wix, wins_v[wix].nvalid_fw, hdist_th);
        if (strands.has_rc) strands.rc.set_u(wix, wins_v[wix].nvalid_rc, hdist_th);
      }
      return strands;
    }

    // The scan counts the misses itself.
    const scan_ctx_t ctx = make_scan_ctx(reference, 0, hdist_th);
    const seq_pack_t& packs = query.get_windows().packs;
    str cseq;
    for (size_t wix = 0; wix < nwins; ++wix) {
      packs.unpack(wix, cseq);
      // Compatible sketches share a window length, so the query's is authoritative.
      const uint64_t nmers = std::min(query.get_config().tau, wins_v[wix].end - wins_v[wix].start);
      const window_agg_t agg{strands, wix, hdist_th};
      if (canonical) {
        scan_mers_range<true>(ctx, cseq.data(), 0, nmers, agg);
      } else {
        scan_mers_range<false>(ctx, cseq.data(), 0, nmers, agg);
      }
    }
    return strands;
  }

  void require_comparable(const Sketch& query, const Sketch& reference)
  {
    if (!compatible_configs(query.get_config(), reference.get_config())) {
      error_exit(concat_msg("Incompatible sketch pair: ",
                            query.get_rname(),
                            " vs ",
                            reference.get_rname(),
                            " (k/w/h/LSH/window length must match; use the same --seed)"));
    }
    if (!reference.has_buckets()) error_exit(concat_msg("Reference sketch has no buckets: ", reference.get_rname()));
    if (!query.has_windows()) error_exit(concat_msg("Query sketch has no sampled windows: ", query.get_rname()));
  }

  double mean_of(const vec<double>& v)
  {
    if (v.empty()) return nanx();
    return std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
  }

  double max_of(const vec<double>& v) { return v.empty() ? nanx() : *std::max_element(v.begin(), v.end()); }

  // A window's reported estimate; strand-aware sketches also record the strand that supplied it.
  mle_t estimate_window(const LLH<double>& llhf, const strands_t& strands, size_t wix, uint32_t hdist_th)
  {
    mle_t m;
    m.d = strands.fw.mle(llhf, wix, hdist_th);
    const evidence_t* ev = &strands.fw;
    if (strands.has_rc) {
      const auto picked = select_strand_distance(m.d, strands.rc.mle(llhf, wix, hdist_th));
      m.d = picked.first;
      m.strand = picked.second;
      if (m.strand == '-') ev = &strands.rc;
    }
    if (is_valid_distance(m.d)) m.lr_ub = compute_lr_ub(llhf, m.d, ev->get_nmers(wix, hdist_th));
    return m;
  }

} // namespace

profile_t obtain_profile(const Sketch& query, const Sketch& reference, uint32_t hdist_th)
{
  require_comparable(query, reference);

  const LLH<double> llhf = make_llhf(reference, hdist_th);
  const strands_t strands = query_windows(query, reference, hdist_th);
  const size_t nwins = query.get_windows().wins_v.size();

  profile_t profile;
  profile.windows_v.reserve(nwins);
  vec<double> d_v;
  d_v.reserve(nwins);
  for (size_t wix = 0; wix < nwins; ++wix) {
    const mle_t m = estimate_window(llhf, strands, wix, hdist_th);
    if (is_valid_distance(m.d)) d_v.push_back(m.d);
    profile.windows_v.push_back(m);
  }

  std::sort(d_v.begin(), d_v.end());
  profile.d_median = linear_quantile(d_v, 0.5);
  profile.d_mean = mean_of(d_v);
  profile.d_highest = max_of(d_v);
  profile.n_valid = d_v.size();
  return profile;
}

namespace {

  bool order_finite_first(const mle_t& lhs, const mle_t& rhs)
  {
    const bool lnan = std::isnan(lhs.d);
    const bool rnan = std::isnan(rhs.d);
    if (lnan || rnan) return !lnan && rnan;
    return lhs.d < rhs.d;
  }

  // Sorts each direction and keeps, at every rank, the smaller distance (or the only finite one).
  // Ranks with no finite distance in either direction are counted as NA.
  std::pair<vec<mle_t>, uint64_t> reconcile_by_rank(vec<mle_t>& ab_v, vec<mle_t>& ba_v)
  {
    std::sort(ab_v.begin(), ab_v.end(), order_finite_first);
    std::sort(ba_v.begin(), ba_v.end(), order_finite_first);

    uint64_t n_na = 0;
    const size_t nranks = std::max(ab_v.size(), ba_v.size());
    vec<mle_t> mle_v;
    mle_v.reserve(nranks);
    for (size_t i = 0; i < nranks; ++i) {
      const bool has_ab = i < ab_v.size() && !std::isnan(ab_v[i].d);
      const bool has_ba = i < ba_v.size() && !std::isnan(ba_v[i].d);
      if (has_ab && has_ba) {
        mle_v.push_back(ab_v[i].d <= ba_v[i].d ? ab_v[i] : ba_v[i]);
      } else if (has_ab) {
        mle_v.push_back(ab_v[i]);
      } else if (has_ba) {
        mle_v.push_back(ba_v[i]);
      } else {
        ++n_na;
      }
    }
    return {std::move(mle_v), n_na};
  }

} // namespace

summary_t summarize_symmetric(vec<mle_t> ab_v, vec<mle_t> ba_v, double lr_th, double min_portion)
{
  summary_t s;
  const auto [mle_v, n_na] = reconcile_by_rank(ab_v, ba_v);
  s.n_na = n_na;

  // A missing likelihood-ratio statistic carries no evidence, so its window is rejected.
  vec<double> all_d_v, kept_d_v;
  all_d_v.reserve(mle_v.size());
  kept_d_v.reserve(mle_v.size());
  for (const mle_t& m : mle_v) {
    if (m.lr_ub == 0.0) ++s.n_ub;
    all_d_v.push_back(m.d);
    if (m.lr_ub > lr_th) kept_d_v.push_back(m.d);
  }
  s.n_filtered = all_d_v.size() - kept_d_v.size();

  const double portion = all_d_v.empty() ? 0.0 : static_cast<double>(kept_d_v.size()) / static_cast<double>(all_d_v.size());
  const bool is_filtered = portion > min_portion;
  s.d = mean_of(is_filtered ? kept_d_v : all_d_v);
  s.d_v = std::move(is_filtered ? kept_d_v : all_d_v);
  return s;
}

void DistSC::write_header(std::ostream& os, bool canonical) const
{
  tsv_row row(os);
  if (params.output_samples) {
    row("config")("genome_a")("genome_b")("seq")("start")("end");
    if (!canonical) row("strand");
    row("direction")("d")("p_ub");
  } else {
    row("genome_a")("genome_b")("d");
    row("d_median_ab")("d_median_ba")("d_mean_ab")("d_mean_ba")("d_highest_ab")("d_highest_ba");
    row("n_ab")("n_ba")("n_ub")("n_na")("n_filtered");
  }
  os << '\n';
}

void DistSC::write_pair_line(std::ostream& os, const pair_t& pr) const
{
  const summary_t s = summarize_symmetric(pr.ab.windows_v, pr.ba.windows_v, params.lr_th, params.min_portion);
  tsv_row row(os);
  row(entries_v[pr.a].rname)(entries_v[pr.b].rname)(s.d);
  row(pr.ab.d_median)(pr.ba.d_median)(pr.ab.d_mean)(pr.ba.d_mean)(pr.ab.d_highest)(pr.ba.d_highest);
  row(pr.ab.n_valid)(pr.ba.n_valid)(s.n_ub)(s.n_na)(s.n_filtered);
  os << '\n';
}

void DistSC::write_samples(std::ostream& os, const pair_t& pr, const vec<Sketch>& targets_v) const
{
  const str& name_a = entries_v[pr.a].rname;
  const str& name_b = entries_v[pr.b].rname;
  const struct
  {
    const Sketch& query;
    const profile_t& profile;
    const char* direction;
  } directions[] = {{targets_v[pr.a], pr.ab, "ab"}, {targets_v[pr.b], pr.ba, "ba"}};

  for (const auto& [query, profile, direction] : directions) {
    const vec<window_t>& wins_v = query.get_windows().wins_v;
    for (size_t wix = 0; wix < wins_v.size(); ++wix) {
      const window_t& win = wins_v[wix];
      const mle_t& m = profile.windows_v[wix];
      tsv_row row(os);
      row("gdiff")(name_a)(name_b)(win.qid)(win.start + 1)(win.end + query.get_k() - 1);
      if (!query.is_canonical()) row(m.strand);
      row(direction)(m.d)(lr_pvalue(m.lr_ub));
      os << '\n';
    }
  }
}

void DistSC::resolve_container(const std::filesystem::path& path, vec<uint32_t>& out)
{
  auto it = std::find_if(loaded_v.begin(), loaded_v.end(), [&](const loaded_t& l) { return l.file->get_path() == path; });
  if (it == loaded_v.end()) {
    loaded_v.push_back({std::make_unique<Container>(path), static_cast<uint32_t>(entries_v.size())});
    it = std::prev(loaded_v.end());
    for (uint32_t rix = 0; rix < it->file->size(); ++rix)
      entries_v.push_back({it->file.get(), rix, str{}});
  }
  for (uint32_t rix = 0; rix < it->file->size(); ++rix)
    out.push_back(it->first + rix);
}

void DistSC::resolve_list(const std::filesystem::path& path, vec<uint32_t>& out)
{
  std::ifstream in(path);
  check_fstream(in, "Cannot open input list", path.string());
  const std::filesystem::path base_dir = path.parent_path();
  str line;
  while (std::getline(in, line)) {
    str name;
    std::filesystem::path entry;
    if (!parse_input_entry(line, name, entry)) continue;
    if (entry.is_relative()) entry = base_dir / entry;
    if (!is_container_file(entry)) error_exit("List entry is not a sketch container: " + entry.string());
    resolve_container(entry, out);
  }
}

void DistSC::resolve_set(const std::filesystem::path& container_path,
                         const std::filesystem::path& list_path,
                         char set,
                         vec<uint32_t>& out)
{
  // The CLI already rejects a set given both ways.
  if (!list_path.empty()) {
    resolve_list(list_path, out);
  } else if (!container_path.empty()) {
    if (!is_container_file(container_path)) {
      error_exit(
        concat_msg("Not a sketch container: ", container_path.string(), "; pass a list of containers with --list-", set));
    }
    resolve_container(container_path, out);
  }
}

// Every container needs windows and must be comparable with the first one.
void DistSC::require_comparable_entries() const
{
  const Container& first = *entries_v.front().file;
  for (const entry_t& m : entries_v) {
    const Container& file = *m.file;
    if (file.get_config().tau == 0) {
      error_exit("Sketch has no sampled windows: " + file.get_path().string() +
                 "; re-run `gdiff sketch` with a window length (-l)");
    }
    if (!compatible_configs(first.get_config(), file.get_config())) {
      error_exit(concat_msg("Incompatible containers: ",
                            first.get_path().string(),
                            " vs ",
                            file.get_path().string(),
                            " (k/w/h/LSH/window length must match; use the same --seed)"));
    }
  }
}

vec<Sketch> DistSC::load_query_sketches(ThreadPool& pool)
{
  vec<Sketch> targets_v(entries_v.size());
  pool.parallel_for(entries_v.size(), 1, [&](uint64_t mi) {
    const entry_t& m = entries_v[static_cast<size_t>(mi)];
    targets_v[static_cast<size_t>(mi)] = m.file->open(m.rix, SketchLoad::Windows);
  });
  for (size_t mi = 0; mi < entries_v.size(); ++mi)
    entries_v[mi].rname = targets_v[mi].get_rname();
  return targets_v;
}

// Within mode compares every unordered pair of set A; cross mode the product of A and B.
vec<DistSC::pair_t> DistSC::plan_pairs() const
{
  const bool is_within = set_b_v.empty();
  const vec<uint32_t>& refs_v = is_within ? set_a_v : set_b_v;
  vec<pair_t> pairs_v;
  pairs_v.reserve(is_within ? set_a_v.size() * (set_a_v.size() - 1) / 2 : set_a_v.size() * set_b_v.size());
  for (size_t i = 0; i < set_a_v.size(); ++i) {
    for (size_t j = is_within ? i + 1 : 0; j < refs_v.size(); ++j)
      pairs_v.push_back({set_a_v[i], refs_v[j], {}, {}});
  }

  for (pair_t& pr : pairs_v) {
    if (entries_v[pr.b].rname < entries_v[pr.a].rname) std::swap(pr.a, pr.b);
  }
  std::stable_sort(pairs_v.begin(), pairs_v.end(), [&](const pair_t& x, const pair_t& y) {
    const str& xa = entries_v[x.a].rname;
    const str& ya = entries_v[y.a].rname;
    if (xa != ya) return xa < ya;
    return entries_v[x.b].rname < entries_v[y.b].rname;
  });
  return pairs_v;
}

// Scores both directions of every pair.
// Jobs sharing a source sketch run together, so each source's buckets are loaded once and released afterwards.
void DistSC::score_pairs(ThreadPool& pool, const vec<Sketch>& targets_v, vec<pair_t>& pairs_v)
{
  struct job_t
  {
    uint32_t pair_ix;
    uint32_t target_ix;
    uint32_t source_ix;
    bool is_ba;
  };
  vec<job_t> jobs_v;
  jobs_v.reserve(pairs_v.size() * 2);
  for (uint32_t p = 0; p < pairs_v.size(); ++p) {
    jobs_v.push_back({p, pairs_v[p].a, pairs_v[p].b, false});
    jobs_v.push_back({p, pairs_v[p].b, pairs_v[p].a, true});
  }
  std::stable_sort(jobs_v.begin(), jobs_v.end(), [](const job_t& x, const job_t& y) { return x.source_ix < y.source_ix; });

  const uint64_t total_jobs = jobs_v.size();
  uint64_t done_jobs = 0;
  size_t batch_start = 0;
  while (batch_start < jobs_v.size()) {
    const uint32_t source_ix = jobs_v[batch_start].source_ix;
    size_t batch_end = batch_start;
    while (batch_end < jobs_v.size() && jobs_v[batch_end].source_ix == source_ix)
      ++batch_end;

    const entry_t& sentry = entries_v[source_ix];
    const Sketch source = sentry.file->open(sentry.rix, SketchLoad::Buckets);

    const uint64_t nvalid_bp = source.get_nvalid_bp();
    const uint32_t hdist_th = hdist_th_for(nvalid_bp, params.hdist_th);

    pool.parallel_for(batch_end - batch_start, 1, [&](uint64_t i) {
      const job_t& job = jobs_v[batch_start + static_cast<size_t>(i)];
      profile_t profile = obtain_profile(targets_v[job.target_ix], source, hdist_th);
      pair_t& pr = pairs_v[job.pair_ix];
      (job.is_ba ? pr.ba : pr.ab) = std::move(profile);
    });

    const scentry& e = sentry.file->get_entry(sentry.rix);
    sentry.file->advise(e.buckets_offset, e.buckets_len, MADV_DONTNEED);

    done_jobs += batch_end - batch_start;
    progress_continue("Obtained profiles", done_jobs, total_jobs);
    batch_start = batch_end;
  }
  progress_done();
}

void DistSC::estimate_distances()
{
  if (set_a_v.empty()) {
    error_exit("Set A is empty; nothing to compare");
  }

  auto dedup = [](vec<uint32_t>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  };
  dedup(set_a_v);
  dedup(set_b_v);

  if (set_b_v.empty()) {
    info_msg("Comparing ", set_a_v.size(), " sketch(es) within set A");
  } else {
    info_msg("Comparing ", set_a_v.size(), " sketch(es) in set A against ", set_b_v.size(), " in set B");
  }
  require_comparable_entries();

  ThreadPool pool(std::max(1u, num_threads));
  const vec<Sketch> targets_v = load_query_sketches(pool);
  vec<pair_t> pairs_v = plan_pairs();

  std::ostream& os = *output_stream;
  write_provenance(os);
  write_header(os, entries_v.front().file->get_config().canonical);
  if (pairs_v.empty()) {
    warn_msg("No pairs to compare");
    os.flush();
    return;
  }

  score_pairs(pool, targets_v, pairs_v);

  if (params.output_samples) {
    set_precision(os, 5);
    for (const pair_t& pr : pairs_v)
      write_samples(os, pr, targets_v);
  } else {
    set_precision(os, 8);
    for (const pair_t& pr : pairs_v)
      write_pair_line(os, pr);
  }
  os.flush();
}

void DistSC::dist()
{
  if (container_a_path.empty() && list_a_path.empty()) {
    error_exit("No set A given: pass a sketch container, or --list-a");
  }
  resolve_set(container_a_path, list_a_path, 'a', set_a_v);
  resolve_set(container_b_path, list_b_path, 'b', set_b_v);
  estimate_distances();
}

DistSC::DistSC(CLI::App& sc)
{
  sc.description("Pairwise distance between two sets AxB or within AxA (when B is not given).");
  sc.add_option("sketch-a", container_a_path, "Sketch container for set A (omit when using --list-a)")
    ->check(CLI::ExistingFile);
  sc.add_option("sketch-b", container_b_path, "Sketch container for set B (omit when using --list-b)")
    ->check(CLI::ExistingFile);
  sc.add_option("--list-a", list_a_path, "Set A as a list file of sketch containers")
    ->excludes("sketch-a")
    ->check(CLI::ExistingFile);
  sc.add_option("--list-b", list_b_path, "Set B as a list file of sketch containers")
    ->excludes("sketch-b")
    ->check(CLI::ExistingFile);
  sc.add_option("--hdist-th", params.hdist_th, "Max. Hamming distance for k-mer search [3 (<20Mbp) or 2 (>20Mb)]");
  sc.add_option("--lr-th", params.lr_th, "Likelihood-ratio cut for the reconciliation filter [6.635]");
  sc.add_option("--min-portion", params.min_portion, "Filter only if at least --min-portion exceeds --lr-th [0.66]");
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.add_flag("--output-samples", params.output_samples, "Write per-window sample lines instead of per-pair summaries");
  sc.callback([&]() {
    if (!validate_configuration()) exit_failure();
    open_output(output_file, output_path, output_stream);
  });
}

bool DistSC::validate_configuration()
{
  bool is_valid = true;
  auto require = [&](bool condition, const auto&... msg) {
    if (condition) return;
    is_valid = false;
    error_msg(msg...);
  };

  const bool is_hdist_valid = !params.hdist_th || *params.hdist_th <= hdist_bound;
  require(is_hdist_valid, "--hdist-th must be in [0, ", hdist_bound, "]; got ", params.hdist_th.value_or(0));
  require(params.lr_th >= 0.0, "--lr-th must be non-negative; got ", params.lr_th);
  require(
    params.min_portion >= 0.0 && params.min_portion <= 1.0, "--min-portion must be in [0, 1]; got ", params.min_portion);
  return is_valid;
}
