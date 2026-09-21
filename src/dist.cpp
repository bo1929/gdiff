#include "dist.hpp"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <sys/mman.h>
#include <utility>
#include <unistd.h>

#include "common.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "records.hpp"
#include "scan.hpp"
#include "tpool.hpp"
#include "windows.hpp"

namespace {

  bool order_finite_first(const mle_t& lhs, const mle_t& rhs)
  {
    const bool lnan = std::isnan(lhs.d);
    const bool rnan = std::isnan(rhs.d);
    if (lnan || rnan) return !lnan && rnan;
    return lhs.d < rhs.d;
  }

  std::pair<vec<mle_t>, uint64_t> merge_samples(vec<mle_t>& ab_v, vec<mle_t>& ba_v)
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
  const auto [mle_v, n_na] = merge_samples(ab_v, ba_v);
  s.n_na = n_na;

  double bf_sum = 0.0, af_sum = 0.0;
  vec<double> bf_d_v, af_d_v;
  bf_d_v.reserve(mle_v.size());
  af_d_v.reserve(mle_v.size());

  for (const mle_t& m : mle_v) {
    if (m.s == 0.0) ++s.n_ub;

    bf_sum += m.d;
    bf_d_v.push_back(m.d);
    if (std::isnan(s.d_highest) || m.d > s.d_highest) s.d_highest = m.d;

    // A missing likelihood-ratio bound carries no evidence, so it is rejected.
    if (!std::isnan(m.s) && m.s > lr_th) {
      af_sum += m.d;
      af_d_v.push_back(m.d);
      if (std::isnan(s.d_upper) || m.d > s.d_upper) s.d_upper = m.d;
    } else {
      ++s.n_filtered;
    }
  }

  const size_t n_total = bf_d_v.size();
  const size_t n_kept = af_d_v.size();
  const double portion = n_total ? static_cast<double>(n_kept) / static_cast<double>(n_total) : 0.0;
  const double bf_mean = n_total ? bf_sum / static_cast<double>(n_total) : nanx();
  const double af_mean = n_kept ? af_sum / static_cast<double>(n_kept) : nanx();

  // The plain mean is always reported, so a caller can compare it with the filtered estimate.
  s.d_mean = bf_mean;
  if (portion > min_portion) {
    s.d = af_mean;
    s.d_v = std::move(af_d_v);
  } else {
    s.d = bf_mean;
    s.d_v = std::move(bf_d_v);
  }
  s.d_median = linear_quantile(s.d_v, 0.5);
  return s;
}

namespace {

  // One strand's match histograms for every window of a side, flattened.
  struct strand_counts_t
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
  };

  // The pool is grouped by bucket, so each bucket is resolved once per sample.
  void accumulate_pool(const Buckets& buckets,
                       const hash_pool_t& pool,
                       uint32_t hdist_th,
                       size_t nwins,
                       strand_counts_t& out) noexcept
  {
    const uint64_t* hashes = pool.hashes_ptr();
    const uint16_t* win_ix = pool.win_ix_ptr();
    const uint64_t n = pool.size();
    for (uint64_t i = 0; i < n;) {
      const uint32_t bix = key_bix(hashes[i]);
      uint64_t j = i + 1;
      while (j < n && key_bix(hashes[j]) == bix)
        ++j;
      const enc_t* beg = nullptr;
      const enc_t* end = nullptr;
      if (!buckets.range(bix, beg, end)) {
        i = j;
        continue;
      }
      __builtin_prefetch(beg, 0, 0);
      for (uint64_t e = i; e < j; ++e) {
        const size_t wix = win_ix[e];
        if (wix >= nwins) continue;
        const uint32_t hd = bucket_hdist_min(beg, end, key_enc(hashes[e]));
        if (hd <= hdist_th) ++out.window_hist(wix)[hd];
      }
      i = j;
    }
  }

  // Per-window counts for one side; rc stays empty in canonical mode.
  void
  query_windows(const Sketch& query, const Sketch& reference, uint32_t hdist_th, strand_counts_t& fw, strand_counts_t& rc)
  {
    const vec<window_t>& wins_v = query.get_windows().wins_v;
    const size_t nwins = wins_v.size();
    const bool canonical = query.is_canonical();
    fw.assign(nwins);
    if (!canonical) rc.assign(nwins);

    if (!query.get_config().keep_seq) {
      const Buckets& buckets = reference.get_buckets();
      accumulate_pool(buckets, query.get_windows().pool_fw, hdist_th, nwins, fw);
      if (!canonical) accumulate_pool(buckets, query.get_windows().pool_rc, hdist_th, nwins, rc);
      // A pool holds only retained k-mers; the rest are recovered as misses.
      for (size_t wix = 0; wix < nwins; ++wix) {
        const uint64_t matched = hist_total(fw.window_hist(wix), hdist_th);
        fw.u_v[wix] = wins_v[wix].nvalid_fw > matched ? wins_v[wix].nvalid_fw - matched : 0;
        if (!canonical) {
          const uint64_t rmatched = hist_total(rc.window_hist(wix), hdist_th);
          rc.u_v[wix] = wins_v[wix].nvalid_rc > rmatched ? wins_v[wix].nvalid_rc - rmatched : 0;
        }
      }
      return;
    }

    // scan_mers_range already counts misses, so u comes out directly.
    const scan_ctx_t ctx = make_scan_ctx(reference, 0, hdist_th);
    const seq_pack_t& packs = query.get_windows().packs;
    str cseq;
    for (size_t wix = 0; wix < nwins; ++wix) {
      packs.unpack(wix, cseq);
      const uint64_t nmers = wins_v[wix].end - wins_v[wix].start;
      if (canonical) {
        window_counts_t agg(hdist_th);
        scan_mers_range<true>(ctx, cseq.data(), 0, nmers, agg);
        std::copy(agg.hist(), agg.hist() + hdist_bound + 1, fw.window_hist(wix));
        fw.u_v[wix] = agg.u;
      } else {
        swindow_counts_t agg(hdist_th);
        scan_mers_range<false>(ctx, cseq.data(), 0, nmers, agg);
        std::copy(agg.hist_fw(), agg.hist_fw() + hdist_bound + 1, fw.window_hist(wix));
        fw.u_v[wix] = agg.u_fw;
        std::copy(agg.hist_rc(), agg.hist_rc() + hdist_bound + 1, rc.window_hist(wix));
        rc.u_v[wix] = agg.u_rc;
      }
    }
  }

} // namespace

samples_t process_samples(const Sketch& query, const Sketch& reference, uint32_t hdist_th, bool output_samples, bool is_ba)
{
  if (!compatible_configs(query.get_config(), reference.get_config())) {
    error_exit(concat_msg("Incompatible sketch pair: ",
                          query.get_rname(),
                          " vs ",
                          reference.get_rname(),
                          " (k/w/h/LSH/window length must match; use the same --seed)"));
  }
  if (!reference.has_buckets()) {
    error_exit(concat_msg("Reference sketch has no buckets: ", reference.get_rname()));
  }
  if (!query.has_windows()) {
    error_exit(concat_msg("Query sketch has no sampled windows: ", query.get_rname()));
  }

  // s depends on the reference's rho, so it differs between sides.
  const LLH<double> llhf = make_llhf(reference, hdist_th);
  const bool canonical = query.is_canonical();
  const vec<window_t>& wins_v = query.get_windows().wins_v;
  const size_t nwins = wins_v.size();

  strand_counts_t fw, rc;
  query_windows(query, reference, hdist_th, fw, rc);

  samples_t out;
  out.windows_v.reserve(nwins);
  vec<double> d_v;
  d_v.reserve(nwins);

  strstream ss;
  if (output_samples) set_precision(ss, 5);

  for (size_t wix = 0; wix < nwins; ++wix) {
    const uint64_t* hfw = fw.window_hist(wix);
    const uint64_t t_fw = hist_total(hfw, hdist_th);
    const double d_fw = t_fw == 0 ? nanx() : llhf.mle(hfw, fw.u_v[wix]);

    mle_t mle;
    mle.d = d_fw;
    char strand = canonical ? '.' : '+';
    const uint64_t* hist = hfw;
    uint64_t u = fw.u_v[wix];
    if (!canonical) {
      const uint64_t* hrc = rc.window_hist(wix);
      const uint64_t t_rc = hist_total(hrc, hdist_th);
      const double d_rc = t_rc == 0 ? nanx() : llhf.mle(hrc, rc.u_v[wix]);
      const auto picked = select_strand_distance(d_fw, d_rc);
      mle.d = picked.first;
      strand = picked.second;
      if (strand == '-') {
        hist = hrc;
        u = rc.u_v[wix];
      }
    }

    if (is_valid_distance(mle.d)) {
      mle.s = compute_lr_ub(llhf, mle.d, u + hist_total(hist, hdist_th));
      d_v.push_back(mle.d);
    }
    out.windows_v.push_back(mle);

    if (output_samples) {
      write_tsv(ss,
                "gdiff",
                is_ba ? reference.get_rname() : query.get_rname(),
                is_ba ? query.get_rname() : reference.get_rname(),
                wins_v[wix].qid,
                wins_v[wix].start + 1,
                wins_v[wix].end + query.get_k() - 1,
                strand,
                is_ba ? "ba" : "ab",
                mle.d,
                mle.s)
        << '\n';
    }
  }

  std::sort(d_v.begin(), d_v.end());
  out.d_median = linear_quantile(d_v, 0.5);
  out.n_valid = d_v.size();
  if (output_samples) out.samples = ss.str();
  return out;
}

constexpr uint64_t hcompl_bp_min = 20ull * 1000 * 1000; // 20 Mbp of valid bp
constexpr uint64_t hcompl_hdist_th = 2u;

namespace {
  void print_pair_progress(uint64_t done_jobs, uint64_t total_jobs)
  {
    if (!stderr_is_tty()) return; // keep redirected logs clean
    constexpr int bar_width = 40;
    const double frac = total_jobs ? static_cast<double>(done_jobs) / static_cast<double>(total_jobs) : 1.0;
    const int filled = static_cast<int>(std::lround(frac * bar_width));
    std::ostringstream os;
    os << "\rprogress: [" << std::setfill('#') << std::setw(static_cast<int>(filled)) << "" << std::setfill('.')
       << std::setw(bar_width - filled) << "" << std::setfill(' ') << "] " << std::setw(3)
       << static_cast<int>(std::lround(frac * 100.0)) << "% (" << std::fixed << std::setprecision(1)
       << 0.5 * static_cast<double>(done_jobs) << "/" << 0.5 * static_cast<double>(total_jobs) << ")" << std::flush;
    std::cerr << os.str();
  }

} // namespace

const Container* DistSC::container_for(const std::filesystem::path& path)
{
  for (const auto& f : containers_v) {
    if (f->get_path() == path) return f.get();
  }
  containers_v.push_back(std::make_unique<Container>(path));
  return containers_v.back().get();
}

void DistSC::resolve_container(const std::filesystem::path& path, vec<uint32_t>& out)
{
  const Container* file = container_for(path);
  size_t first = entries_v.size();
  for (size_t i = 0; i < entries_v.size(); ++i) {
    if (entries_v[i].file == file) {
      first = i;
      break;
    }
  }
  if (first == entries_v.size()) {
    for (uint32_t rix = 0; rix < file->size(); ++rix)
      entries_v.push_back({file, rix, str{}});
  }
  for (uint32_t rix = 0; rix < file->size(); ++rix)
    out.push_back(static_cast<uint32_t>(first + rix));
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

void DistSC::emit_header(std::ostream& os) const
{
  if (params.output_samples) {
    write_tsv(os, "config", "genome_a", "genome_b", "seq", "start", "end", "strand", "direction", "d", "lr_ub") << '\n';
  } else {
    write_tsv(os,
              "genome_a",
              "genome_b",
              "d",
              "d_median",
              "d_mean",
              "d_upper",
              "d_highest",
              "d_ab",
              "d_ba",
              "n_ab",
              "n_ba",
              "n_ub",
              "n_na",
              "n_filtered")
      << '\n';
  }
}

void DistSC::write_pair_line(std::ostream& os, const pair_t& pr, const str& name_a, const str& name_b)
{
  const summary_t s = summarize_symmetric(pr.ab.windows_v, pr.ba.windows_v, params.lr_th, params.min_portion);
  write_tsv(os,
            name_a,
            name_b,
            s.d,
            s.d_median,
            s.d_mean,
            s.d_upper,
            s.d_highest,
            pr.ab.d_median,
            pr.ba.d_median,
            pr.ab.n_valid,
            pr.ba.n_valid,
            s.n_ub,
            s.n_na,
            s.n_filtered)
    << '\n';
}

void DistSC::estimate_distances()
{
  if (set_a_v.empty()) error_exit("Set A is empty; nothing to compare");

  auto dedup = [](vec<uint32_t>& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  };
  dedup(set_a_v);
  dedup(set_b_v);

  const bool is_within = set_b_v.empty();
  if (is_within) {
    cerr_msg("entries=", entries_v.size(), " within-mode a=", set_a_v.size());
  } else {
    cerr_msg("entries=", entries_v.size(), " cross-mode a=", set_a_v.size(), " b=", set_b_v.size());
  }

  if (entries_v.front().file->get_config().tau == 0) {
    error_exit("Sketch has no sampled windows; re-run `gdiff sketch` appropriately");
  }

  // Every entry is a query in some pair, so every entry needs its window payload.
  ThreadPool pool(std::max(1u, num_threads));
  vec<Sketch> targets_v(entries_v.size());
  pool.parallel_for(entries_v.size(), 1, [&](uint64_t mi) {
    const entry_t& m = entries_v[static_cast<size_t>(mi)];
    targets_v[static_cast<size_t>(mi)] = m.file->open(m.rix, SketchLoad::Windows);
  });
  for (size_t mi = 0; mi < entries_v.size(); ++mi)
    entries_v[mi].rname = targets_v[mi].get_rname();

  // Within mode compares every unordered pair of one set; cross mode the full product.
  const vec<uint32_t>& refs_v = is_within ? set_a_v : set_b_v;
  vec<pair_t> pairs_v;
  pairs_v.reserve(is_within ? set_a_v.size() * (set_a_v.size() - 1) / 2 : set_a_v.size() * set_b_v.size());
  for (size_t i = 0; i < set_a_v.size(); ++i) {
    for (size_t j = is_within ? i + 1 : 0; j < refs_v.size(); ++j) {
      pairs_v.push_back({set_a_v[i], refs_v[j], {}, {}});
    }
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
  std::ostream& os = *output_stream;
  write_provenance(os);
  emit_header(os);
  if (pairs_v.empty()) {
    cerr_msg("No pairs to compare");
    os.flush();
    return;
  }

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
    while (batch_end < jobs_v.size() && jobs_v[batch_end].source_ix == source_ix) {
      ++batch_end;
    }

    const entry_t& sentry = entries_v[source_ix];
    const Sketch source = sentry.file->open(sentry.rix, SketchLoad::Buckets);

    // One threshold per source: its length caps whatever was asked for.
    const uint64_t nvalid_bp = source.get_nvalid_bp();
    const uint32_t hdist_cap = nvalid_bp >= hcompl_bp_min ? hcompl_hdist_th : params.hdist_th;
    const uint32_t hdist_th = std::min(params.hdist_th, hdist_cap);
    if (hdist_th != params.hdist_th) {
      warn_pmsg(source.get_rname(), nvalid_bp, " valid bp: --hdist-th is set to ", hdist_th, " (was ", params.hdist_th, ")");
    }

    pool.parallel_for(batch_end - batch_start, 1, [&](uint64_t i) {
      const job_t& job = jobs_v[batch_start + static_cast<size_t>(i)];
      samples_t r = process_samples(targets_v[job.target_ix], source, hdist_th, params.output_samples, job.is_ba);
      pair_t& pr = pairs_v[job.pair_ix];
      if (job.is_ba)
        pr.ba = std::move(r);
      else
        pr.ab = std::move(r);
    });

    const scentry& e = sentry.file->get_entry(sentry.rix);
    sentry.file->advise(e.buckets_offset, e.buckets_len, MADV_DONTNEED);

    done_jobs += batch_end - batch_start;
    print_pair_progress(done_jobs, total_jobs);
    batch_start = batch_end;
  }
  if (stderr_is_tty()) std::cerr << '\n';

  if (params.output_samples) {
    for (const pair_t& pr : pairs_v) {
      os << pr.ab.samples << pr.ba.samples;
    }
  } else {
    set_precision(os, 8);
    for (const pair_t& pr : pairs_v)
      write_pair_line(os, pr, entries_v[pr.a].rname, entries_v[pr.b].rname);
  }
  os.flush();
}

void DistSC::dist()
{
  const bool is_a_positional = !container_a_path.empty();
  const bool is_b_positional = !container_b_path.empty();
  const bool is_a_list = !list_a_path.empty();
  const bool is_b_list = !list_b_path.empty();

  if (is_a_positional && is_a_list) {
    error_exit("Set A given both positionally and via --list-a; pick one");
  }
  if (is_b_positional && is_b_list) {
    error_exit("Set B given both positionally and via --list-b; pick one");
  }
  if (!is_a_positional && !is_a_list) {
    error_exit("No set A given: pass a sketch container, or --list-a");
  }

  if (is_a_list) {
    resolve_list(list_a_path, set_a_v);
  } else if (is_container_file(container_a_path)) {
    resolve_container(container_a_path, set_a_v);
  } else {
    error_exit("Not a sketch container: " + container_a_path.string() + "; pass a list of containers with --list-a");
  }

  if (is_b_list) {
    resolve_list(list_b_path, set_b_v);
  } else if (is_b_positional) {
    if (is_container_file(container_b_path)) {
      resolve_container(container_b_path, set_b_v);
    } else {
      error_exit("Not a sketch container: " + container_b_path.string() + "; pass a list of containers with --list-b");
    }
  }
  estimate_distances();
}

DistSC::DistSC(CLI::App& sc)
{
  sc.add_option("sketch-a", container_a_path, "Sketch container for set A (omit when using --list-a)")
    ->check(CLI::ExistingFile);
  sc.add_option("sketch-b", container_b_path, "Sketch container for set B (omit for every pair within set A)")
    ->check(CLI::ExistingFile);
  sc.add_option("--list-a", list_a_path, "Set A as a list file of sketch containers")
    ->excludes("sketch-a")
    ->check(CLI::ExistingFile);
  sc.add_option("--list-b", list_b_path, "Set B as a list file of sketch containers")
    ->excludes("sketch-b")
    ->check(CLI::ExistingFile);
  sc.add_option(
      "--hdist-th", params.hdist_th, "Maximum Hamming distance for k-mer search; capped at 2 for inputs >20 Mbp [3]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("--lr-th", params.lr_th, "Likelihood-ratio cut for the reconciliation filter [10.828]")
    ->check(CLI::NonNegativeNumber);
  sc.add_option("--min-portion",
                params.min_portion,
                "Apply the filter iff at least --min-portion of the windows exceeds --lr-th [0.66]")
    ->check(CLI::Range(0.0, 1.0));
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.add_flag("--output-samples", params.output_samples, "Write per-window sample lines instead of per-pair summaries");
  sc.callback([&]() {
    if (!validate_configuration()) error_exit("Invalid configuration!");
    open_output(output_file, output_path, output_stream);
  });
}

bool DistSC::validate_configuration()
{
  if (params.hdist_th > hdist_bound) {
    cerr_msg("--hdist-th must be in [0, ", hdist_bound, "]; got ", params.hdist_th);
    return false;
  }
  return true;
}
