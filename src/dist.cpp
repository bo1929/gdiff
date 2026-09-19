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
#include "tsv.hpp"

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
    write_tsv(os, "config", "genome_a", "genome_b", "seq", "start", "end", "strand", "direction", "d", "lr_bg", "lr_ub")
      << '\n';
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

void DistSC::write_pair_row(std::ostream& os, const pair_t& pr, const str& name_a, const str& name_b)
{
  const summary_t s = summarize_symmetric(pr.ab.rows_v, pr.ba.rows_v, params.lr_th, params.min_portion);
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
      direction_t r = run_direction(targets_v[job.target_ix], source, hdist_th, params.output_samples, job.is_ba);
      pair_t& pr = pairs_v[job.pair_ix];
      if (job.is_ba)
        pr.ba = std::move(r);
      else
        pr.ab = std::move(r);
    });

    const sketch_entry_t& e = sentry.file->get_entry(sentry.rix);
    sentry.file->advise(e.buckets_off, e.buckets_len, MADV_DONTNEED);

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
      write_pair_row(os, pr, entries_v[pr.a].rname, entries_v[pr.b].rname);
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
      "--hdist-th", params.hdist_th, "Hamming distance limit for a k-mer match; if input is >20 Mbp, it is capped at 2 [3]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("--lr-th", params.lr_th, "Likelihood-ratio cut for the reconciliation filter [3.841]")
    ->check(CLI::NonNegativeNumber);
  sc.add_option("--min-portion",
                params.min_portion,
                "Apply the filter only if at least this fraction of the windows exceeds --lr-th [0.66]")
    ->check(CLI::Range(0.0, 1.0));
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.add_flag("--output-samples", params.output_samples, "Write per-window sample rows instead of per-pair summaries");
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
