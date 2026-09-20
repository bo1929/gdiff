#include "roll.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <sys/mman.h>

#include "common.hpp"
#include "records.hpp"
#include "distance.hpp"
#include "llh.hpp"
#include "sketch.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "rqseq.hpp"
#include "scan.hpp"
#include "tpool.hpp"

namespace {

  // Record every k-mer's Hamming distance once, so windows can then slide.
  // Indices are relative to `base`, the first k-mer of the scanned range.
  struct block_scan_t
  {
    uint8_t* fw = nullptr;
    uint8_t* rc = nullptr;
    uint64_t base = 0;
    uint32_t hdist_th = 0;

    inline void operator()(uint64_t j, uint32_t hd, bool is_rc) const noexcept
    {
      const uint8_t v = static_cast<uint8_t>(hd <= hdist_th ? hd : hdist_th + 1);
      if (is_rc) {
        rc[j - base] = v;
      } else {
        fw[j - base] = v;
      }
    }

    inline void skip_mer(uint64_t) const noexcept {}
  };

  // Sliding per-window counts: hits 0..hdist_th plus misses, with incremental remove.
  struct sliding_counts_t
  {
    vec<uint64_t> hist_v;
    uint64_t u = 0;
    uint64_t t = 0;
    uint32_t hdist_th = 0;

    explicit sliding_counts_t(uint32_t hdist_th)
      : hist_v(hdist_bound + 1, 0)
      , hdist_th(hdist_th)
    {
    }

    void clear() noexcept
    {
      std::fill(hist_v.begin(), hist_v.end(), uint64_t(0));
      u = 0;
      t = 0;
    }

    void add(uint8_t v) noexcept
    {
      if (v == 0xFF) return;
      if (v <= hdist_th) {
        ++hist_v[v];
        ++t;
      } else {
        ++u;
      }
    }

    void remove(uint8_t v) noexcept
    {
      if (v == 0xFF) return;
      if (v <= hdist_th) {
        --hist_v[v];
        --t;
      } else {
        --u;
      }
    }

    double mle(const LLH<double>& llhf) const { return t == 0 ? nanx() : llhf.mle(hist_v.data(), u); }
  };

  // Everything one reference sketch needs to roll a batch of query sequences.
  struct roll_ctx_t
  {
    const Sketch& reference;
    const LLH<double>& llhf;
    uint32_t hdist_th;
    uint64_t tau;
    uint64_t step;
    uint64_t k;
  };

  // Fill hd_fw_v/hd_rc_v with the per-k-mer HD of the k-mers [r0, r1) (bin_shift = 0).
  template<bool Canonical>
  void
  scan_block(const scan_ctx_t& ctx, const str& seq, uint64_t r0, uint64_t r1, vec<uint8_t>& hd_fw_v, vec<uint8_t>& hd_rc_v)
  {
    block_scan_t scanner;
    scanner.fw = hd_fw_v.data();
    scanner.rc = Canonical ? nullptr : hd_rc_v.data();
    scanner.base = r0;
    scanner.hdist_th = ctx.hdist_th;
    if constexpr (Canonical) {
      scan_mers_range<true>(ctx, seq.c_str(), r0, r1, scanner);
    } else {
      scan_mers_range<false>(ctx, seq.c_str(), r0, r1, scanner);
    }
  }

  // Roll windows [w0, w1) of one sequence; window w starts at k-mer w * step.
  // A block therefore only needs the k-mers [w0 * step, (w1 - 1) * step + tau).
  template<bool Canonical>
  void roll_block(const roll_ctx_t& rctx, const qseq_t& q, uint64_t w0, uint64_t w1, std::ostream& os)
  {
    const uint64_t r0 = w0 * rctx.step;
    const uint64_t r1 = (w1 - 1) * rctx.step + rctx.tau;
    const uint64_t seg = r1 - r0;

    const scan_ctx_t sctx = make_scan_ctx(rctx.reference, 0, rctx.hdist_th);
    vec<uint8_t> hd_fw_v(seg, 0xFF);
    vec<uint8_t> hd_rc_v;
    if constexpr (!Canonical) hd_rc_v.assign(seg, 0xFF);
    scan_block<Canonical>(sctx, q.seq, r0, r1, hd_fw_v, hd_rc_v);

    sliding_counts_t counts_fw(rctx.hdist_th);
    sliding_counts_t counts_rc(rctx.hdist_th);
    for (uint64_t w = w0; w < w1; ++w) {
      const uint64_t off = (w - w0) * rctx.step; // window start relative to r0
      if (w == w0 || rctx.step >= rctx.tau) {
        // Disjoint (or the first) window: rebuilding beats remove-then-add.
        counts_fw.clear();
        counts_rc.clear();
        for (uint64_t i = off; i < off + rctx.tau; ++i) {
          counts_fw.add(hd_fw_v[i]);
          if constexpr (!Canonical) counts_rc.add(hd_rc_v[i]);
        }
      } else {
        // Overlapping windows share the middle, so only step k-mers change.
        for (uint64_t i = off - rctx.step; i < off; ++i) {
          counts_fw.remove(hd_fw_v[i]);
          if constexpr (!Canonical) counts_rc.remove(hd_rc_v[i]);
        }
        for (uint64_t i = off + rctx.tau - rctx.step; i < off + rctx.tau; ++i) {
          counts_fw.add(hd_fw_v[i]);
          if constexpr (!Canonical) counts_rc.add(hd_rc_v[i]);
        }
      }

      const uint64_t j = r0 + off;
      const double d_fw = counts_fw.mle(rctx.llhf);
      if constexpr (Canonical) {
        write_tsv(os, q.qid, j + 1, j + rctx.tau + rctx.k - 1, '.', rctx.reference.get_rname(), d_fw) << '\n';
      } else {
        const double d_rc = counts_rc.mle(rctx.llhf);
        write_tsv(os, q.qid, j + 1, j + rctx.tau + rctx.k - 1, rctx.reference.get_rname(), d_fw, d_rc) << '\n';
      }
    }
  }

  // Split every sequence's window range into contiguous, non-overlapping blocks.
  // Only the tau overlap between adjacent blocks is re-scanned (only when step < tau).
  template<bool Canonical>
  void dispatch_batch(const roll_ctx_t& ctx, const vec<qseq_t>& batch_v, ThreadPool& pool, std::ostream& os)
  {
    vec<uint64_t> nwins_v(batch_v.size(), 0);
    uint64_t total_nwins = 0;
    for (size_t b = 0; b < batch_v.size(); ++b) {
      const uint64_t len = batch_v[b].seq.size();
      if (len < ctx.k) continue;
      const uint64_t enmers = len - ctx.k + 1;
      if (enmers < ctx.tau) continue;
      nwins_v[b] = (enmers - ctx.tau) / ctx.step + 1;
      total_nwins += nwins_v[b];
    }
    if (total_nwins == 0) return;

    // Oversubscribe ~4x for load balance.
    // Keep each block long enough that its tau overlap stays a relatively small.
    const uint64_t target_tasks = 4 * static_cast<uint64_t>(pool.size());
    uint64_t nwins_per_block = (total_nwins + target_tasks - 1) / target_tasks;
    nwins_per_block = std::max(nwins_per_block, std::max<uint64_t>(1, (4 * ctx.tau) / ctx.step));

    struct block_t
    {
      uint64_t bix;
      uint64_t w0;
      uint64_t w1;
    };
    vec<block_t> blocks_v;
    blocks_v.reserve(batch_v.size() + total_nwins / nwins_per_block);
    for (size_t b = 0; b < batch_v.size(); ++b) {
      for (uint64_t w = 0; w < nwins_v[b]; w += nwins_per_block)
        blocks_v.push_back({b, w, std::min(w + nwins_per_block, nwins_v[b])});
    }

    // One buffer per block keeps output in sequence and coordinate order.
    vec<str> chunks_v(blocks_v.size());
    pool.parallel_for(blocks_v.size(), 1, [&](uint64_t bi) {
      const block_t& blk = blocks_v[bi];
      strstream ss;
      set_precision(ss, std::numeric_limits<double>::max_digits10);
      roll_block<Canonical>(ctx, batch_v[blk.bix], blk.w0, blk.w1, ss);
      chunks_v[bi] = ss.str();
    });
    for (const str& chunk : chunks_v)
      os << chunk;
  }

} // namespace

void RollSC::roll()
{
  if (is_container_file(query_path)) {
    error_exit("roll needs a query FASTA/FASTQ to roll over, not a sketch: " + query_path.string());
  }
  const Container file(sketch_path);
  if (file.size() == 0) error_exit("Container has no sketches: " + sketch_path.string());

  // One shared config per container, so the header shape is known before the loop.
  const bool canonical = file.get_config().canonical;
  std::ostream& os = *output_stream;
  write_provenance(os);
  if (canonical) {
    write_tsv(os, "seq", "start", "end", "strand", "reference", "d") << '\n';
  } else {
    write_tsv(os, "seq", "start", "end", "reference", "d_fw", "d_rc") << '\n';
  }

  const uint32_t nsketches = file.size();
  ThreadPool pool(std::max(1u, num_threads));
  cerr_msg("Rolling ",
           params.tau,
           "-mer windows in steps of ",
           params.step,
           " over ",
           query_path.filename().string(),
           " against ",
           nsketches,
           " reference sketch(es)...");

  for (uint32_t i = 0; i < nsketches; ++i) {
    const Sketch reference = file.open(i, SketchLoad::Buckets);
    const LLH<double> llhf = make_llhf(reference, params.hdist_th);
    const roll_ctx_t ctx{reference, llhf, params.hdist_th, params.tau, params.step, reference.get_k()};

    // Stream the query in batches so it never becomes fully resident.
    QSeq qs(query_path.string());
    bool more = true;
    while (more) {
      qs.clear();
      more = qs.read_next_batch();
      const vec<qseq_t>& batch_v = qs.get_batch_v();
      if (batch_v.empty()) break;
      if (canonical) {
        dispatch_batch<true>(ctx, batch_v, pool, os);
      } else {
        dispatch_batch<false>(ctx, batch_v, pool, os);
      }
    }

    const scentry& e = file.get_entry(i);
    file.advise(e.buckets_offset, e.buckets_len, MADV_DONTNEED);
    progress("Rolled", i + 1, nsketches);
  }
  progress_done();
  os.flush();
}

bool RollSC::validate_configuration()
{
  bool is_valid = true;
  if (params.hdist_th > hdist_bound) {
    cerr_msg("--hdist-th must be in [0, ", hdist_bound, "]; got ", params.hdist_th);
    is_valid = false;
  }
  if (params.step == 0) {
    cerr_msg("-s must be positive; got 0");
    is_valid = false;
  }
  if (params.tau == 0) {
    cerr_msg("-l must be positive; got 0");
    is_valid = false;
  }
  return is_valid;
}

RollSC::RollSC(CLI::App& sc)
{
  sc.add_option("query-path", query_path, "Query FASTA/FASTQ to roll the window over (gzip ok)")
    ->required()
    ->check(CLI::ExistingFile);
  sc.add_option("sketch-path", sketch_path, "Reference container <path>")->required()->check(CLI::ExistingFile);
  sc.add_option("-l", params.tau, "Window length in k-mers")->required()->check(CLI::PositiveNumber);
  sc.add_option("-s", params.step, "Step between consecutive window starts [-l]")->check(CLI::PositiveNumber);
  sc.add_option("--hdist-th", params.hdist_th, "Maximum Hamming distance for a k-mer to match [3]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.callback([&]() {
    if (!sc.count("-s")) params.step = params.tau;
    if (!validate_configuration()) error_exit("Invalid configuration!");
    open_output(output_file, output_path, output_stream);
  });
}
