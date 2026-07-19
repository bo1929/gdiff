#include "dist.hpp"

#include "common.hpp"
#include "enc.hpp"
#include "msg.hpp"
#include "random.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>

namespace {

  const auto dist_url_validator = CLI::Validator(
    [](std::string& input) { return match_url(input) ? std::string("") : "Given URL is not valid: " + input; },
    "URL",
    "URL validator");

  double linear_quantile(const vec<double>& values, const double probability)
  {
    if (values.empty()) return nanx();
    if (values.size() == 1) return values.front();

    const double index = probability * static_cast<double>(values.size() - 1);
    const size_t lo = static_cast<size_t>(std::floor(index));
    const size_t hi = static_cast<size_t>(std::ceil(index));
    const double fraction = index - static_cast<double>(lo);
    return values[lo] + fraction * (values[hi] - values[lo]);
  }

} // namespace

HDHistogram::HDHistogram(const uint64_t nmers, const uint32_t hdist_th)
  : nmers(nmers)
  , hdist_th(hdist_th)
  , width(hdist_th + 1)
  , hist_v((nmers + 1) * width, 0)
{
}

void HDHistogram::aggregate_mer(const uint32_t hdist, const uint64_t pos)
{
  if (hdist <= hdist_th && pos < nmers) {
    ++hist_v[((pos + 1) * width) + hdist];
  }
}

void HDHistogram::inclusive_scan()
{
  for (uint64_t pos = 1; pos <= nmers; ++pos) {
    for (uint32_t d = 0; d < width; ++d) {
      hist_v[(pos * width) + d] += hist_v[((pos - 1) * width) + d];
    }
  }
}

void HDHistogram::extract_histogram(const uint64_t a,
                                    const uint64_t b,
                                    vec<uint64_t>& hist,
                                    uint64_t& misses,
                                    uint64_t& hits) const
{
  assert(a <= b);
  assert(b <= nmers);

  hist.assign(hdist_bound + 1, 0);
  hits = 0;
  for (uint32_t d = 0; d < width; ++d) {
    hist[d] = hist_v[(b * width) + d] - hist_v[(a * width) + d];
    hits += hist[d];
  }
  misses = (b - a) - hits;
}

dist_summary_t summarize_distances(vec<double> distances)
{
  dist_summary_t summary;
  distances.erase(std::remove_if(distances.begin(), distances.end(), [](const double d) { return !std::isfinite(d); }),
                  distances.end());
  if (distances.empty()) {
    summary.quantiles.fill(nanx());
    return summary;
  }

  std::sort(distances.begin(), distances.end());
  summary.n = distances.size();
  summary.mean = std::accumulate(distances.begin(), distances.end(), 0.0) / static_cast<double>(summary.n);

  if (summary.n > 1) {
    double sum_sq = 0.0;
    for (const double d : distances) {
      const double delta = d - summary.mean;
      sum_sq += delta * delta;
    }
    summary.sd = std::sqrt(sum_sq / static_cast<double>(summary.n - 1));
  } else {
    summary.sd = 0.0;
  }

  constexpr arr<double, 7> probabilities{0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99};
  for (size_t i = 0; i < probabilities.size(); ++i) {
    summary.quantiles[i] = linear_quantile(distances, probabilities[i]);
  }
  return summary;
}

vec<uint64_t>
sample_region_starts(const uint64_t seq_len, const uint64_t region_len, const uint64_t sample_size, std::mt19937& rng)
{
  vec<uint64_t> starts;
  if (region_len > seq_len) return starts;

  starts.reserve(sample_size);
  std::uniform_int_distribution<uint64_t> start_dist(0, seq_len - region_len);
  for (uint64_t i = 0; i < sample_size; ++i) {
    starts.push_back(start_dist(rng));
  }
  return starts;
}

std::pair<double, char> select_strand_distance(const double d_fw, const double d_rc)
{
  const bool fw_valid = std::isfinite(d_fw);
  const bool rc_valid = std::isfinite(d_rc);
  if (!fw_valid && !rc_valid) return {nanx(), '.'};
  if (!rc_valid || (fw_valid && d_fw <= d_rc)) return {d_fw, '+'};
  return {d_rc, '-'};
}

DistSC::DistSC(CLI::App& sc)
{
  sc.add_option("-q,--query-path", query_path, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible)")
    ->required()
    ->check(dist_url_validator | CLI::ExistingFile);
  sc.add_option("-i,--sketch-path", sketch_path, "Sketch file at <path> to query")->required()->check(CLI::ExistingFile);
  sc.add_option("-l,--length", region_len, "Length of sampled query regions")->required()->check(CLI::PositiveNumber);
  sc.add_option("--sample-size", sample_size, "Number of regions sampled per query/reference pair [200]")
    ->check(CLI::PositiveNumber);
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match [4]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("-o,--output-path", output_path, "Write summary output to a file at <path> [stdout]");
  sc.add_option("--samples-output", samples_output_path, "Write sampled regions and distances to a TSV file");
  sc.callback([&]() {
    if (!output_path.empty() && !samples_output_path.empty() && output_path == samples_output_path) {
      error_exit("--output-path and --samples-output must be different files");
    }
    if (!output_path.empty()) {
      output_file.open(output_path);
      check_fstream(output_file, "Cannot open output file", output_path.string());
      output_stream = &output_file;
    }
    if (!samples_output_path.empty()) {
      samples_output_file.open(samples_output_path);
      check_fstream(samples_output_file, "Cannot open samples output file", samples_output_path.string());
      samples_output_stream = &samples_output_file;
    }
  });
}

void DistSC::dist()
{
  qseq_sptr_t queries = std::make_shared<QSeq>(query_path);
  while (queries->read_next_batch()) {
  }

  std::ifstream sketch_stream(sketch_path, std::ifstream::binary);
  check_fstream(sketch_stream, "Cannot open sketch file", sketch_path.string());

  uint32_t nsketches = 0;
  sketch_stream.read(reinterpret_cast<char*>(&nsketches), sizeof(uint32_t));
  for (uint32_t i = 0; i < nsketches; ++i) {
    const uint64_t offset = static_cast<uint64_t>(sketch_stream.tellg());
    Sketch::seek_past(sketch_stream);
    const uint64_t next_offset = static_cast<uint64_t>(sketch_stream.tellg());

    sketch_sptr_t sketch = std::make_shared<Sketch>(sketch_path);
    sketch->load_from_offset(sketch_stream, offset);
    sketch->make_rho_partial();

    const auto& seq_batch = queries->get_seq_batch();
    const auto& qid_batch = queries->get_qid_batch();
    for (size_t qix = 0; qix < seq_batch.size(); ++qix) {
      process_pair(sketch, seq_batch[qix], qid_batch[qix]);
    }
    sketch_stream.clear();
    sketch_stream.seekg(static_cast<std::streamoff>(next_offset));
  }
}

void DistSC::process_pair(const sketch_sptr_t& sketch, const str& seq, const str& qid)
{
  const lshf_sptr_t lshf = sketch->get_lshf();
  const uint32_t k = lshf->get_k();
  if (region_len < k) {
    error_exit(concat_msg("--length must be at least the sketch k-mer length (length=", region_len, ", k=", k, ")"));
  }
  if (seq.size() < region_len) {
    warn_pmsg(qid, "skipped for ", sketch->get_rid(), ": sequence shorter than requested region length");
    return;
  }

  const uint64_t nmers = seq.size() - k + 1;
  HDHistogram hist_fw(nmers, hdist_th);
  std::unique_ptr<HDHistogram> hist_rc;
  if (!sketch->is_canonical()) hist_rc = std::make_unique<HDHistogram>(nmers, hdist_th);
  search_mers(sketch, seq.data(), seq.size(), hist_fw, hist_rc.get());
  hist_fw.inclusive_scan();
  if (hist_rc) hist_rc->inclusive_scan();

  LLH<double> llhf(k, lshf->get_h(), sketch->get_rho(), hdist_th, 0.0, false);
  const uint64_t region_nmers = region_len - k + 1;
  const vec<uint64_t> starts = sample_region_starts(seq.size(), region_len, sample_size, gen);
  vec<double> distances;
  distances.reserve(starts.size());
  vec<uint64_t> hist;

  for (const uint64_t start : starts) {
    uint64_t misses = 0, hits = 0;
    hist_fw.extract_histogram(start, start + region_nmers, hist, misses, hits);
    const double d_fw = validate_distance(llhf.mle(hist.data(), misses));

    double distance = d_fw;
    char strand = '+';
    if (hist_rc) {
      hist_rc->extract_histogram(start, start + region_nmers, hist, misses, hits);
      const double d_rc = validate_distance(llhf.mle(hist.data(), misses));
      std::tie(distance, strand) = select_strand_distance(d_fw, d_rc);
    }
    if (!std::isfinite(distance)) continue;

    distances.push_back(distance);
    write_sample(qid, seq.size(), start, strand, sketch->get_rid(), distance);
  }

  write_summary(qid, sketch->get_rid(), summarize_distances(std::move(distances)));
}

void DistSC::search_mers(const sketch_sptr_t& sketch,
                         const char* cseq,
                         const uint64_t len,
                         HDHistogram& hist_fw,
                         HDHistogram* hist_rc) const
{
  const lshf_sptr_t lshf = sketch->get_lshf();
  const uint32_t k = lshf->get_k();
  const uint64_t u64m = std::numeric_limits<uint64_t>::max();
  const uint64_t mask_lr = ((u64m >> (64 - k)) << 32) + ((u64m << 32) >> (64 - k));
  const uint64_t mask_bp = u64m >> ((32 - k) * 2);

  uint64_t valid_len = 0;
  uint64_t enc_bp = 0, enc_lr = 0;
  for (uint64_t i = 0; i < len; ++i) {
    if (__builtin_expect(SEQ_NT4_TABLE[cseq[i]] >= 4, 0)) {
      valid_len = 0;
      continue;
    }
    ++valid_len;
    if (valid_len < k) continue;

    const uint64_t pos = i - k + 1;
    if (valid_len == k) {
      compute_encoding(cseq + pos, cseq + i + 1, enc_lr, enc_bp);
    } else {
      update_encoding(cseq + i, enc_lr, enc_bp);
    }
    enc_bp &= mask_bp;
    enc_lr &= mask_lr;
    const uint64_t rc_bp = revcomp_bp64(enc_bp, k);

    auto search = [&](const uint64_t bp, const uint64_t lr, HDHistogram& histogram) {
      const uint32_t rix = lshf->compute_hash(bp);
      const uint32_t offset = sketch->partial_offset(rix);
      sketch->prefetch_offset_inc(offset);
      const enc_t enc_lr_compact = lshf->drop_ppos_lr(lr);
      sketch->prefetch_offset_enc(offset);
      uint32_t hdist = 0;
      if (sketch->scan_bucket(offset, enc_lr_compact, hdist)) histogram.aggregate_mer(hdist, pos);
    };

    if (hist_rc) {
      search(enc_bp, enc_lr, hist_fw);
      search(rc_bp, bp64_to_lr64(rc_bp), *hist_rc);
    } else if (rc_bp < enc_bp) {
      search(enc_bp, enc_lr, hist_fw);
    } else {
      search(rc_bp, bp64_to_lr64(rc_bp), hist_fw);
    }
  }
}

void DistSC::write_summary(const str& qid, const str& rid, const dist_summary_t& summary)
{
  *output_stream << std::setprecision(8);
  write_tsv(*output_stream, qid, rid, summary.n, summary.mean, summary.sd);
  for (const double quantile : summary.quantiles) {
    *output_stream << '\t' << quantile;
  }
  *output_stream << '\n';
}

void DistSC::write_sample(const str& qid,
                          const uint64_t seq_len,
                          const uint64_t start,
                          const char strand,
                          const str& rid,
                          const double distance)
{
  if (!samples_output_stream) return;
  *samples_output_stream << std::setprecision(8);
  write_tsv(*samples_output_stream, qid, seq_len, start + 1, start + region_len, strand, rid, distance) << '\n';
}
