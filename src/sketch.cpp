#include "sketch.hpp"

#include <atomic>
#include <chrono>
#include <mutex>
#include <thread>
#include <utility>

#include "msg.hpp"
#include "random.hpp"
#include "rqseq.hpp"

extern uint32_t num_threads;

Sketch::Sketch(std::filesystem::path sketch_path)
  : sketch_path(std::move(sketch_path))
{
}

void Sketch::load_from_offset(std::ifstream& stream, uint64_t offset)
{
  if (offset > 0) {
    stream.clear();
    stream.seekg(offset);
    check_fstream(stream, "Failed to seek in the sketch file", sketch_path);
  }

  uint64_t rid_len;
  stream.read(reinterpret_cast<char*>(&rid_len), sizeof(uint64_t));
  rid.resize(rid_len);
  stream.read(&rid[0], rid_len);
  stream.read(reinterpret_cast<char*>(&timestamp), sizeof(uint64_t));

  stream.read(reinterpret_cast<char*>(&k), sizeof(uint8_t));
  stream.read(reinterpret_cast<char*>(&w), sizeof(uint8_t));
  stream.read(reinterpret_cast<char*>(&h), sizeof(uint8_t));
  stream.read(reinterpret_cast<char*>(&canonical), sizeof(bool));
  stream.read(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));

  vec<uint8_t> ppos_v(h), npos_v(k - h);
  stream.read(reinterpret_cast<char*>(ppos_v.data()), h * sizeof(uint8_t));
  stream.read(reinterpret_cast<char*>(npos_v.data()), (k - h) * sizeof(uint8_t));

  lshf = std::make_shared<LSHF>(ppos_v, npos_v);

  stream.read(reinterpret_cast<char*>(&rho), sizeof(double));

  sfhm = std::make_shared<SFHM>();
  sfhm->load(stream);

  check_fstream(stream, "Failed to read the sketch file!", sketch_path);
}

void Sketch::seek_past(std::ifstream& stream)
{
  uint64_t rid_len;
  stream.read(reinterpret_cast<char*>(&rid_len), sizeof(uint64_t));
  // Skip rid (rid_len bytes) + timestamp (8 bytes)
  stream.seekg(static_cast<std::streamoff>(rid_len) + static_cast<std::streamoff>(sizeof(uint64_t)), std::ios::cur);

  uint8_t k, h;
  stream.read(reinterpret_cast<char*>(&k), sizeof(uint8_t));
  // skip: w (1), h is read next, then canonical (1), nrows (4) = 6 bytes
  stream.seekg(sizeof(uint8_t), std::ios::cur);
  stream.read(reinterpret_cast<char*>(&h), sizeof(uint8_t));
  // skip: canonical (1) + nrows (4) = 5 bytes
  stream.seekg(sizeof(bool) + sizeof(uint32_t), std::ios::cur);
  // skip: ppos_v (h bytes) + npos_v ((k-h) bytes) = k bytes total
  stream.seekg(static_cast<std::streamoff>(k), std::ios::cur);
  // skip: rho (8 bytes)
  stream.seekg(static_cast<std::streamoff>(sizeof(double)), std::ios::cur);

  uint64_t nkmers;
  stream.read(reinterpret_cast<char*>(&nkmers), sizeof(uint64_t));
  stream.seekg(static_cast<std::streamoff>(nkmers) * static_cast<std::streamoff>(sizeof(enc_t)), std::ios::cur);
  uint32_t sfhm_nrows;
  stream.read(reinterpret_cast<char*>(&sfhm_nrows), sizeof(uint32_t));
  stream.seekg(static_cast<std::streamoff>(sfhm_nrows) * static_cast<std::streamoff>(sizeof(inc_t)), std::ios::cur);
}

double Sketch::get_rho() const { return rho; }

void Sketch::prefetch_offset_inc(uint32_t offset) const noexcept
{
  if (offset != OFF_INVALID) {
    sfhm->prefetch_inc(offset);
  }
}

void Sketch::prefetch_offset_enc(uint32_t offset) const noexcept
{
  // inc_v must already be in cache for this to be effective
  if (offset != OFF_INVALID) {
    sfhm->prefetch_enc(offset);
  }
}

bool Sketch::scan_bucket(uint32_t offset, enc_t enc_lr, uint32_t& hdist_min) const noexcept
{
  if (offset == OFF_INVALID) return false;
  const enc_t* ix1 = sfhm->bucket_ptr_start(offset);
  const enc_t* ix2 = sfhm->bucket_ptr_next(offset);
  uint32_t hmin = std::numeric_limits<uint32_t>::max();
  for (; ix1 < ix2; ++ix1) {
    const uint32_t hd = popcount_lr32((*ix1) ^ enc_lr);
    hmin = hd < hmin ? hd : hmin;
  }
  hdist_min = hmin;
  return true;
}

void Sketch::canonicalize()
{
  if (canonical) {
    return;
  }
  const uint64_t mask_bp = std::numeric_limits<uint64_t>::max() >> ((32 - k) * 2);
  sdhm_sptr_t sdhm = std::make_shared<SDHM>();
  sdhm->enc_vvec.resize(nrows);
  for (uint32_t off = 0; off < nrows; ++off) {
    const enc_t* ix1 = sfhm->bucket_ptr_start(off);
    const enc_t* ix2 = sfhm->bucket_ptr_next(off);
    for (; ix1 < ix2; ++ix1) {
      const uint64_t bp_ppos = lshf->inv_ppos_bp(off);
      const uint64_t bp_npos = lr64_to_bp64(lshf->inv_ppos_lr(*ix1));
      const uint64_t fw_bp = (bp_ppos | bp_npos) & mask_bp;
      const uint64_t rc_bp = revcomp_bp64(fw_bp, k);
      const uint64_t can_bp = std::max(fw_bp, rc_bp);
      const uint32_t rixn = lshf->compute_hash(can_bp);
      if (rixn >= nrows) continue;
      const enc_t new_enc = lshf->drop_ppos_lr(bp64_to_lr64(can_bp));
      sdhm->enc_vvec[rixn].push_back(new_enc);
    }
  }
  sdhm->sort_columns();
  sdhm->make_unique();
  sfhm = std::make_shared<SFHM>(sdhm);
  canonical = true;
}

void BaseLSH::set_lshf() { lshf = std::make_shared<LSHF>(k, h); }

void BaseLSH::set_nrows()
{
  // Flat sampling: keep the prefix [0, T) of the 2^(2h) LSH space.
  const uint64_t hash_size = uint64_t(1) << (2 * h);
  const uint64_t t = static_cast<uint64_t>(hash_size * frac + 0.5);
  nrows = static_cast<uint32_t>(std::max<uint64_t>(1, std::min(t, hash_size)));
}

bool SketchSC::validate_configuration()
{
  bool is_invalid = false;
  if (frac <= 0.0 || frac > 1.0) {
    is_invalid = true;
    cerr_msg("--frac must be in (0, 1]; got ", frac);
  }
  if (w < k) {
    is_invalid = true;
    cerr_msg("The minimum minimizer window size (-w) is k (-k)!");
  }
  if (h < 3) {
    is_invalid = true;
    cerr_msg("The minimum number of LSH positions (-h) is 3!");
  }
  if (h > 15) {
    is_invalid = true;
    cerr_msg("The maximum number of LSH positions (-h) is 15!");
  }
  if (k > 31) {
    is_invalid = true;
    cerr_msg("The maximum allowed k-mer length (-k) is 31!");
  }
  if (k < 19) {
    is_invalid = true;
    cerr_msg("The minimum allowed k-mer length (-k) is 19!");
  }
  if ((k - h) > 16) {
    is_invalid = true;
    cerr_msg("For compact k-mer encodings, h must be >= k-16!");
  }
  return !is_invalid;
}

void SketchSC::process()
{
  if (input_paths.empty()) error_exit("No input files provided!");

  const uint32_t nsketches = static_cast<uint32_t>(input_paths.size());
  const uint32_t nthreads = std::max(1u, std::min(num_threads, nsketches));
  cerr_msg("Preparing to sketch ", nsketches, " file(s) w/ ", nthreads, " thread(s)");

  std::ofstream sketch_stream(sketch_path, std::ofstream::binary);
  sketch_stream.write(reinterpret_cast<const char*>(&nsketches), sizeof(uint32_t));
  rho_v.assign(nsketches, 0.0);

  std::atomic<uint32_t> next_idx{0};
  std::atomic<uint32_t> count_p{0};
  std::mutex write_mtx;
  std::mutex cerr_mtx;

  auto worker = [&](const uint32_t tseed) {
    init_thread_rng(tseed);
    uint32_t i;
    while ((i = next_idx.fetch_add(1, std::memory_order_relaxed)) < nsketches) {
      const str& input_path = input_paths[i];
      rseq_sptr_t rs = std::make_shared<RSeq>(input_path, lshf, w, nrows, canonical);
      sdhm_sptr_t sdhm = std::make_shared<SDHM>();
      sdhm->fill_table(nrows, rs);
      sfhm_sptr_t sketch_sfhm = std::make_shared<SFHM>(sdhm);
      rho_v[i] = rs->get_rho();

      {
        std::lock_guard<std::mutex> lock(write_mtx);
        write_header(sketch_stream, i);
        write_config(sketch_stream, i);
        sketch_sfhm->save(sketch_stream);
      }

      const uint32_t num_p = count_p.fetch_add(1, std::memory_order_relaxed) + 1;
      {
        std::lock_guard<std::mutex> lock(cerr_mtx);
        std::cerr << "\rCreated a sketch for [" << num_p << "/" << nsketches << "] "
                  << "..." << std::flush;
        if (num_p == nsketches) std::cerr << std::endl;
      }
    }
  };

  std::vector<std::thread> threads;
  threads.reserve(nthreads);
  for (uint32_t t = 0; t < nthreads; ++t) {
    threads.emplace_back([&, t]() { worker(t + 1); });
  }
  for (auto& t : threads) {
    t.join();
  }

  check_fstream(sketch_stream, std::string("Failed to write the sketch!"), sketch_path.string());
  sketch_stream.close();
  cerr_msg("Sketch file saved to ", sketch_path.string(), " with ", nsketches, " sketch(es)");
}

void SketchSC::write_header(std::ofstream& sout, uint32_t i)
{
  const str rid = std::filesystem::path(input_paths[i]).filename().string();
  uint64_t timestamp =
    std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  uint64_t rid_len = rid.length();
  sout.write(reinterpret_cast<char*>(&rid_len), sizeof(uint64_t));
  sout.write(rid.c_str(), rid_len);
  sout.write(reinterpret_cast<char*>(&timestamp), sizeof(uint64_t));
}

void SketchSC::write_config(std::ofstream& sout, uint32_t i)
{
  // Keeps a minimizer iff LSH(x) < nrows; rho already includes that keep rate.
  sout.write(reinterpret_cast<const char*>(&k), sizeof(uint8_t));
  sout.write(reinterpret_cast<const char*>(&w), sizeof(uint8_t));
  sout.write(reinterpret_cast<const char*>(&h), sizeof(uint8_t));
  sout.write(reinterpret_cast<char*>(&canonical), sizeof(bool));
  sout.write(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));
  sout.write(reinterpret_cast<char*>(lshf->ppos_data()), h * sizeof(uint8_t));
  sout.write(reinterpret_cast<char*>(lshf->npos_data()), (k - h) * sizeof(uint8_t));
  sout.write(reinterpret_cast<const char*>(&rho_v[i]), sizeof(double));
}

SketchSC::SketchSC(CLI::App& sc)
{
  set_sketch_defaults();
  sc.add_option("-i,--input-path", input_paths, "Input FASTA/FASTQ file(s) <path> (or URL) (gzip compatible)")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("-o,--output-path", sketch_path, "Path to store the resulting binary sketch file")->required();
  sc.add_option("-k,--mer-len", k, "Length of k-mers [27]")->check(CLI::Range(19, 31))->check(CLI::PositiveNumber);
  sc.add_option("-w,--win-len", w, "Length of the minimizer window (w>=k) [k+6]")->check(CLI::PositiveNumber);
  sc.add_option("-h,--num-positions", h, "Number of positions for the LSH [k-16]")->check(CLI::PositiveNumber);
  sc.add_option("--frac", frac, "Keep a k-mer if LSH(x) < frac * 2^(2h); i.e., subsampling ratio [1.0]")
    ->check(CLI::Range(std::numeric_limits<double>::min(), 1.0));
  sc.add_flag(
    "--strand-agnostic,!--strand-aware", canonical, "A (canonical) strand-agnostic (default) or strand-aware sketch");
  sc.callback([&]() {
    if (!(sc.count("-w") + sc.count("--win-len"))) {
      w = k + 6;
      h = k - 16;
    }
    if (!validate_configuration()) {
      error_exit("Invalid configuration!");
    }
  });
}
