#include "gdiff.hpp"
#include "random.hpp"

void BaseLSH::set_lshf() { lshf = std::make_shared<LSHF>(k, h, m); }

void BaseLSH::set_nrows()
{
  uint32_t hash_size = 1u << (2 * h);
  uint32_t full_residue = hash_size % m;
  if (frac) {
    nrows = (hash_size / m) * (r + 1);
    nrows = full_residue > r ? nrows + (r + 1) : nrows + full_residue;
  } else {
    nrows = (hash_size / m);
    nrows = full_residue > r ? nrows + 1 : nrows;
  }
}

bool MapSC::validate_configuration()
{
  bool is_invalid = false;
  if (dist_th.size() != 1 && dist_th.size() != 8) {
    is_invalid = true;
    cerr_msg("--dist-th requires exactly 1 or 8 thresholds; got ", dist_th.size());
  }
  for (size_t i = 0; i < dist_th.size(); ++i) {
    if (dist_th[i] <= 0.0) {
      is_invalid = true;
      cerr_msg("--dist-th[", i, "] must be positive: ", dist_th[i]);
    }
  }
  {
    auto sdist_th = dist_th;
    std::sort(sdist_th.begin(), sdist_th.end());
    if (const auto it = std::adjacent_find(sdist_th.begin(), sdist_th.end()); it != sdist_th.end()) {
      is_invalid = true;
      cerr_msg("--dist-th values must be unique; duplicate: ", *it);
    }
  }
  if (hdist_th > hdist_bound) {
    is_invalid = true;
    cerr_msg("--hdist-th must be in [0, ", hdist_bound, "] with the current SIMD histogram layout; got ", hdist_th);
  }
  if (bin_shift >= 63) {
    is_invalid = true;
    cerr_msg("--bin-shift must be less than 63; got ", bin_shift);
  }
  const uint64_t bin_size = (bin_shift < 63) ? (uint64_t(1) << bin_shift) : 0;
  if (bin_size > tau) {
    is_invalid = true;
    cerr_msg("--bin-shift gives bin_size=", bin_size, ", which exceeds --min-length=", tau);
  }
  const uint64_t tau_bin = (bin_size > 0) ? ((tau + bin_size - 1) >> bin_shift) : 0;
  if (tau_bin < 2) {
    is_invalid = true;
    cerr_msg("--min-length must span at least two bins after binning ", "(tau=", tau, ", bin_size=", bin_size, ")");
  }
  return !is_invalid;
}

void MapSC::write_header()
{
  (*output_stream)
    << "QUERY_ID\tSEQ_LEN\tINTERVAL_START\tINTERVAL_END\tSTRAND\tIS_RC\tREF_ID\tDIST\tMASK\tD_INTERVAL\tDIST_CONTIG\tSTRAND_DIFF\tDIST_GENOME\tPERCENTILE\tFOLD\tQVALUE\n";
}

void MapSC::map()
{
  *(output_stream) << std::setprecision(5);
  write_header();

  // Load all query sequences once? Might be inefficient
  qseq_sptr_t qs = std::make_shared<QSeq>(query_path);

  bool cont_reading;
  while ((cont_reading = qs->read_next_batch())) {
    total_qseq += qs->get_cbatch_size();
  }
  total_qseq += qs->get_cbatch_size();

  std::ifstream sketch_stream(sketch_path, std::ifstream::binary);
  check_fstream(sketch_stream, std::string("Cannot open sketch file: "), sketch_path.string());

  uint32_t nsketches;
  sketch_stream.read(reinterpret_cast<char*>(&nsketches), sizeof(uint32_t));
  const uint32_t nthreads = std::max(1u, std::min(num_threads, nsketches));
  cerr_msg("Processing ", nsketches, " sketches w/ ", nthreads, " thread(s)...");

  std::vector<uint64_t> sketch_offsets(nsketches);
  for (uint32_t i = 0; i < nsketches; ++i) {
    sketch_offsets[i] = static_cast<uint64_t>(sketch_stream.tellg());
    Sketch::seek_past(sketch_stream); // reads headers, seeks over data
  }
  sketch_stream.close();

  size_t n = dist_th.size();
  params_t<double> params_single(n, dist_th.front(), hdist_th, tau, chisq, bin_shift, sample_size, enum_only);
  params_t<cm512_t> params_multiple(n, {0}, hdist_th, tau, chisq, bin_shift, sample_size, enum_only);
  std::copy(dist_th.begin(), dist_th.end(), params_multiple.dist_th.begin());

  // Per-sketch result buffers
  std::vector<strstream> results(nsketches);
  std::atomic<uint32_t> next_idx{0};
  std::atomic<uint32_t> done_count{0};
  std::mutex cerr_mtx;

  auto worker = [&](const uint32_t tseed) {
    init_thread_rng(tseed);
    uint32_t i;
    while ((i = next_idx.fetch_add(1, std::memory_order_relaxed)) < nsketches) {
      // Each worker opens its own file handle so no stream sharing occurs
      std::ifstream sketch_stream(sketch_path, std::ifstream::binary);
      sketch_sptr_t sketch = std::make_shared<Sketch>(sketch_path);
      sketch->load_from_offset(sketch_stream, sketch_offsets[i]);
      sketch_stream.close();
      sketch->make_rho_partial();

      strstream sout;
      sout << std::setprecision(5);
      if (dist_th.size() == 1) {
        QIE<double> qie(params_single, sketch, sketch->get_lshf(), qs->get_seq_batch(), qs->get_qid_batch());
        qie.map_sequences(sout, sketch->get_rid());
      } else {
        QIE<cm512_t> qie(params_multiple, sketch, sketch->get_lshf(), qs->get_seq_batch(), qs->get_qid_batch());
        qie.map_sequences(sout, sketch->get_rid());
      }

      // Store result at its reserved slot (no aliasing between threads)
      results[i] = std::move(sout);

      uint32_t done = done_count.fetch_add(1, std::memory_order_relaxed) + 1;
      {
        std::lock_guard<std::mutex> lock(cerr_mtx);
        std::cerr << "\rProcessed sketch " << done << "/" << nsketches << "..." << std::flush;
        if (done == nsketches) std::cerr << std::endl;
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

  for (uint32_t i = 0; i < nsketches; ++i) {
    if (results[i].tellp() > 0) *(output_stream) << results[i].rdbuf();
  }
}

bool SketchSC::validate_configuration()
{
  bool is_invalid = false;
  if (w < k) {
    is_invalid = true;
    cerr_msg("The minimum minimizer window size (-w) is k (-k)!");
  }
  if (h < 3) {
    is_invalid = true;
    cerr_msg("The minimum number of LSH positions (-h) is 3!");
  }
  if (h > 16) {
    is_invalid = true;
    cerr_msg("The maximum number of LSH positions (-h) is 16!");
  }
  if (k > 32) {
    is_invalid = true;
    cerr_msg("The maximum allowed k-mer length (-k) is 32!");
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

void SketchSC::create()
{
  rseq_sptr_t rs = std::make_shared<RSeq>(input_path, lshf, w, r, frac);
  sdhm_sptr_t sdhm = std::make_shared<SDHM>();
  sdhm->fill_table(nrows, rs);
  sketch_sfhm = std::make_shared<SFHM>(sdhm);
  rho = rs->get_rho();

  cerr_msg("Total number of k-mers included in the sketch: ", sdhm->get_nmers());
  cerr_msg("Subsampling rate (rho) is: ", rho);
}

void SketchSC::write_header(std::ofstream& sout)
{
  const str& rid = sketch_path.filename();
  uint64_t timestamp =
    std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
  uint64_t rid_len = rid.length();
  sout.write(reinterpret_cast<char*>(&rid_len), sizeof(uint64_t));
  sout.write(rid.c_str(), rid_len);
  sout.write(reinterpret_cast<char*>(&timestamp), sizeof(uint64_t));
}

void SketchSC::write_config(std::ofstream& sout)
{
  sout.write(reinterpret_cast<char*>(&k), sizeof(uint8_t));
  sout.write(reinterpret_cast<char*>(&w), sizeof(uint8_t));
  sout.write(reinterpret_cast<char*>(&h), sizeof(uint8_t));
  sout.write(reinterpret_cast<char*>(&m), sizeof(uint32_t));
  sout.write(reinterpret_cast<char*>(&r), sizeof(uint32_t));
  sout.write(reinterpret_cast<char*>(&frac), sizeof(bool));
  sout.write(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));
  sout.write(reinterpret_cast<char*>(lshf->ppos_data()), h * sizeof(uint8_t));
  sout.write(reinterpret_cast<char*>(lshf->npos_data()), (k - h) * sizeof(uint8_t));
  sout.write(reinterpret_cast<char*>(&rho), sizeof(double));
}

void SketchSC::save()
{
  std::ofstream sketch_stream(sketch_path, std::ofstream::binary);
  uint32_t nsketches = 1;
  sketch_stream.write(reinterpret_cast<char*>(&nsketches), sizeof(uint32_t));

  write_header(sketch_stream);
  write_config(sketch_stream);
  sketch_sfhm->save(sketch_stream);

  check_fstream(sketch_stream, std::string("Failed to write the sketch!"), sketch_path.string());
  sketch_stream.close();
}

SketchSC::SketchSC(CLI::App& sc)
{
  set_sketch_defaults();
  sc.add_option("-i,--input-path", input_path, "Input FASTA/FASTQ file <path> (or URL) (gzip compatible)")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("-o,--output-path", sketch_path, "Path to store the resulting binary sketch file")->required();
  sc.add_option("-k,--mer-len", k, "Length of k-mers [27]")->check(CLI::Range(19, 32))->check(CLI::PositiveNumber);
  sc.add_option("-w,--win-len", w, "Length of the minimizer window (w>=k) [k+6]")->check(CLI::PositiveNumber);
  sc.add_option("-h,--num-positions", h, "Number of positions for the LSH [k-16]")->check(CLI::PositiveNumber);
  sc.add_option("-m,--modulo-lsh", m, "Modulo value to partition LSH space [2]")->check(CLI::PositiveNumber);
  sc.add_option("-r,--residue-lsh", r, "A k-mer x will be included only if r = LSH(x) mod m [1]")
    ->check(CLI::NonNegativeNumber);
  sc.add_flag("--frac,!--no-frac", frac, "Include k-mers with r <= LSH(x) mod m [true]");
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

MapSC::MapSC(CLI::App& sc)
{
  sc.add_option("-q,--query-path", query_path, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible)")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("-i,--sketch-path", sketch_path, "Sketch file at <path> to query")->required()->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path> [stdout]");
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match [4]")
    ->check(CLI::Range(0, static_cast<int>(hdist_bound)));
  sc.add_option("--chisq", chisq, "Chi-square threshold [33.00051]")->check(CLI::NonNegativeNumber);
  sc.add_option("-d,--dist-th", dist_th, "Distance threshold(s) - provide exactly 1 or 8 values")->required()->expected(1, 8);
  sc.add_option("-l,--min-length", tau, "Minimum interval length")->required()->check(CLI::PositiveNumber);
  sc.add_option("-b,--bin-shift", bin_shift, "Group consecutive k-mers into bins of size 2^b [0]")->check(CLI::Range(0, 62));
  sc.add_flag("--enum-only,!--no-enum-only", enum_only, "Enumerate intervals without MLE distance estimation [false]");
  sc.add_option("--sample-size", sample_size, "Samples for significance test (0: skip) [200]")->check(CLI::NonNegativeNumber);
  sc.callback([&]() {
    if (!validate_configuration()) {
      error_exit("Invalid configuration!");
    }
    if (!output_path.empty()) {
      output_file.open(output_path);
      output_stream = &output_file;
    }
  });
}

void MergeSC::merge()
{
  cerr_msg("Preparing to merge ", sketch_paths.size(), " sketch file(s)");

  std::ofstream sout(output_path, std::ofstream::binary);
  check_fstream(sout, "Cannot open output file", output_path);

  // Writing a placeholder
  uint32_t total_sketches = 0;
  sout.write(reinterpret_cast<const char*>(&total_sketches), sizeof(uint32_t));

  constexpr size_t buffer_size = 10 * 1024 * 1024;
  std::vector<char> buffer(buffer_size);

  for (size_t i = 0; i < sketch_paths.size(); ++i) {
    std::ifstream sin;
    sin.rdbuf()->pubsetbuf(buffer.data(), buffer_size);
    sin.open(sketch_paths[i], std::ifstream::binary);
    check_fstream(sin, "Cannot open sketch file", sketch_paths[i]);

    uint32_t nsketches = 0;
    sin.read(reinterpret_cast<char*>(&nsketches), sizeof(uint32_t));
    total_sketches += nsketches;

    sout << sin.rdbuf();
    sin.close();
  }

  // Seek back and patch the real count into the header
  sout.seekp(0, std::ios::beg);
  sout.write(reinterpret_cast<const char*>(&total_sketches), sizeof(uint32_t));

  check_fstream(sout, "Failed to write the merged sketch file!", output_path);
  sout.close();

  cerr_msg("Merged sketch saved to ", output_path.string(), " with ", total_sketches, " sketch(es)");
}

MergeSC::MergeSC(CLI::App& sc)
{
  sc.add_option("-i,--sketch-paths", sketch_paths, "Input sketch files to merge")->required()->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", output_path, "Path to store the merged sketch file")->required();
}

void InfoSC::info()
{
  std::ifstream stream(sketch_path, std::ifstream::binary);
  check_fstream(stream, "Cannot open sketch file: ", sketch_path.string());

  uint32_t nsketches = 0;
  stream.read(reinterpret_cast<char*>(&nsketches), sizeof(uint32_t));

  std::cout << "File:             " << sketch_path.string() << "\n";
  std::cout << "Number of sketches: " << nsketches << "\n";

  for (uint32_t i = 0; i < nsketches; ++i) {
    uint64_t rid_len = 0;
    stream.read(reinterpret_cast<char*>(&rid_len), sizeof(uint64_t));
    std::string rid(rid_len, '\0');
    stream.read(&rid[0], static_cast<std::streamsize>(rid_len));

    uint64_t timestamp = 0;
    stream.read(reinterpret_cast<char*>(&timestamp), sizeof(uint64_t));

    uint8_t k = 0, w = 0, h = 0;
    uint32_t m = 0, r = 0, nrows = 0;
    bool frac = false;
    stream.read(reinterpret_cast<char*>(&k), sizeof(uint8_t));
    stream.read(reinterpret_cast<char*>(&w), sizeof(uint8_t));
    stream.read(reinterpret_cast<char*>(&h), sizeof(uint8_t));
    stream.read(reinterpret_cast<char*>(&m), sizeof(uint32_t));
    stream.read(reinterpret_cast<char*>(&r), sizeof(uint32_t));
    stream.read(reinterpret_cast<char*>(&frac), sizeof(bool));
    stream.read(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));

    stream.seekg(static_cast<std::streamoff>(h) + static_cast<std::streamoff>(k - h), std::ios::cur);

    double rho = 0.0;
    stream.read(reinterpret_cast<char*>(&rho), sizeof(double));

    uint64_t nkmers = 0;
    stream.read(reinterpret_cast<char*>(&nkmers), sizeof(uint64_t));
    stream.seekg(static_cast<std::streamoff>(nkmers) * static_cast<std::streamoff>(sizeof(enc_t)), std::ios::cur);
    uint32_t sfhm_nrows = 0;
    stream.read(reinterpret_cast<char*>(&sfhm_nrows), sizeof(uint32_t));
    stream.seekg(static_cast<std::streamoff>(sfhm_nrows) * static_cast<std::streamoff>(sizeof(inc_t)), std::ios::cur);

    std::time_t ts = static_cast<std::time_t>(timestamp);
    std::string ts_str = std::ctime(&ts);
    if (!ts_str.empty() && ts_str.back() == '\n') ts_str.pop_back();

    std::cout << "\n[Sketch " << (i + 1) << "/" << nsketches << "]\n";
    std::cout << "  Name:        " << rid << "\n";
    std::cout << "  Date:        " << ts_str << "\n";
    std::cout << "  k (mer len): " << static_cast<int>(k) << "\n";
    std::cout << "  w (win len): " << static_cast<int>(w) << "\n";
    std::cout << "  h (LSH pos): " << static_cast<int>(h) << "\n";
    std::cout << "  m (modulo):  " << m << "\n";
    std::cout << "  r (residue): " << r << "\n";
    std::cout << "  frac:        " << (frac ? "true" : "false") << "\n";
    std::cout << "  nrows:       " << nrows << "\n";
    std::cout << "  rho:         " << rho << "\n";
    std::cout << "  k-mers:      " << nkmers << "\n";
  }

  stream.close();
}

InfoSC::InfoSC(CLI::App& sc)
{
  sc.add_option("-i,--sketch-path", sketch_path, "Sketch file (single or multi) to inspect")
    ->required()
    ->check(CLI::ExistingFile);
}

int main(int argc, char** argv)
{
  PRINT_VERSION
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  CLI::App app{"gdiff"};
  app.set_help_flag("--help");
  app.fallthrough();

  bool verbose = false;
  app.add_flag("--verbose,!--no-verbose", verbose, "Increased verbosity and progress report");
  app.require_subcommand();
  app.add_option("--seed", seed, "Random seed for the LSH and other parts that require randomness [0]");
  app.callback([&]() { init_thread_rng(0); });
  app.add_option("--num-threads", num_threads, "Number of threads for parallel sketch processing [1]");

  auto& sc_sketch = *app.add_subcommand("sketch", "Create sketches from FASTA/FASTQ files");
  auto& sc_map = *app.add_subcommand("map", "Map queries and extract distance-based patterns from sketches");
  auto& sc_merge = *app.add_subcommand("merge", "Merge multiple sketches into a single sketch file");
  auto& sc_info = *app.add_subcommand("info", "Show metadata for all sketches in a sketch file");

  SketchSC gdiff_sketch(sc_sketch);
  MapSC gdiff_map(sc_map);
  MergeSC gdiff_merge(sc_merge);
  InfoSC gdiff_info(sc_info);

  CLI11_PARSE(app, argc, argv);
  for (int i = 0; i < argc; ++i) {
    invocation += str(argv[i]) + " ";
  }
  if (!invocation.empty()) {
    invocation.pop_back();
  }

  auto tstart = std::chrono::system_clock::now();
  std::time_t tstart_f = std::chrono::system_clock::to_time_t(tstart);
  str invocation_str = "Invocation: " + invocation + "\n";
  std::cerr << invocation_str;
  std::cerr << std::ctime(&tstart_f);

  if (sc_sketch.parsed()) {
    cerr_msg("Initializing the sketch...");
    gdiff_sketch.set_nrows();
    gdiff_sketch.set_lshf();
    std::chrono::duration<float> es_b = std::chrono::system_clock::now() - tstart;
    gdiff_sketch.create();
    gdiff_sketch.save();
    std::chrono::duration<float> es_s = std::chrono::system_clock::now() - tstart - es_b;
    cerr_msg("Done sketching & saving, elapsed: ", es_s.count(), " sec");
  }

  if (sc_merge.parsed()) {
    cerr_msg("Merging sketches...");
    std::chrono::duration<float> es_b = std::chrono::system_clock::now() - tstart;
    gdiff_merge.merge();
    std::chrono::duration<float> es_s = std::chrono::system_clock::now() - tstart - es_b;
    cerr_msg("Done merging sketches, elapsed: ", es_s.count(), " sec");
  }

  if (sc_map.parsed()) {
    cerr_msg("Loading the sketch...");
    cerr_msg("Seeking query sequences in the sketch...");
    std::chrono::duration<float> es_b = std::chrono::system_clock::now() - tstart;
    gdiff_map.map();
    std::chrono::duration<float> es_s = std::chrono::system_clock::now() - tstart - es_b;
    cerr_msg("Done mapping sequences, elapsed: ", es_s.count(), " sec");
    cerr_msg("Total number of sequences queried: ", gdiff_map.get_total_qseq());
  }

  if (sc_info.parsed()) {
    gdiff_info.info();
  }

  auto tend = std::chrono::system_clock::now();
  std::time_t tend_f = std::chrono::system_clock::to_time_t(tend);
  std::cerr << std::ctime(&tend_f);

  return 0;
}
