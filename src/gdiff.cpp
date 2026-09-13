#include "gdiff.hpp"

void MergeSC::merge()
{
  cerr_msg("Merging ", paths_v.size(), " sketch file(s)");

  // Re-index rather than concatenate: records move, absolute offsets are rewritten.
  vec<std::unique_ptr<SketchFile>> inputs;
  uint64_t nsketches = 0;
  for (const str& path : paths_v) {
    inputs.push_back(std::make_unique<SketchFile>(path));
    const SketchFile& in = *inputs.back();
    if (!in.get_config().compatible_with(inputs.front()->get_config())) {
      error_exit("Cannot merge sketches with different configurations: " + path);
    }
    if (in.get_config().win_repr != inputs.front()->get_config().win_repr) {
      error_exit("Cannot merge sketches with different --window-repr: " + path);
    }
    nsketches += in.size();
  }
  if (nsketches == 0) error_exit("Nothing to merge: the inputs hold no records");

  std::ofstream sout(sketch_path, std::ofstream::binary);
  check_fstream(sout, "Cannot open output file", sketch_path.string());
  write_sketch_header(sout, inputs.front()->get_config(), nsketches);
  const std::streampos index_pos = sout.tellp();
  vec<record_entry_t> entries(nsketches);
  sout.write(reinterpret_cast<const char*>(entries.data()),
             static_cast<std::streamsize>(sizeof(record_entry_t) * nsketches));

  vec<char> buffer;
  uint64_t pos = static_cast<uint64_t>(sout.tellp());
  uint64_t out_ix = 0;
  for (size_t fi = 0; fi < inputs.size(); ++fi) {
    std::ifstream sin(paths_v[fi], std::ifstream::binary);
    check_fstream(sin, "Cannot open sketch file", paths_v[fi]);
    for (const record_entry_t& e : inputs[fi]->get_index().records) {
      buffer.resize(static_cast<size_t>(e.len));
      sin.seekg(static_cast<std::streamoff>(e.offset));
      sin.read(buffer.data(), static_cast<std::streamsize>(e.len));
      check_fstream(sin, "Failed to read a sketch record", paths_v[fi]);

      record_entry_t& out = entries[out_ix++];
      out.offset = pos;
      out.len = e.len;
      out.buckets_off = pos + (e.buckets_off - e.offset);
      out.buckets_len = e.buckets_len;
      if (e.windows_len) {
        out.windows_off = pos + (e.windows_off - e.offset);
        out.windows_len = e.windows_len;
      }
      sout.write(buffer.data(), static_cast<std::streamsize>(e.len));
      pos += e.len;
    }
  }

  sout.seekp(index_pos);
  sout.write(reinterpret_cast<const char*>(entries.data()),
             static_cast<std::streamsize>(sizeof(record_entry_t) * nsketches));
  sout.seekp(0, std::ios::end);
  check_fstream(sout, "Failed to write the merged sketch file", sketch_path.string());
  sout.close();

  cerr_msg("Merged sketch saved to ", sketch_path.string(), " with ", nsketches, " record(s)");
}

MergeSC::MergeSC(CLI::App& sc)
{
  sc.add_option("-i,--sketch-paths", paths_v, "Input sketch files to merge")->required()->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", sketch_path, "Path to store the merged sketch file")->required();
}

void InfoSC::info()
{
  const SketchFile file(sketch_path);
  const sketch_config_t& cfg = file.get_config();

  std::cout << "File:               " << sketch_path.string() << "\n";
  std::cout << "Records:            " << file.size() << "\n";
  std::cout << "k (mer len):        " << static_cast<int>(cfg.k) << "\n";
  std::cout << "w (win len):        " << static_cast<int>(cfg.w) << "\n";
  std::cout << "h (LSH pos):        " << static_cast<int>(cfg.h) << "\n";
  std::cout << "canonical:          " << (cfg.canonical ? "true" : "false") << "\n";
  std::cout << "nrows:              " << cfg.nrows << "\n";
  std::cout << "frac:               " << cfg.frac << "\n";
  std::cout << "-l (window len):    " << cfg.tau << "\n";
  std::cout << "--sample-size:      " << cfg.sample_size << "\n";
  std::cout << "--window-repr:      " << (cfg.win_repr == WinRepr::Pool ? "pool" : "seq") << "\n";
  std::cout << "seed:               " << cfg.seed << "\n";

  for (uint32_t i = 0; i < file.size(); ++i) {
    const Sketch sk = file.open(i, SketchPart::All);
    std::time_t ts = static_cast<std::time_t>(sk.get_timestamp());
    str ts_str = std::ctime(&ts);
    if (!ts_str.empty() && ts_str.back() == '\n') ts_str.pop_back();

    std::cout << "\n[Record " << (i + 1) << "/" << file.size() << "]\n";
    std::cout << "  Name:             " << sk.get_rname() << "\n";
    std::cout << "  Date:             " << ts_str << "\n";
    std::cout << "  Genome length:    " << sk.get_genome_bp() << "\n";
    std::cout << "  Valid bases:      " << sk.get_nvalid_bases() << "\n";
    std::cout << "  k-mers:           " << sk.get_nkmers() << "\n";
    std::cout << "  card_est:         " << sk.get_card_est() << "\n";
    std::cout << "  rho:              " << sk.get_rho() << "\n";
    std::cout << "  Nonempty buckets: " << sk.get_buckets().get_nnonempty() << "\n";
    std::cout << "  Windows:          " << sk.get_wins().size() << "\n";
    const record_entry_t& e = file.get_index().records[i];
    std::cout << "  Bucket bytes:     " << e.buckets_len << "\n";
    std::cout << "  Window bytes:     " << e.windows_len << "\n";
  }
}

InfoSC::InfoSC(CLI::App& sc)
{
  sc.add_option("-i,--sketch-path", sketch_path, "Sketch file to inspect")
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
  app.add_option("--num-threads", num_threads, "Number of threads for parallel sketch/map/dist processing [1]");

  auto& sc_sketch = *app.add_subcommand("sketch", "Create sketches from FASTA/FASTQ files");
  auto& sc_map = *app.add_subcommand("map", "Map queries and extract distance-based patterns from sketches");
  auto& sc_dist = *app.add_subcommand("dist", "Sample query regions and summarize MLE distances");
  auto& sc_detect = *app.add_subcommand("detect", "Fit a background distance distribution and detect outlier regions");
  auto& sc_merge = *app.add_subcommand("merge", "Merge multiple sketches into a single sketch file");
  auto& sc_info = *app.add_subcommand("info", "Show metadata for all sketches in a sketch file");

  SketchSC gdiff_sketch(sc_sketch);
  MapSC gdiff_map(sc_map);
  DistSC gdiff_dist(sc_dist);
  DetectSC gdiff_detect(sc_detect);
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

  auto run_timed = [&](const char* done_msg, auto&& work) {
    const auto t0 = std::chrono::system_clock::now();
    work();
    std::chrono::duration<float> es = std::chrono::system_clock::now() - t0;
    cerr_msg(done_msg, es.count(), " sec");
  };

  if (sc_sketch.parsed()) {
    cerr_msg("Initializing the sketch...");
    gdiff_sketch.set_nrows();
    gdiff_sketch.set_lshf();
    run_timed("Done sketching & saving, elapsed: ", [&]() { gdiff_sketch.process(); });
  }

  if (sc_merge.parsed()) {
    cerr_msg("Merging sketches...");
    run_timed("Done merging sketches, elapsed: ", [&]() { gdiff_merge.merge(); });
  }

  if (sc_map.parsed()) {
    cerr_msg("Loading the sketch...");
    cerr_msg("Seeking query sequences in the sketch...");
    run_timed("Done mapping sequences, elapsed: ", [&]() { gdiff_map.map(); });
    cerr_msg("Total number of sequences queried: ", gdiff_map.get_total_qseq());
  }

  if (sc_dist.parsed()) {
    cerr_msg("Sampling query regions and calculating distances...");
    run_timed("Done calculating distances, elapsed: ", [&]() { gdiff_dist.dist(); });
  }

  if (sc_detect.parsed()) {
    cerr_msg("Sampling distances and detecting outlier regions...");
    run_timed("Done detecting outlier regions, elapsed: ", [&]() { gdiff_detect.detect(); });
  }

  if (sc_info.parsed()) {
    gdiff_info.info();
  }

  auto tend = std::chrono::system_clock::now();
  std::time_t tend_f = std::chrono::system_clock::to_time_t(tend);
  std::cerr << std::ctime(&tend_f);

  return 0;
}
