#include "gdiff.hpp"

void MergeSC::merge()
{
  cerr_msg("Preparing to merge ", paths_v.size(), " sketch file(s)");

  std::ofstream sout(sketch_path, std::ofstream::binary);
  check_fstream(sout, "Cannot open output file", sketch_path);

  // Writing a placeholder
  uint32_t total_sketches = 0;
  sout.write(reinterpret_cast<const char*>(&total_sketches), sizeof(uint32_t));

  constexpr size_t buffer_size = 10 * 1024 * 1024;
  std::vector<char> buffer(buffer_size);

  for (size_t i = 0; i < paths_v.size(); ++i) {
    std::ifstream sin;
    sin.rdbuf()->pubsetbuf(buffer.data(), buffer_size);
    sin.open(paths_v[i], std::ifstream::binary);
    check_fstream(sin, "Cannot open sketch file", paths_v[i]);

    uint32_t nsketches = 0;
    sin.read(reinterpret_cast<char*>(&nsketches), sizeof(uint32_t));
    total_sketches += nsketches;

    sout << sin.rdbuf();
    sin.close();
  }

  // Seek back and patch the real count into the header
  sout.seekp(0, std::ios::beg);
  sout.write(reinterpret_cast<const char*>(&total_sketches), sizeof(uint32_t));

  check_fstream(sout, "Failed to write the merged sketch file!", sketch_path);
  sout.close();

  cerr_msg("Merged sketch saved to ", sketch_path.string(), " with ", total_sketches, " sketch(es)");
}

MergeSC::MergeSC(CLI::App& sc)
{
  sc.add_option("-i,--sketch-paths", paths_v, "Input sketch files to merge")->required()->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", sketch_path, "Path to store the merged sketch file")->required();
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
    std::string rname(rid_len, '\0');
    stream.read(&rname[0], static_cast<std::streamsize>(rid_len));

    uint64_t timestamp = 0;
    stream.read(reinterpret_cast<char*>(&timestamp), sizeof(uint64_t));

    uint8_t k = 0, w = 0, h = 0;
    uint32_t nrows = 0;
    bool canonical = false;
    stream.read(reinterpret_cast<char*>(&k), sizeof(uint8_t));
    stream.read(reinterpret_cast<char*>(&w), sizeof(uint8_t));
    stream.read(reinterpret_cast<char*>(&h), sizeof(uint8_t));
    stream.read(reinterpret_cast<char*>(&canonical), sizeof(bool));
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
    std::cout << "  Name:        " << rname << "\n";
    std::cout << "  Date:        " << ts_str << "\n";
    std::cout << "  k (mer len): " << static_cast<int>(k) << "\n";
    std::cout << "  w (win len): " << static_cast<int>(w) << "\n";
    std::cout << "  h (LSH pos): " << static_cast<int>(h) << "\n";
    std::cout << "  canonical:   " << (canonical ? "true" : "false") << "\n";
    std::cout << "  nrows:       " << nrows << "\n";
    std::cout << "  frac:        " << static_cast<double>(nrows) / static_cast<double>(uint64_t(1) << (2 * h)) << "\n";
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
  app.add_option("--num-threads", num_threads, "Number of threads for parallel sketch/map/dist processing [1]");

  auto& sc_sketch = *app.add_subcommand("sketch", "Create sketches from FASTA/FASTQ files");
  auto& sc_map = *app.add_subcommand("map", "Map queries and extract distance-based patterns from sketches");
  auto& sc_dist = *app.add_subcommand("dist", "Sample query regions and summarize MLE distances");
  auto& sc_detect = *app.add_subcommand("detect", "Fit a null distance distribution and detect outlier regions");
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
