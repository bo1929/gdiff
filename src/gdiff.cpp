#include "gdiff.hpp"
#include "roll.hpp"

void MergeSC::merge()
{
  cerr_msg("Merging ", paths_v.size(), " container(s)");

  // Re-index rather than concatenate: sketches move, absolute offsets are rewritten.
  vec<std::unique_ptr<Container>> inputs_v;
  uint64_t nsketches = 0;
  for (const str& path : paths_v) {
    inputs_v.push_back(std::make_unique<Container>(path));
    const Container& in = *inputs_v.back();
    if (!compatible_configs(in.get_config(), inputs_v.front()->get_config())) {
      error_exit("Cannot merge sketches with different configurations: " + path);
    }
    if (in.get_config().keep_seq != inputs_v.front()->get_config().keep_seq) {
      error_exit("Cannot merge sketches with different --keep-seq: " + path);
    }
    nsketches += in.size();
  }
  if (nsketches == 0) error_exit("Nothing to merge: the inputs hold no sketches");

  std::ofstream sout(sketch_path, std::ofstream::binary);
  check_fstream(sout, "Cannot open output file", sketch_path.string());
  write_container_header(sout, inputs_v.front()->get_config(), nsketches);
  const std::streampos index_pos = sout.tellp();
  vec<sketch_entry_t> entries_v(nsketches);
  sout.write(reinterpret_cast<const char*>(entries_v.data()),
             static_cast<std::streamsize>(sizeof(sketch_entry_t) * nsketches));

  vec<char> buffer_v;
  uint64_t pos = static_cast<uint64_t>(sout.tellp());
  uint64_t out_ix = 0;
  for (size_t fi = 0; fi < inputs_v.size(); ++fi) {
    std::ifstream sin(paths_v[fi], std::ifstream::binary);
    check_fstream(sin, "Cannot open container", paths_v[fi]);
    for (const sketch_entry_t& e : inputs_v[fi]->get_entries()) {
      buffer_v.resize(static_cast<size_t>(e.len));
      sin.seekg(static_cast<std::streamoff>(e.offset));
      sin.read(buffer_v.data(), static_cast<std::streamsize>(e.len));
      check_fstream(sin, "Failed to read a sketch", paths_v[fi]);

      sketch_entry_t& out = entries_v[out_ix++];
      out.offset = pos;
      out.len = e.len;
      out.buckets_off = pos + (e.buckets_off - e.offset);
      out.buckets_len = e.buckets_len;
      if (e.windows_len) {
        out.windows_off = pos + (e.windows_off - e.offset);
        out.windows_len = e.windows_len;
      }
      sout.write(buffer_v.data(), static_cast<std::streamsize>(e.len));
      pos += e.len;
    }
  }

  sout.seekp(index_pos);
  sout.write(reinterpret_cast<const char*>(entries_v.data()),
             static_cast<std::streamsize>(sizeof(sketch_entry_t) * nsketches));
  sout.seekp(0, std::ios::end);
  check_fstream(sout, "Failed to write the merged container", sketch_path.string());
  sout.close();

  cerr_msg("Merged sketch saved to ", sketch_path.string(), " with ", nsketches, " sketch(es)");
}

MergeSC::MergeSC(CLI::App& sc)
{
  sc.add_option("-i,--sketch-paths", paths_v, "Input containers to merge")->required()->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", sketch_path, "Path to store the merged container")->required();
}

void InfoSC::info()
{
  const Container file(sketch_path);
  const sketch_config_t& cfg = file.get_config();

  std::cout << "File:               " << sketch_path.string() << "\n";
  std::cout << "Sketches:           " << file.size() << "\n";
  std::cout << "k (mer len):        " << static_cast<int>(cfg.k) << "\n";
  std::cout << "w (win len):        " << static_cast<int>(cfg.w) << "\n";
  std::cout << "h (LSH pos):        " << static_cast<int>(cfg.h) << "\n";
  std::cout << "canonical:          " << (cfg.canonical ? "true" : "false") << "\n";
  std::cout << "nrows:              " << cfg.nrows << "\n";
  std::cout << "frac:               " << cfg.frac << "\n";
  std::cout << "-l (window len):    " << cfg.tau << "\n";
  std::cout << "--sample-size:      " << cfg.sample_size << "\n";
  std::cout << "--keep-seq:         " << (cfg.keep_seq ? "true" : "false") << "\n";
  std::cout << "seed:               " << cfg.seed << "\n";

  for (uint32_t i = 0; i < file.size(); ++i) {
    const Sketch sk = file.open(i, SketchLoad::All);
    std::time_t ts = static_cast<std::time_t>(sk.get_timestamp());
    str ts_str = std::ctime(&ts);
    if (!ts_str.empty() && ts_str.back() == '\n') ts_str.pop_back();

    std::cout << "\n[Sketch " << (i + 1) << "/" << file.size() << "]\n";
    std::cout << "  Name:             " << sk.get_rname() << "\n";
    std::cout << "  Date:             " << ts_str << "\n";
    std::cout << "  Genome length:    " << sk.get_ntotal_bp() << "\n";
    std::cout << "  Valid bases:      " << sk.get_nvalid_bp() << "\n";
    std::cout << "  k-mers:           " << sk.get_nkmers() << "\n";
    std::cout << "  card:              " << sk.get_card() << "\n";
    std::cout << "  rho:              " << sk.get_rho() << "\n";
    std::cout << "  Nonempty buckets: " << sk.get_buckets().get_nnonempty() << "\n";
    std::cout << "  Windows:          " << sk.get_windows().wins_v.size() << "\n";
    const sketch_entry_t& e = file.get_entry(i);
    std::cout << "  Bucket bytes:     " << e.buckets_len << "\n";
    std::cout << "  Window bytes:     " << e.windows_len << "\n";
  }
}

InfoSC::InfoSC(CLI::App& sc)
{
  sc.add_option("-i,--sketch-path", sketch_path, "Container to inspect")->required()->check(CLI::ExistingFile);
}

int main(int argc, char** argv)
{
  PRINT_VERSION
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  CLI::App app{"gdiff"};
  app.set_help_flag("--help");
  app.fallthrough();

  app.add_flag("--verbose,!--no-verbose", verbose, "Report progress even when stderr is not a terminal");
  app.require_subcommand();
  app.add_option("--seed", seed, "Random seed for the LSH and other parts that require randomness [0]");
  app.callback([&]() { init_thread_rng(0); });
  app.add_option("--num-threads", num_threads, "Number of threads for parallel sketch/map/dist processing [1]");

  auto& sc_sketch = *app.add_subcommand("sketch", "Create sketches from FASTA/FASTQ files");
  auto& sc_map = *app.add_subcommand("map", "Map queries and extract distance-based patterns from sketches");
  auto& sc_dist = *app.add_subcommand("dist", "Summarize MLE distances between sketches");
  auto& sc_detect = *app.add_subcommand("detect", "Fit a background distance distribution and detect outlier regions");
  auto& sc_merge = *app.add_subcommand("merge", "Merge containers into one");
  auto& sc_info = *app.add_subcommand("info", "Show metadata for all sketches in a container");
  auto& sc_roll = *app.add_subcommand("roll", "Roll a window over query sequences and report per-window MLE distances");

  SketchSC gdiff_sketch(sc_sketch);
  MapSC gdiff_map(sc_map);
  DistSC gdiff_dist(sc_dist);
  DetectSC gdiff_detect(sc_detect);
  MergeSC gdiff_merge(sc_merge);
  InfoSC gdiff_info(sc_info);
  RollSC gdiff_roll(sc_roll);

  CLI11_PARSE(app, argc, argv);
  str invocation;
  for (int i = 0; i < argc; ++i) {
    invocation += str(argv[i]) + " ";
  }
  if (!invocation.empty()) {
    invocation.pop_back();
  }

  const auto tstart = std::chrono::system_clock::now();
  const std::time_t tstart_f = std::chrono::system_clock::to_time_t(tstart);
  std::cerr << "Invocation: " << invocation << '\n' << std::ctime(&tstart_f);

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
    cerr_msg("Comparing sketches and calculating distances...");
    run_timed("Done calculating distances, elapsed: ", [&]() { gdiff_dist.dist(); });
  }

  if (sc_detect.parsed()) {
    cerr_msg("Sampling distances and detecting outlier regions...");
    run_timed("Done detecting outlier regions, elapsed: ", [&]() { gdiff_detect.detect(); });
  }

  if (sc_info.parsed()) {
    gdiff_info.info();
  }

  if (sc_roll.parsed()) {
    cerr_msg("Rolling windows and calculating local distances...");
    run_timed("Done rolling windows, elapsed: ", [&]() { gdiff_roll.roll(); });
  }

  auto tend = std::chrono::system_clock::now();
  std::time_t tend_f = std::chrono::system_clock::to_time_t(tend);
  std::cerr << std::ctime(&tend_f);

  return 0;
}
