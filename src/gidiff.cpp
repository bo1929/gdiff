#include "gidiff.hpp"

void BaseLSH::set_lshf() { lshf = std::make_shared<LSHF>(k, h, m); }

void BaseLSH::set_nrows()
{
  uint32_t hash_size = pow(2, 2 * h);
  uint32_t full_residue = hash_size % m;
  if (frac) {
    nrows = (hash_size / m) * (r + 1);
    nrows = full_residue > r ? nrows + (r + 1) : nrows + full_residue;
  } else {
    nrows = (hash_size / m);
    nrows = full_residue > r ? nrows + 1 : nrows;
  }
}

bool BaseLSH::validate_configuration()
{
  bool is_invalid = true;
  if ((is_invalid = w < k)) {
    std::cerr << "The minimum minimizer window size (-w) is k (-k)!" << std::endl;
  }
  if ((is_invalid = h < 3)) {
    std::cerr << "The minimum number of LSH positions (-h) is 3!" << std::endl;
  }
  if ((is_invalid = h > 15)) {
    std::cerr << "The maximum number of LSH positions (-h) is 15!" << std::endl;
  }
  if ((is_invalid = k > 31)) {
    std::cerr << "The maximum allowed k-mer length (-k) is 31!" << std::endl;
  }
  if ((is_invalid = k < 19)) {
    std::cerr << "The minimum allowed k-mer length (-k) is 19!" << std::endl;
  }
  if ((is_invalid = (k - h) > 16)) {
    std::cerr << "For compact k-mer encodings, h must be >= k-16!" << std::endl;
  }
  // if (sdust_t == 0 || sdust_w == 0) {
  //   std::cerr << "Setting --sdust-w or --sdust-t to 0 will disable dustmasker!" << std::endl;
  // }
  return !is_invalid;
}

void BaseLSH::save_configuration(std::ofstream& cfg_stream)
{
  cfg_stream.write(reinterpret_cast<char*>(&k), sizeof(uint8_t));
  cfg_stream.write(reinterpret_cast<char*>(&w), sizeof(uint8_t));
  cfg_stream.write(reinterpret_cast<char*>(&h), sizeof(uint8_t));
  cfg_stream.write(reinterpret_cast<char*>(&m), sizeof(uint32_t));
  cfg_stream.write(reinterpret_cast<char*>(&r), sizeof(uint32_t));
  cfg_stream.write(reinterpret_cast<char*>(&frac), sizeof(bool));
  cfg_stream.write(reinterpret_cast<char*>(&nrows), sizeof(uint32_t));
  cfg_stream.write(reinterpret_cast<char*>(lshf->ppos_data()), (h) * sizeof(uint8_t));
  cfg_stream.write(reinterpret_cast<char*>(lshf->npos_data()), (k - h) * sizeof(uint8_t));
}

void MapSketch::load_sketch()
{
  sketch->load_full_sketch();
  sketch->make_rho_partial();
}

void SketchTarget::create_sketch()
{
  rseq_sptr_t rs = std::make_shared<RSeq>(input, lshf, w, r, frac);
  sdhm_sptr_t sdhm = std::make_shared<SDHM>();
  sdhm->fill_table(nrows, rs);
  std::cout << sdhm->get_nkmers() << std::endl;
  sketch_sfhm = std::make_shared<SFHM>(sdhm);

  rho = rs->get_rho();
  std::cerr << "Total number of k-mers included in the sketch: " << sdhm->get_nkmers() << std::endl;
  std::cerr << "Subsampling rate (rho) is: " << rho << std::endl;
}

void SketchTarget::save_sketch()
{
  std::ofstream sketch_stream(sketch_path, std::ofstream::binary);
  sketch_sfhm->save(sketch_stream);
  save_configuration(sketch_stream);
  sketch_stream.write(reinterpret_cast<char*>(&rho), sizeof(double));
  check_fstream(sketch_stream, "Failed to write the sketch!", sketch_path);
  sketch_stream.close();
}

void MapSketch::header_dreport(strstream& dreport_stream) {}

void MapSketch::map_sequences()
{
  strstream dreport_stream;
  header_dreport(dreport_stream);
  qseq_sptr_t qs = std::make_shared<QSeq>(query);
  bool cont_reading = false;
  while ((cont_reading = qs->read_next_batch()) || !qs->is_batch_finished()) {
    total_qseq += qs->get_cbatch_size();
    SBatch sb(sketch, qs, hdist_th, dist_th, min_length, chisq);
    {
      sb.map_sequences(*output_stream);
    }
  }
}

SketchTarget::SketchTarget(CLI::App& sc)
{
  set_sketch_defaults();
  sc.add_option("-i,--input-file", input, "Input FASTA/FASTQ file <path> (or URL) (gzip compatible).")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("-o,--output-path", sketch_path, "Path to store the resulting binary sketch file.")->required();
  sc.add_option("-k,--kmer-len", k, "Length of k-mers. [26]")->check(CLI::Range(19, 31));
  sc.add_option("-w,--win-len", w, "Length of minimizer window (w>=k). [k+6]");
  sc.add_option("-h,--num-positions", h, "Number of positions for the LSH. [k-16]");
  sc.add_option("-m,--modulo-lsh", m, "Modulo value to partition LSH space. [4]")->check(CLI::PositiveNumber);
  sc.add_option("-r,--residue-lsh", r, "A k-mer x will be included only if r = LSH(x) mod m. [1]")
    ->check(CLI::NonNegativeNumber);
  sc.add_flag("--frac,!--no-frac", frac, "Include k-mers with r <= LSH(x) mod m. [true]");
  // sc.add_option("--sdust-t", sdust_t, "SDUST threshold (NCBI dustmasker: 20). [0]")->check(CLI::NonNegativeNumber);
  // sc.add_option("--sdust-w", sdust_w, "SDUST window (NCBI dustmasker: 64). [0]")->check(CLI::NonNegativeNumber);
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

MapSketch::MapSketch(CLI::App& sc)
{
  sc.add_option("-q,--query", query, "Query FASTA/FASTQ file <path> (or URL) (gzip compatible).")
    ->required()
    ->check(url_validator | CLI::ExistingFile);
  sc.add_option("-i,--sketch-path", sketch_path, "Sketch file at <path> to query.")->required()->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", output_path, "Write output to a file at <path>. [stdout]");
  sc.add_option("--hdist-th", hdist_th, "Maximum Hamming distance for a k-mer to match. [4]")->check(CLI::NonNegativeNumber);
  sc.add_option("--chi-sq", chisq, "Maximum Hamming distance for a k-mer to match. [4]")->check(CLI::NonNegativeNumber);
  sc.add_option("-d,--dist-th", dist_th, "Maximum (or minimum) distance for an interval to match.")->required();
  sc.add_option("-l,--min-length", min_length, "Maximum (or minimum) length for an interval to match.")->required();
  sc.callback([&]() {
    if (!output_path.empty()) {
      output_file.open(output_path);
      output_stream = &output_file;
    }
    sketch = std::make_shared<Sketch>(sketch_path);
  });
}

int main(int argc, char** argv)
{
  PRINT_VERSION
  std::ios::sync_with_stdio(false);
  CLI::App app{"gidiff"};
  app.set_help_flag("--help");
  app.fallthrough();

  bool verbose = false;
  app.add_flag("--verbose,!--no-verbose", verbose, "Increased verbosity and progress report.");
  app.require_subcommand();
  app.add_option("--seed", seed, "Random seed for the LSH and other parts that require randomness. [0]");
  app.callback([&]() {
    if (app.count("--seed")) {
      gen.seed(seed);
    }
  });
  app.add_option("--num-threads", num_threads, "Number of threads to use in OpenMP-based parallelism. [1]");

  auto& sc_sketch = *app.add_subcommand("sketch", "Create a sketch from k-mers in a single FASTA/FASTQ file.");
  auto& sc_map = *app.add_subcommand("map", "Seek query sequences in a sketch and estimate distances.");

  SketchTarget krepp_sketch(sc_sketch);
  MapSketch krepp_map(sc_map);

  CLI11_PARSE(app, argc, argv);
  for (int i = 0; i < argc; ++i) {
    invocation += std::string(argv[i]) + " ";
  }
  invocation.pop_back();

  auto tstart = std::chrono::system_clock::now();
  std::time_t tstart_f = std::chrono::system_clock::to_time_t(tstart);
  std::cerr << "Invocation: " << invocation << "\n";
  std::cerr << std::ctime(&tstart_f);

  if (sc_sketch.parsed()) {
    std::cerr << "Initializing the sketch..." << std::endl;
    krepp_sketch.set_nrows();
    krepp_sketch.set_lshf();
    std::chrono::duration<float> es_b = std::chrono::system_clock::now() - tstart;
    krepp_sketch.create_sketch();
    krepp_sketch.save_sketch();
    std::chrono::duration<float> es_s = std::chrono::system_clock::now() - tstart - es_b;
    std::cerr << "Done sketching & saving, elapsed: " << es_s.count() << " sec" << std::endl;
  }

  if (sc_map.parsed()) {
    std::cerr << "Loading the sketch..." << std::endl;
    krepp_map.load_sketch();
    std::cerr << "Seeking query sequences in the sketch..." << std::endl;
    std::chrono::duration<float> es_b = std::chrono::system_clock::now() - tstart;
    krepp_map.map_sequences();
    std::chrono::duration<float> es_s = std::chrono::system_clock::now() - tstart - es_b;
    std::cerr << "Done maping sequences, elapsed: " << es_s.count() << " sec" << std::endl;
    std::cerr << "Total number of sequences queried: " << krepp_map.get_total_qseq() << std::endl;
  }

  auto tend = std::chrono::system_clock::now();
  std::time_t tend_f = std::chrono::system_clock::to_time_t(tend);
  std::cerr << std::ctime(&tend_f);

  return 0;
}
