#include "info.hpp"

#include "utils.hpp"
#include "container.hpp"
#include "msg.hpp"
#include "sketch.hpp"

#include <ctime>
#include <iostream>

namespace {

  void print_sketch_info(const Container& file, uint32_t i)
  {
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
    if (!sk.get_spurious_v().empty()) {
      uint64_t sp_total = 0;
      for (uint64_t c : sk.get_spurious_v())
        sp_total += c;
      std::cout << "  spurious:         " << sp_total << " reference k-mers\n";
    }
    std::cout << "  Nonempty buckets: " << sk.get_buckets().get_nnonempty() << "\n";
    std::cout << "  Windows:          " << sk.get_windows().wins_v.size() << "\n";
    const scentry& e = file.get_entry(i);
    std::cout << "  Bucket bytes:     " << e.buckets_len << "\n";
    std::cout << "  Window bytes:     " << e.windows_len << "\n";
  }

} // namespace

void InfoSC::info()
{
  if (!check_input_type(sketch_path, input_type::container, "info -i/--sketch-path")) exit_failure();

  write_provenance(std::cout);

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

  for (uint32_t i = 0; i < file.size(); ++i)
    print_sketch_info(file, i);
}

InfoSC::InfoSC(CLI::App& sc)
{
  sc.add_option("-i,--sketch-path", sketch_path, "Container to inspect")->required()->check(CLI::ExistingFile);
}
