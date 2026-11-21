#include "utils.hpp"
#include "dist.hpp"
#include "info.hpp"
#include "map.hpp"
#include "merge.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "roll.hpp"
#include "sketch.hpp"

#include <chrono>
#include <cmath>
#include <iostream>

static int run(int argc, char** argv)
{
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  CLI::App app{"gdiff"};
  app.set_help_flag("--help");
  app.fallthrough();

  app.add_flag("--verbose,!--no-verbose", verbose, "Report progress even when stderr is not a terminal and details");
  app.require_subcommand(1);
  app.add_option("--seed", seed, "Random seed for the LSH and other parts that require randomness [0]");
  app.callback([&]() { init_thread_rng(0); });
  app.add_option("-t,--num-threads", num_threads, "Worker threads per subcommand [1]");

  auto& sc_sketch = *app.add_subcommand("sketch", "Create sketches from FASTA/FASTQ files");
  auto& sc_map = *app.add_subcommand("map", "Map queries and extract distance-based patterns from sketches");
  auto& sc_dist = *app.add_subcommand("dist", "Summarize MLE distances between sketches");
  auto& sc_merge = *app.add_subcommand("merge", "Merge containers into one");
  auto& sc_info = *app.add_subcommand("info", "Show metadata for all sketches in a container");
  auto& sc_roll = *app.add_subcommand("roll", "Roll a window over query sequences and report per-window MLE distances");

  SketchSC gdiff_sketch(sc_sketch);
  MapSC gdiff_map(sc_map);
  DistSC gdiff_dist(sc_dist);
  MergeSC gdiff_merge(sc_merge);
  InfoSC gdiff_info(sc_info);
  RollSC gdiff_roll(sc_roll);

  CLI11_PARSE(app, argc, argv);
  invocation.clear();
  for (int i = 0; i < argc; ++i) {
    invocation += str(argv[i]) + " ";
  }
  if (!invocation.empty()) {
    invocation.pop_back();
  }

  auto get_running_time = [&](auto&& work) {
    const auto t0 = std::chrono::steady_clock::now();
    work();
    const std::chrono::duration<double> es = std::chrono::steady_clock::now() - t0;
    info_msg("Done in ", std::round(es.count() * 100.0) / 100.0, " s");
  };

  if (sc_sketch.parsed()) {
    write_provenance(std::cerr);
    gdiff_sketch.set_nrows();
    gdiff_sketch.set_lshf();
    get_running_time([&]() { gdiff_sketch.process(); });
  }

  if (sc_merge.parsed()) {
    write_provenance(std::cerr);
    get_running_time([&]() { gdiff_merge.merge(); });
  }

  if (sc_map.parsed()) get_running_time([&]() { gdiff_map.map(); });
  if (sc_dist.parsed()) get_running_time([&]() { gdiff_dist.dist(); });
  if (sc_info.parsed()) gdiff_info.info();
  if (sc_roll.parsed()) get_running_time([&]() { gdiff_roll.roll(); });

  return 0;
}

int main(int argc, char** argv)
{
  try {
    return run(argc, argv);
  } catch (const fatal_error& e) {
    error_msg(e.what());
  } catch (const std::exception& e) {
    error_msg("Unexpected failure: ", e.what());
  }
  return EXIT_FAILURE;
}
