// Temporary sketch2/dist2 driver. Safe to delete with sketch2.* / dist2.* / makefile gdiff2 block.
#include "CLI11.hpp"
#include "dist2.hpp"
#include "msg.hpp"
#include "random.hpp"
#include "sketch2.hpp"

#include <chrono>
#include <ctime>
#include <iostream>

#define GDIFF2_VERSION "v0.0.1-sketch2"
#define PRINT_GDIFF2_VERSION std::cerr << "gdiff2 version: " << GDIFF2_VERSION << std::endl;

extern uint32_t num_threads;
extern str invocation;

int main(int argc, char** argv)
{
  PRINT_GDIFF2_VERSION
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  CLI::App app{"gdiff2 (sketch2/dist2 experiment)"};
  app.set_help_flag("--help");
  app.fallthrough();

  app.require_subcommand();
  app.add_option("--seed", seed, "Random seed for the LSH and sampling [0]");
  app.callback([&]() { init_thread_rng(0); });
  app.add_option("--num-threads", num_threads, "Number of worker threads for sketch2/dist2 [1]");

  auto& sc_sketch2 = *app.add_subcommand("sketch2", "Create sketches with sampled window hash values");
  auto& sc_dist2 = *app.add_subcommand("dist2", "Distance between two sketch2 files using buckets and window hashes");

  Sketch2SC gdiff_sketch2(sc_sketch2);
  Dist2SC gdiff_dist2(sc_dist2);

  CLI11_PARSE(app, argc, argv);
  invocation.clear();
  for (int i = 0; i < argc; ++i) {
    if (i) invocation += ' ';
    invocation += str(argv[i]);
  }

  auto tstart = std::chrono::system_clock::now();
  std::time_t tstart_f = std::chrono::system_clock::to_time_t(tstart);
  std::cerr << "Invocation: " << invocation << '\n';
  std::cerr << std::ctime(&tstart_f);

  auto run_timed = [&](const char* done_msg, auto&& work) {
    const auto t0 = std::chrono::system_clock::now();
    work();
    std::chrono::duration<float> es = std::chrono::system_clock::now() - t0;
    cerr_msg(done_msg, es.count(), " sec");
  };

  if (sc_sketch2.parsed()) {
    cerr_msg("Initializing the sketch2...");
    gdiff_sketch2.set_nrows();
    gdiff_sketch2.set_lshf();
    run_timed("Done sketch2 & saving, elapsed: ", [&]() { gdiff_sketch2.process(); });
  }

  if (sc_dist2.parsed()) {
    cerr_msg("Comparing sketch2 files...");
    run_timed("Done dist2, elapsed: ", [&]() { gdiff_dist2.dist(); });
  }

  auto tend = std::chrono::system_clock::now();
  std::time_t tend_f = std::chrono::system_clock::to_time_t(tend);
  std::cerr << std::ctime(&tend_f);

  return 0;
}
