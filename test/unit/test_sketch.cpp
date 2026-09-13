// GDSK write/read round-trip for buckets and both window representations.
#include "doctest/doctest.h"
#include "hm.hpp"
#include "sketch.hpp"
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <sstream>

namespace {

  const std::string GDIFF_BIN = "./gdiff";

  bool gdiff_available() { return std::filesystem::exists(GDIFF_BIN); }

  std::filesystem::path tmp_path(const std::string& name)
  {
    return std::filesystem::temp_directory_path() / name;
  }

  // Random sequence: tests assert container structure, not which k-mers survive.
  std::filesystem::path write_fasta(const std::string& name, uint64_t len, uint32_t seed)
  {
    const auto path = tmp_path(name + ".fa");
    std::mt19937_64 rng(seed);
    const char bases[] = "ACGT";
    std::ofstream out(path);
    out << ">" << name << "\n";
    for (uint64_t i = 0; i < len; ++i) {
      out << bases[rng() & 3];
      if ((i + 1) % 80 == 0) out << "\n";
    }
    out << "\n";
    return path;
  }

  // Returns an empty path when the sketch call fails, so callers can REQUIRE.
  std::filesystem::path run_sketch(const std::string& out_name, const std::string& args)
  {
    const auto path = tmp_path(out_name + ".gdsk");
    const std::string cmd = GDIFF_BIN + " sketch " + args + " -o " + path.string() + " >/dev/null 2>&1";
    if (std::system(cmd.c_str()) != 0) return {};
    return path;
  }

} // namespace

TEST_SUITE("GDSK container")
{

  TEST_CASE("record fields survive a write/read round trip" * doctest::skip(!gdiff_available()))
  {
    const auto fa = write_fasta("gdsk_one", 200000, 1);
    const auto skc = run_sketch("gdsk_one", "-i " + fa.string() + " --frac 0.5");
    REQUIRE_FALSE(skc.empty());

    const SketchFile sf(skc);
    REQUIRE(sf.size() == 1);

    const sketch_config_t& cfg = sf.get_config();
    CHECK(cfg.k == 27);
    CHECK(cfg.h == 11);
    CHECK(cfg.canonical);
    CHECK(cfg.frac == doctest::Approx(0.5));
    CHECK(cfg.ppos.size() == cfg.h);
    CHECK(cfg.npos.size() == size_t(cfg.k - cfg.h));

    const Sketch sk = sf.open(0);
    CHECK(sk.get_rname() == fa.filename().string());
    CHECK(sk.get_nkmers() > 0);
    CHECK(sk.get_card_est() > 0.0);
    // Random sequence has almost no duplicate k-mers, so rho is in (0, 1].
    CHECK(sk.get_rho() > 0.0);
    CHECK(sk.get_rho() <= 1.0);
    CHECK(sk.get_timestamp() > 0);
    CHECK(sk.get_genome_bp() == 200000);
    CHECK(sk.has_buckets());
    CHECK(sk.has_windows());
    CHECK(sk.get_buckets().get_nkmers() == sk.get_nkmers());

    std::filesystem::remove(fa);
    std::filesystem::remove(skc);
  }

  TEST_CASE("every record of a multi-record container is independently addressable" *
            doctest::skip(!gdiff_available()))
  {
    const auto fa1 = write_fasta("gdsk_multi1", 120000, 2);
    const auto fa2 = write_fasta("gdsk_multi2", 90000, 3);
    const auto fa3 = write_fasta("gdsk_multi3", 60000, 4);
    const auto skc =
      run_sketch("gdsk_multi", "-i " + fa1.string() + " -i " + fa2.string() + " -i " + fa3.string());
    REQUIRE_FALSE(skc.empty());

    const SketchFile sf(skc);
    REQUIRE(sf.size() == 3);

    // Records are indexed in input order, each at its own offset.
    const std::vector<std::string> expect{
      fa1.filename().string(), fa2.filename().string(), fa3.filename().string()};
    const std::vector<uint64_t> expect_bp{120000, 90000, 60000};
    std::set<uint64_t> offsets;
    for (uint32_t r = 0; r < sf.size(); ++r) {
      const Sketch sk = sf.open(r);
      CHECK(sk.get_rname() == expect[r]);
      CHECK(sk.get_genome_bp() == expect_bp[r]);
      CHECK(sk.get_buckets().get_nkmers() > 0);
      offsets.insert(sf.get_index().records[r].offset);
    }
    CHECK(offsets.size() == sf.size());

    std::filesystem::remove(fa1);
    std::filesystem::remove(fa2);
    std::filesystem::remove(fa3);
    std::filesystem::remove(skc);
  }

  TEST_CASE("SketchPart loads only the requested sections" * doctest::skip(!gdiff_available()))
  {
    const auto fa = write_fasta("gdsk_part", 150000, 5);
    const auto skc = run_sketch("gdsk_part", "-i " + fa.string());
    REQUIRE_FALSE(skc.empty());

    const SketchFile sf(skc);
    const Sketch buckets_only = sf.open(0, SketchPart::Buckets);
    CHECK(buckets_only.has_buckets());
    CHECK_FALSE(buckets_only.has_windows());

    const Sketch windows_only = sf.open(0, SketchPart::Windows);
    CHECK_FALSE(windows_only.has_buckets());
    CHECK(windows_only.has_windows());

    // Metadata lives in the record preamble, so it is present either way.
    CHECK(buckets_only.get_rname() == windows_only.get_rname());
    CHECK(buckets_only.get_rho() == windows_only.get_rho());

    std::filesystem::remove(fa);
    std::filesystem::remove(skc);
  }

  TEST_CASE("pool and seq records agree on window geometry" * doctest::skip(!gdiff_available()))
  {
    const auto fa = write_fasta("gdsk_repr", 200000, 6);
    const auto pool = run_sketch("gdsk_pool", "-i " + fa.string() + " --window-repr pool");
    const auto seq = run_sketch("gdsk_seq", "-i " + fa.string() + " --window-repr seq");
    REQUIRE_FALSE(pool.empty());
    REQUIRE_FALSE(seq.empty());

    const SketchFile pf(pool), qf(seq);
    CHECK(pf.get_config().win_repr == WinRepr::Pool);
    CHECK(qf.get_config().win_repr == WinRepr::Seq);
    // win_repr is excluded from compatibility: Pool and Seq of one geometry compare.
    CHECK(pf.get_config().compatible_with(qf.get_config()));

    const Sketch ps = pf.open(0), qs = qf.open(0);
    const vec<window_t>& pw = ps.get_wins();
    const vec<window_t>& qw = qs.get_wins();
    REQUIRE(pw.size() == qw.size());
    REQUIRE_FALSE(pw.empty());
    // Same seed and genome: same sampled coordinates.
    for (size_t i = 0; i < pw.size(); ++i) {
      CHECK(pw[i].qid == qw[i].qid);
      CHECK(pw[i].start == qw[i].start);
      CHECK(pw[i].end == qw[i].end);
    }

    std::filesystem::remove(fa);
    std::filesystem::remove(pool);
    std::filesystem::remove(seq);
  }

  TEST_CASE("seq windows unpack to the sampled length" * doctest::skip(!gdiff_available()))
  {
    const auto fa = write_fasta("gdsk_unpack", 150000, 7);
    const auto skc = run_sketch("gdsk_unpack", "-i " + fa.string() + " --window-repr seq -l 500");
    REQUIRE_FALSE(skc.empty());

    const SketchFile sf(skc);
    const Sketch sk = sf.open(0, SketchPart::Windows);
    const window_sample_t& ws = sk.get_windows();
    REQUIRE_FALSE(ws.wins.empty());

    str cseq;
    for (uint64_t wi = 0; wi < std::min<uint64_t>(ws.wins.size(), 16); ++wi) {
      ws.packs.unpack(wi, cseq);
      // A window of n k-mers spans n + k - 1 bases.
      CHECK(cseq.size() == ws.wins[wi].end - ws.wins[wi].start + sk.get_k() - 1);
      CHECK(cseq.find_first_not_of("ACGTN") == str::npos);
    }

    std::filesystem::remove(fa);
    std::filesystem::remove(skc);
  }

  TEST_CASE("canonicalize is a no-op on a canonical sketch" * doctest::skip(!gdiff_available()))
  {
    const auto fa = write_fasta("gdsk_canon", 120000, 8);
    const auto skc = run_sketch("gdsk_canon", "-i " + fa.string());
    REQUIRE_FALSE(skc.empty());

    const SketchFile sf(skc);
    Sketch sk = sf.open(0, SketchPart::Buckets);
    REQUIRE(sk.is_canonical());

    const uint64_t before = sk.get_buckets().get_nkmers();
    sk.canonicalize();
    CHECK(sk.is_canonical());
    CHECK(sk.get_buckets().get_nkmers() == before);

    std::filesystem::remove(fa);
    std::filesystem::remove(skc);
  }

  TEST_CASE("strand-aware containers are flagged as such" * doctest::skip(!gdiff_available()))
  {
    const auto fa = write_fasta("gdsk_sa", 120000, 9);
    const auto ag = run_sketch("gdsk_ag", "-i " + fa.string());
    const auto sa = run_sketch("gdsk_sa", "-i " + fa.string() + " --strand-aware");
    REQUIRE_FALSE(ag.empty());
    REQUIRE_FALSE(sa.empty());

    CHECK(SketchFile(ag).open(0).is_canonical());
    CHECK_FALSE(SketchFile(sa).open(0).is_canonical());
    // Strandedness changes what a bucket means, so the two must not compare.
    CHECK_FALSE(SketchFile(ag).get_config().compatible_with(SketchFile(sa).get_config()));

    std::filesystem::remove(fa);
    std::filesystem::remove(ag);
    std::filesystem::remove(sa);
  }

  TEST_CASE("is_sketch_file accepts containers and rejects anything else" *
            doctest::skip(!gdiff_available()))
  {
    const auto fa = write_fasta("gdsk_probe", 60000, 10);
    const auto skc = run_sketch("gdsk_probe", "-i " + fa.string());
    REQUIRE_FALSE(skc.empty());

    CHECK(is_sketch_file(skc));
    CHECK_FALSE(is_sketch_file(fa));
    CHECK_FALSE(is_sketch_file(tmp_path("gdsk_does_not_exist.gdsk")));

    const auto empty = tmp_path("gdsk_empty.bin");
    std::ofstream(empty).close();
    CHECK_FALSE(is_sketch_file(empty));

    std::filesystem::remove(fa);
    std::filesystem::remove(skc);
    std::filesystem::remove(empty);
  }

  TEST_CASE("a corrupt magic is rejected rather than misread" * doctest::skip(!gdiff_available()))
  {
    const auto fa = write_fasta("gdsk_bad", 60000, 11);
    const auto skc = run_sketch("gdsk_bad", "-i " + fa.string());
    REQUIRE_FALSE(skc.empty());

    const auto bad = tmp_path("gdsk_bad_magic.gdsk");
    std::filesystem::copy_file(skc, bad, std::filesystem::copy_options::overwrite_existing);
    {
      std::fstream f(bad, std::ios::binary | std::ios::in | std::ios::out);
      const uint32_t junk = 0xdeadbeefu;
      f.write(reinterpret_cast<const char*>(&junk), sizeof(junk));
    }
    CHECK_FALSE(is_sketch_file(bad));
    // Opening a bad file calls error_exit; observe it out of process.
    const std::string cmd = GDIFF_BIN + " info -i " + bad.string() + " >/dev/null 2>&1";
    CHECK(std::system(cmd.c_str()) != 0);

    std::filesystem::remove(fa);
    std::filesystem::remove(skc);
    std::filesystem::remove(bad);
  }

} // TEST_SUITE

TEST_SUITE("Buckets")
{

  TEST_CASE("range resolves owned entries and reports empty buckets")
  {
    const uint32_t nrows = 256;
    vec<uint64_t> keys{
      pack_key(0, 100), pack_key(0, 50), pack_key(3, 200), pack_key(255, 7), pack_key(3, 200)};
    Buckets b;
    b.build(nrows, std::move(keys));

    // build() uniques, so the duplicate in bucket 3 collapses.
    CHECK(b.get_nkmers() == 4);
    CHECK(b.get_nnonempty() == 3);
    CHECK(b.get_nrows() == nrows);
    CHECK_FALSE(b.empty());

    const enc_t *beg = nullptr, *end = nullptr;
    REQUIRE(b.range(0, beg, end));
    REQUIRE(end - beg == 2);
    // Entries are sorted within a bucket because the packed keys are.
    CHECK(beg[0] == 50);
    CHECK(beg[1] == 100);

    REQUIRE(b.range(3, beg, end));
    CHECK(end - beg == 1);
    CHECK(beg[0] == 200);

    REQUIRE(b.range(255, beg, end));
    CHECK(end - beg == 1);
    CHECK(beg[0] == 7);

    CHECK_FALSE(b.range(1, beg, end));
    CHECK_FALSE(b.range(254, beg, end));
    CHECK_FALSE(b.range(nrows, beg, end));      // out of range
    CHECK_FALSE(b.range(nrows + 1000, beg, end));

    CHECK(b.nonempty(0));
    CHECK_FALSE(b.nonempty(1));
    CHECK_FALSE(b.nonempty(nrows));
  }

  TEST_CASE("a saved table views back identically")
  {
    const uint32_t nrows = 1024;
    std::mt19937_64 rng(42);
    vec<uint64_t> keys;
    for (int i = 0; i < 500; ++i)
      keys.push_back(pack_key(static_cast<uint32_t>(rng() % nrows), static_cast<enc_t>(rng())));

    Buckets owned;
    owned.build(nrows, vec<uint64_t>(keys));

    std::ostringstream os(std::ios::binary);
    owned.save(os);
    const std::string bytes = os.str();
    CHECK(bytes.size() == Buckets::byte_size(owned.get_nkmers(), owned.get_nnonempty(), nrows));

    Buckets viewed;
    const char* p = bytes.data();
    viewed.view(p, bytes.data() + bytes.size(), nrows);
    CHECK(p == bytes.data() + bytes.size()); // view consumed exactly the payload

    REQUIRE(viewed.get_nkmers() == owned.get_nkmers());
    REQUIRE(viewed.get_nnonempty() == owned.get_nnonempty());
    for (uint32_t rix = 0; rix < nrows; ++rix) {
      const enc_t *a1 = nullptr, *a2 = nullptr, *b1 = nullptr, *b2 = nullptr;
      const bool ha = owned.range(rix, a1, a2);
      const bool hb = viewed.range(rix, b1, b2);
      REQUIRE(ha == hb);
      if (!ha) continue;
      REQUIRE((a2 - a1) == (b2 - b1));
      for (; a1 < a2; ++a1, ++b1) CHECK(*a1 == *b1);
    }

    // skip() must land where view() did.
    CHECK(Buckets::skip(bytes.data(), bytes.data() + bytes.size()) ==
          bytes.data() + bytes.size());
  }

  TEST_CASE("an empty table answers every probe with false")
  {
    Buckets b;
    CHECK(b.empty());
    b.build(64, {});
    CHECK(b.get_nkmers() == 0);
    CHECK(b.get_nnonempty() == 0);
    const enc_t *beg = nullptr, *end = nullptr;
    for (uint32_t rix = 0; rix < 64; ++rix)
      CHECK_FALSE(b.range(rix, beg, end));
  }

  TEST_CASE("move leaves the source empty and the target intact")
  {
    Buckets src;
    src.build(128, {pack_key(5, 11), pack_key(5, 12), pack_key(9, 13)});
    const uint64_t nkmers = src.get_nkmers();

    Buckets dst(std::move(src));
    CHECK(dst.get_nkmers() == nkmers);
    CHECK(src.empty());

    const enc_t *beg = nullptr, *end = nullptr;
    REQUIRE(dst.range(5, beg, end));
    CHECK(end - beg == 2);
  }

} // TEST_SUITE
