// Container write/read round-trip for buckets and both window representations.
#include "doctest/doctest.h"
#include "buckets.hpp"
#include "sketch.hpp"
#include "test_util.hpp"
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <set>
#include <sstream>

namespace {

  const std::string GDIFF_BIN = "./gdiff";

  bool gdiff_available() { return std::filesystem::exists(GDIFF_BIN); }

  // Returns an empty path when the sketch call fails, so callers can REQUIRE.
  std::filesystem::path run_sketch(const std::string& out_name, const std::string& args)
  {
    const auto path = test_util::path(out_name + ".gs");
    const std::string cmd = GDIFF_BIN + " sketch " + args + " -o " + path.string() + " >/dev/null 2>&1";
    if (std::system(cmd.c_str()) != 0) return {};
    return path;
  }

  // The i-th (0-based) tab-separated field of a TSV line; empty when the line is too short.
  std::string field(const std::string& line, size_t i)
  {
    size_t begin = 0;
    for (size_t f = 0; f < i; ++f) {
      begin = line.find('\t', begin);
      if (begin == std::string::npos) return {};
      ++begin;
    }
    const size_t end = line.find('\t', begin);
    return line.substr(begin, end == std::string::npos ? std::string::npos : end - begin);
  }

  // Returns the dist command's stdout, or an empty string when it fails.
  std::string run_dist(const std::string& args)
  {
    const auto path = test_util::path("dist_cli.tsv");
    const std::string cmd = GDIFF_BIN + " dist " + args + " > " + path.string() + " 2>/dev/null";
    if (std::system(cmd.c_str()) != 0) return {};
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    std::filesystem::remove(path);
    return ss.str();
  }

  // Returns the merged container's path, or an empty path when merge fails.
  std::filesystem::path run_merge(const std::string& out_name, const std::vector<std::string>& inputs)
  {
    const auto path = test_util::path(out_name + ".gs");
    std::string cmd = GDIFF_BIN + " merge";
    for (const std::string& in : inputs)
      cmd += " -i " + in;
    cmd += " -o " + path.string() + " >/dev/null 2>&1";
    if (std::system(cmd.c_str()) != 0) return {};
    return path;
  }

  // Returns the info report's stdout, or an empty string when info fails.
  std::string run_info(const std::filesystem::path& container)
  {
    const auto path = test_util::path("info_report.txt");
    const std::string cmd = GDIFF_BIN + " info -i " + container.string() + " > " + path.string() + " 2>/dev/null";
    if (std::system(cmd.c_str()) != 0) return {};
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
  }

  // The value on the line labelled `label`, alignment padding stripped.
  std::string info_value(const std::string& report, const std::string& label)
  {
    const size_t at = report.find(label);
    if (at == std::string::npos) return {};
    size_t b = at + label.size();
    while (b < report.size() && report[b] == ' ') ++b;
    const size_t e = report.find('\n', b);
    return report.substr(b, e == std::string::npos ? std::string::npos : e - b);
  }

  size_t count_of(const std::string& hay, const std::string& needle)
  {
    size_t n = 0;
    for (size_t at = hay.find(needle); at != std::string::npos; at = hay.find(needle, at + needle.size())) ++n;
    return n;
  }

  // What `--seed` controls in a sketch: the sampled window starts, in order.
  vec<uint64_t> window_starts(const Sketch& sk)
  {
    vec<uint64_t> v;
    for (const window_t& w : sk.get_windows().wins_v)
      v.push_back(w.start);
    return v;
  }

} // namespace

TEST_SUITE("container")
{

  TEST_CASE("record fields survive a write/read round trip" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("cont_one", 200000, 1);
    const auto skc = run_sketch("cont_one", "-i " + fa.string() + " --frac 0.5");
    REQUIRE_FALSE(skc.empty());

    const Container sf(skc);
    REQUIRE(sf.size() == 1);

    const sketch_config_t& cfg = sf.get_config();
    CHECK(cfg.k == 27);
    CHECK(cfg.h == 11);
    CHECK(cfg.canonical);
    CHECK(cfg.frac == doctest::Approx(0.5));
    CHECK(cfg.ppos_v.size() == cfg.h);
    CHECK(cfg.npos_v.size() == size_t(cfg.k - cfg.h));

    const Sketch sk = sf.open(0);
    CHECK(sk.get_rname() == fa.filename().string());
    CHECK(sk.get_nkmers() > 0);
    CHECK(sk.get_card() > 0.0);
    // Random sequence has almost no duplicate k-mers, so rho is in (0, 1].
    CHECK(sk.get_rho() > 0.0);
    CHECK(sk.get_rho() <= 1.0);
    CHECK(sk.get_timestamp() > 0);
    CHECK(sk.get_ntotal_bp() == 200000);
    CHECK(sk.has_buckets());
    CHECK(sk.has_windows());
    CHECK(sk.get_buckets().get_nkmers() == sk.get_nkmers());

    std::filesystem::remove(fa);
    std::filesystem::remove(skc);
  }

  TEST_CASE("h defaults from k independently of -w" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("cont_hdef", 100000, 9);

    // floor(k/2) - 2, whether or not -w was given.
    const auto with_w = run_sketch("cont_hdef_w", "-i " + fa.string() + " -k 23 -w 23");
    REQUIRE_FALSE(with_w.empty());
    const Container sf_w(with_w);
    CHECK(sf_w.get_config().h == 9);
    CHECK(sf_w.get_config().w == 23);

    const auto no_w = run_sketch("cont_hdef_now", "-i " + fa.string() + " -k 23");
    REQUIRE_FALSE(no_w.empty());
    const Container sf_now(no_w);
    CHECK(sf_now.get_config().h == 9);
    CHECK(sf_now.get_config().w == 29); // k + 6

    // k - h <= 16 (enc_t holds 2(k - h) bits), so h clamps up at large k.
    const auto k31 = run_sketch("cont_hdef_k31", "-i " + fa.string() + " -k 31");
    REQUIRE_FALSE(k31.empty());
    CHECK(Container(k31).get_config().h == 15);

    std::filesystem::remove(fa);
    std::filesystem::remove(with_w);
    std::filesystem::remove(no_w);
    std::filesystem::remove(k31);
  }

  TEST_CASE("dist cross mode reports both records" * doctest::skip(!gdiff_available()))
  {
    const auto fa_a = test_util::write_fasta("cont_ns_a", 120000, 21);
    const auto fa_b = test_util::write_fasta("cont_ns_b", 120000, 22);
    const auto skc_a = run_sketch("cont_ns_a", "-i " + fa_a.string() + " -l 500");
    const auto skc_b = run_sketch("cont_ns_b", "-i " + fa_b.string() + " -l 500");
    REQUIRE_FALSE(skc_a.empty());
    REQUIRE_FALSE(skc_b.empty());

    // Both members are opened as query windows, so both sketch names have to reach the output,
    // under the header that is now always written.
    const auto out = test_util::path("cont_ns.tsv");
    const std::string cmd =
      GDIFF_BIN + " dist " + skc_a.string() + " " + skc_b.string() + " > " + out.string() + " 2>/dev/null";
    REQUIRE(std::system(cmd.c_str()) == 0);

    std::ifstream in(out);
    std::string line;
    REQUIRE(std::getline(in, line)); // header
    CHECK(line.rfind("genome_a\tgenome_b\t", 0) == 0);

    REQUIRE(std::getline(in, line)); // the pair's reconciled line
    size_t tabs = 0;
    for (const char c : line)
      tabs += (c == '\t');
    CHECK(tabs == 13); // genome_a, genome_b, then the 12 reconciled columns
    CHECK(line.find(fa_a.filename().string()) == 0);
    CHECK(line.find(fa_b.filename().string()) != std::string::npos);

    std::filesystem::remove(fa_a);
    std::filesystem::remove(fa_b);
    std::filesystem::remove(skc_a);
    std::filesystem::remove(skc_b);
    std::filesystem::remove(out);
  }

  TEST_CASE("dist names a pair canonically whatever the input order" * doctest::skip(!gdiff_available()))
  {
    const auto fa_a = test_util::write_fasta("cont_ord_a", 120000, 31);
    const auto fa_b = test_util::write_fasta("cont_ord_b", 120000, 32);
    const auto skc_a = run_sketch("cont_ord_a", "-i " + fa_a.string() + " -l 500");
    const auto skc_b = run_sketch("cont_ord_b", "-i " + fa_b.string() + " -l 500");
    REQUIRE_FALSE(skc_a.empty());
    REQUIRE_FALSE(skc_b.empty());

    // `cont_ord_a.fa` sorts first, so it is genome_a either way round and swapping the operands
    // has to reproduce the output byte for byte.
    const std::string fwd = run_dist(skc_a.string() + " " + skc_b.string());
    const std::string rev = run_dist(skc_b.string() + " " + skc_a.string());
    REQUIRE_FALSE(fwd.empty());
    CHECK(fwd == rev);

    std::istringstream in(fwd);
    std::string line;
    REQUIRE(std::getline(in, line)); // header
    REQUIRE(std::getline(in, line)); // the pair's reconciled line, keyed canonically
    CHECK(line.rfind(fa_a.filename().string(), 0) == 0);
    CHECK(line.find(fa_b.filename().string()) != std::string::npos);

    // Every sample line names the pair the same way; only the direction column says which of the
    // two is the query.
    const std::string samples = run_dist("--output-samples " + skc_b.string() + " " + skc_a.string());
    std::istringstream sin(samples);
    REQUIRE(std::getline(sin, line)); // header
    size_t n_ab = 0, n_ba = 0;
    while (std::getline(sin, line)) {
      if (line.empty()) continue;
      CHECK(field(line, 0) == "gdiff"); // the config column precedes the pair's key
      CHECK(field(line, 1) == fa_a.filename().string());
      CHECK(field(line, 2) == fa_b.filename().string());
      const std::string dir = field(line, 7);
      CHECK((dir == "ab" || dir == "ba"));
      n_ab += dir == "ab";
      n_ba += dir == "ba";
    }
    CHECK(n_ab > 0);
    CHECK(n_ba > 0);

    std::filesystem::remove(fa_a);
    std::filesystem::remove(fa_b);
    std::filesystem::remove(skc_a);
    std::filesystem::remove(skc_b);
  }

  TEST_CASE("every record of a multi-record container is independently addressable" *
            doctest::skip(!gdiff_available()))
  {
    const auto fa1 = test_util::write_fasta("cont_multi1", 120000, 2);
    const auto fa2 = test_util::write_fasta("cont_multi2", 90000, 3);
    const auto fa3 = test_util::write_fasta("cont_multi3", 60000, 4);
    const auto skc =
      run_sketch("cont_multi", "-i " + fa1.string() + " -i " + fa2.string() + " -i " + fa3.string());
    REQUIRE_FALSE(skc.empty());

    const Container sf(skc);
    REQUIRE(sf.size() == 3);

    // Records are indexed in input order, each at its own offset.
    const std::vector<std::string> expect{
      fa1.filename().string(), fa2.filename().string(), fa3.filename().string()};
    const std::vector<uint64_t> expect_bp{120000, 90000, 60000};
    std::set<uint64_t> offsets;
    for (uint32_t r = 0; r < sf.size(); ++r) {
      const Sketch sk = sf.open(r);
      CHECK(sk.get_rname() == expect[r]);
      CHECK(sk.get_ntotal_bp() == expect_bp[r]);
      CHECK(sk.get_buckets().get_nkmers() > 0);
      offsets.insert(sf.get_entry(r).offset);
    }
    CHECK(offsets.size() == sf.size());

    std::filesystem::remove(fa1);
    std::filesystem::remove(fa2);
    std::filesystem::remove(fa3);
    std::filesystem::remove(skc);
  }

  TEST_CASE("SketchLoad loads only the requested sections" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("cont_part", 150000, 5);
    const auto skc = run_sketch("cont_part", "-i " + fa.string());
    REQUIRE_FALSE(skc.empty());

    const Container sf(skc);
    const Sketch buckets_only = sf.open(0, SketchLoad::Buckets);
    CHECK(buckets_only.has_buckets());
    CHECK_FALSE(buckets_only.has_windows());

    const Sketch windows_only = sf.open(0, SketchLoad::Windows);
    CHECK_FALSE(windows_only.has_buckets());
    CHECK(windows_only.has_windows());

    // Metadata lives in the sketch preamble, so it is present either way.
    CHECK(buckets_only.get_rname() == windows_only.get_rname());
    CHECK(buckets_only.get_rho() == windows_only.get_rho());

    std::filesystem::remove(fa);
    std::filesystem::remove(skc);
  }

  TEST_CASE("pool and seq records agree on window geometry" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("cont_repr", 200000, 6);
    const auto pool = run_sketch("cont_pool", "-i " + fa.string() + "");
    const auto seq = run_sketch("cont_seq", "-i " + fa.string() + " --keep-seq");
    REQUIRE_FALSE(pool.empty());
    REQUIRE_FALSE(seq.empty());

    const Container pf(pool), qf(seq);
    CHECK_FALSE(pf.get_config().keep_seq);
    CHECK(qf.get_config().keep_seq);
    // keep_seq is excluded from compatibility: default and seq of one geometry compare.
    CHECK(compatible_configs(pf.get_config(), qf.get_config()));

    const Sketch ps = pf.open(0), qs = qf.open(0);
    const vec<window_t>& pw = ps.get_windows().wins_v;
    const vec<window_t>& qw = qs.get_windows().wins_v;
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
    const auto fa = test_util::write_fasta("cont_unpack", 150000, 7);
    const auto skc = run_sketch("cont_unpack", "-i " + fa.string() + " --keep-seq -l 500");
    REQUIRE_FALSE(skc.empty());

    const Container sf(skc);
    const Sketch sk = sf.open(0, SketchLoad::Windows);
    const window_sample_t& ws = sk.get_windows();
    REQUIRE_FALSE(ws.wins_v.empty());

    str cseq;
    for (uint64_t wi = 0; wi < std::min<uint64_t>(ws.wins_v.size(), 16); ++wi) {
      ws.packs.unpack(wi, cseq);
      // A window of n k-mers spans n + k - 1 bases.
      CHECK(cseq.size() == ws.wins_v[wi].end - ws.wins_v[wi].start + sk.get_k() - 1);
      CHECK(cseq.find_first_not_of("ACGTN") == str::npos);
    }

    std::filesystem::remove(fa);
    std::filesystem::remove(skc);
  }

  TEST_CASE("canonicalize is a no-op on a canonical sketch" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("cont_canon", 120000, 8);
    const auto skc = run_sketch("cont_canon", "-i " + fa.string());
    REQUIRE_FALSE(skc.empty());

    const Container sf(skc);
    Sketch sk = sf.open(0, SketchLoad::Buckets);
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
    const auto fa = test_util::write_fasta("cont_sa", 120000, 9);
    const auto ag = run_sketch("cont_ag", "-i " + fa.string());
    const auto sa = run_sketch("cont_sa", "-i " + fa.string() + " --strand-aware");
    REQUIRE_FALSE(ag.empty());
    REQUIRE_FALSE(sa.empty());

    CHECK(Container(ag).open(0).is_canonical());
    CHECK_FALSE(Container(sa).open(0).is_canonical());
    // Strandedness changes what a bucket means, so the two must not compare.
    CHECK_FALSE(compatible_configs(Container(ag).get_config(), Container(sa).get_config()));

    std::filesystem::remove(fa);
    std::filesystem::remove(ag);
    std::filesystem::remove(sa);
  }

  TEST_CASE("is_container_file accepts containers and rejects anything else" *
            doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("cont_probe", 60000, 10);
    const auto skc = run_sketch("cont_probe", "-i " + fa.string());
    REQUIRE_FALSE(skc.empty());

    CHECK(is_container_file(skc));
    CHECK_FALSE(is_container_file(fa));
    CHECK_FALSE(is_container_file(test_util::path("cont_does_not_exist.gs")));

    const auto empty = test_util::path("cont_empty.bin");
    std::ofstream(empty).close();
    CHECK_FALSE(is_container_file(empty));

    std::filesystem::remove(fa);
    std::filesystem::remove(skc);
    std::filesystem::remove(empty);
  }

  TEST_CASE("a corrupt magic is rejected rather than misread" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("cont_bad", 60000, 11);
    const auto skc = run_sketch("cont_bad", "-i " + fa.string());
    REQUIRE_FALSE(skc.empty());

    const auto bad = test_util::path("cont_bad_magic.gs");
    std::filesystem::copy_file(skc, bad, std::filesystem::copy_options::overwrite_existing);
    {
      std::fstream f(bad, std::ios::binary | std::ios::in | std::ios::out);
      const uint32_t junk = 0xdeadbeefu;
      f.write(reinterpret_cast<const char*>(&junk), sizeof(junk));
    }
    CHECK_FALSE(is_container_file(bad));
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
    vec<uint64_t> keys_v{
      pack_key(0, 100), pack_key(0, 50), pack_key(3, 200), pack_key(255, 7), pack_key(3, 200)};
    Buckets b;
    b.build(nrows, std::move(keys_v));

    // build() uniques, so the duplicate in bucket 3 collapses.
    CHECK(b.get_nkmers() == 4);
    CHECK(b.get_nnonempty() == 3);
    CHECK(b.get_nrows() == nrows);
    CHECK_FALSE(b.is_empty());

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

    CHECK(b.is_nonempty(0));
    CHECK_FALSE(b.is_nonempty(1));
    CHECK_FALSE(b.is_nonempty(nrows));
  }

  TEST_CASE("a saved table views back identically")
  {
    const uint32_t nrows = 1024;
    std::mt19937_64 rng(42);
    vec<uint64_t> keys_v;
    for (int i = 0; i < 500; ++i)
      keys_v.push_back(pack_key(static_cast<uint32_t>(rng() % nrows), static_cast<enc_t>(rng())));

    Buckets owned;
    owned.build(nrows, vec<uint64_t>(keys_v));

    std::ostringstream os(std::ios::binary);
    owned.save(os);
    const std::string bytes = os.str();
    CHECK(bytes.size() == Buckets::get_byte_size(owned.get_nkmers(), owned.get_nnonempty(), nrows));

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
    CHECK(b.is_empty());
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
    CHECK(src.is_empty());

    const enc_t *beg = nullptr, *end = nullptr;
    REQUIRE(dst.range(5, beg, end));
    CHECK(end - beg == 2);
  }

} // TEST_SUITE

TEST_CASE("-l 0 stores a buckets-only sketch that map accepts" * doctest::skip(!gdiff_available()))
{
  const auto fa = test_util::write_fasta("bonly_ref", 200000, 5);
  const auto skc = run_sketch("bonly_ref", "-i " + fa.string() + " -l 0");
  REQUIRE_FALSE(skc.empty());

  const Container sf(skc);
  CHECK(sf.get_config().tau == 0);
  CHECK(sf.get_config().sample_size == 1000);

  const Sketch sk = sf.open(0, SketchLoad::All);
  CHECK(sk.has_buckets());
  CHECK_FALSE(sk.has_windows());
  CHECK(sk.get_windows().wins_v.empty());

  // map reads only the bucket index, so the buckets-only container is enough.
  const auto qfa = test_util::write_fasta("bonly_query", 50000, 6);
  const std::string cmd =
    GDIFF_BIN + " map -d 0.10 -l 300 --sample-size 0 " + qfa.string() + " " + skc.string() + " >/dev/null 2>&1";
  CHECK(std::system(cmd.c_str()) == 0);

  // dist needs sampled windows and must reject it.
  const std::string dist_cmd = GDIFF_BIN + " dist " + skc.string() + " >/dev/null 2>&1";
  CHECK(std::system(dist_cmd.c_str()) != 0);
}

TEST_SUITE("merge and info")
{

  TEST_CASE("merge concatenates sketches in input order and keeps the configuration" *
            doctest::skip(!gdiff_available()))
  {
    const auto fa_a = test_util::write_fasta("mg_a", 60000, 41);
    const auto fa_b = test_util::write_fasta("mg_b", 60000, 42);
    const auto skc_a = run_sketch("mg_a", "-i " + fa_a.string());
    const auto skc_b = run_sketch("mg_b", "-i " + fa_b.string());
    REQUIRE_FALSE(skc_a.empty());
    REQUIRE_FALSE(skc_b.empty());

    const auto merged = run_merge("mg_ab", {skc_a.string(), skc_b.string()});
    REQUIRE_FALSE(merged.empty());

    const Container ca(skc_a);
    const Container mc(merged);
    CHECK(mc.size() == 2);
    CHECK(mc.get_config().k == ca.get_config().k);
    CHECK(mc.get_config().w == ca.get_config().w);
    CHECK(mc.get_config().h == ca.get_config().h);
    CHECK(mc.get_config().canonical == ca.get_config().canonical);
    CHECK(mc.get_config().tau == ca.get_config().tau);
    CHECK(mc.get_config().sample_size == ca.get_config().sample_size);
    CHECK(mc.get_config().keep_seq == ca.get_config().keep_seq);

    // Both members stay addressable, in the order they were given, with their windows.
    CHECK(mc.open(0).get_rname() == fa_a.filename().string());
    CHECK(mc.open(1).get_rname() == fa_b.filename().string());
    CHECK(mc.open(0).has_windows());

    // A merged container is an ordinary input for the other subcommands.
    CHECK_FALSE(run_dist(merged.string()).empty());
  }

  TEST_CASE("merge accepts one container and keeps repeats as separate sketches" *
            doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("mg_one", 60000, 43);
    const auto skc = run_sketch("mg_one", "-i " + fa.string());
    REQUIRE_FALSE(skc.empty());

    const auto one = run_merge("mg_single", {skc.string()});
    REQUIRE_FALSE(one.empty());
    CHECK(Container(one).size() == 1);

    // Repeats are not deduplicated: the same container twice is two sketches.
    const auto twice = run_merge("mg_twice", {skc.string(), skc.string()});
    REQUIRE_FALSE(twice.empty());
    CHECK(Container(twice).size() == 2);
  }

  TEST_CASE("merge refuses containers whose configuration differs" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("mg_cfg", 60000, 44);
    const auto k27 = run_sketch("mg_k27", "-i " + fa.string());
    const auto k23 = run_sketch("mg_k23", "-i " + fa.string() + " -k 23");
    REQUIRE_FALSE(k27.empty());
    REQUIRE_FALSE(k23.empty());

    CHECK(run_merge("mg_bad", {k27.string(), k23.string()}).empty());
    // A refused merge writes nothing at all.
    CHECK_FALSE(std::filesystem::exists(test_util::path("mg_bad.gs")));
  }

  TEST_CASE("merge refuses to write over one of its inputs" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("mg_over", 60000, 49);
    const auto skc = run_sketch("mg_over", "-i " + fa.string());
    REQUIRE_FALSE(skc.empty());

    // The output stream truncates its path, so a merge into an input would destroy it before the
    // copy finished. The container must come through the refusal byte for byte.
    const std::string before = test_util::read_file(skc);
    const std::string cmd = GDIFF_BIN + " merge -i " + skc.string() + " -o " + skc.string() + " >/dev/null 2>&1";
    CHECK(std::system(cmd.c_str()) != 0);
    CHECK(test_util::read_file(skc) == before);

    // The same input named by a different spelling is still the same file.
    const std::filesystem::path dotted = skc.parent_path() / "." / skc.filename();
    const std::string dotted_cmd =
      GDIFF_BIN + " merge -i " + dotted.string() + " -o " + skc.string() + " >/dev/null 2>&1";
    CHECK(std::system(dotted_cmd.c_str()) != 0);
    CHECK(test_util::read_file(skc) == before);
  }

  TEST_CASE("info reports the configuration and every sketch" * doctest::skip(!gdiff_available()))
  {
    const auto fa_a = test_util::write_fasta("info_a", 60000, 45);
    const auto fa_b = test_util::write_fasta("info_b", 60000, 46);
    const auto skc_a = run_sketch("info_a", "-i " + fa_a.string());
    const auto skc_b = run_sketch("info_b", "-i " + fa_b.string());
    REQUIRE_FALSE(skc_a.empty());
    REQUIRE_FALSE(skc_b.empty());

    const auto merged = run_merge("info_ab", {skc_a.string(), skc_b.string()});
    REQUIRE_FALSE(merged.empty());
    const std::string report = run_info(merged);
    REQUIRE_FALSE(report.empty());

    CHECK(info_value(report, "Sketches:") == "2");
    CHECK(info_value(report, "k (mer len):") == "27");
    CHECK(info_value(report, "w (win len):") == "33");
    CHECK(info_value(report, "h (LSH pos):") == "11");
    CHECK(info_value(report, "canonical:") == "true");
    CHECK(info_value(report, "-l (window len):") == "500");
    CHECK(info_value(report, "--sample-size:") == "1000");
    CHECK(info_value(report, "--keep-seq:") == "false");
    CHECK(info_value(report, "seed:") == "0");

    // One block per sketch, and both names are shown.
    CHECK(count_of(report, "[Sketch ") == 2);
    CHECK(report.find(fa_a.filename().string()) != std::string::npos);
    CHECK(report.find(fa_b.filename().string()) != std::string::npos);
  }

  TEST_CASE("info reports a buckets-only container as having no windows" *
            doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("info_bonly", 60000, 47);
    const auto skc = run_sketch("info_bonly", "-i " + fa.string() + " -l 0");
    REQUIRE_FALSE(skc.empty());

    const std::string report = run_info(skc);
    REQUIRE_FALSE(report.empty());
    CHECK(info_value(report, "Sketches:") == "1");
    CHECK(info_value(report, "-l (window len):") == "0");
    CHECK(info_value(report, "Windows:") == "0");
  }

  TEST_CASE("info rejects a second container and a missing file" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("info_bad", 60000, 48);
    const auto skc = run_sketch("info_bad", "-i " + fa.string());
    REQUIRE_FALSE(skc.empty());

    const auto sink = test_util::path("info_sink.txt");
    const std::string two_inputs =
      GDIFF_BIN + " info -i " + skc.string() + " " + skc.string() + " > " + sink.string() + " 2>&1";
    CHECK(std::system(two_inputs.c_str()) != 0);

    const std::string missing =
      GDIFF_BIN + " info -i " + test_util::path("info_missing.gs").string() + " > " + sink.string() + " 2>&1";
    CHECK(std::system(missing.c_str()) != 0);
  }

} // TEST_SUITE

TEST_SUITE("CLI surface")
{

  TEST_CASE("the same --seed reproduces the sampled windows" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("seed_src", 120000, 51);
    const auto a = run_sketch("seed_a", "-i " + fa.string() + " --seed 7");
    const auto b = run_sketch("seed_b", "-i " + fa.string() + " --seed 7");
    const auto c = run_sketch("seed_c", "-i " + fa.string() + " --seed 8");
    REQUIRE_FALSE(a.empty());
    REQUIRE_FALSE(b.empty());
    REQUIRE_FALSE(c.empty());

    // A container stores a wall-clock timestamp, so compare what the seed actually controls:
    // the retained k-mer set and the sampled window starts.
    const Container ca(a), cb(b), cc(c);
    const vec<uint64_t> wa = window_starts(ca.open(0));
    CHECK_FALSE(wa.empty());
    CHECK(wa == window_starts(cb.open(0)));
    CHECK(wa != window_starts(cc.open(0)));
    CHECK(ca.open(0).get_buckets().get_nkmers() == cb.open(0).get_buckets().get_nkmers());
  }

  TEST_CASE("sketch --frac subsamples the retained k-mers" * doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("frac_src", 200000, 52);
    const auto all = run_sketch("frac_all", "-i " + fa.string());
    const auto half = run_sketch("frac_half", "-i " + fa.string() + " --frac 0.5");
    REQUIRE_FALSE(all.empty());
    REQUIRE_FALSE(half.empty());

    const Container c_all(all), c_half(half);
    CHECK(c_half.get_config().frac == doctest::Approx(0.5));
    const uint64_t n_all = c_all.open(0).get_buckets().get_nkmers();
    const uint64_t n_half = c_half.open(0).get_buckets().get_nkmers();
    CHECK(n_half > 0);
    CHECK(n_half < n_all);
  }

  TEST_CASE("sketch --input-list reads optional name and path columns" *
            doctest::skip(!gdiff_available()))
  {
    const auto fa = test_util::write_fasta("il_src", 60000, 53);
    const auto list = test_util::path("il_list.txt");
    std::ofstream(list) << "custom_name\t" << fa.string() << "\n";

    const auto skc = run_sketch("il_out", "--input-list " + list.string());
    REQUIRE_FALSE(skc.empty());
    const Container c(skc);
    CHECK(c.size() == 1);
    CHECK(c.open(0).get_rname() == "custom_name"); // the name column wins over the file name
  }

  TEST_CASE("dist --list-a and --list-b select the same sets as positionals" *
            doctest::skip(!gdiff_available()))
  {
    const auto fa_a = test_util::write_fasta("dl_a", 120000, 54);
    const auto fa_b = test_util::write_fasta("dl_b", 120000, 55);
    const auto skc_a = run_sketch("dl_a", "-i " + fa_a.string());
    const auto skc_b = run_sketch("dl_b", "-i " + fa_b.string());
    REQUIRE_FALSE(skc_a.empty());
    REQUIRE_FALSE(skc_b.empty());

    const auto list_a = test_util::path("dl_a.txt");
    const auto list_b = test_util::path("dl_b.txt");
    std::ofstream(list_a) << skc_a.string() << "\n" << skc_b.string() << "\n";
    std::ofstream(list_b) << skc_a.string() << "\n";

    // Two containers in one list is the same within-set comparison as two positionals.
    const std::string positional = run_dist(skc_a.string() + " " + skc_b.string());
    REQUIRE_FALSE(positional.empty());
    CHECK(run_dist("--list-a " + list_a.string()) == positional);

    // Two lists take the cross product, so set B restricts which pairs are reported.
    const std::string cross = run_dist("--list-a " + list_a.string() + " --list-b " + list_b.string());
    REQUIRE_FALSE(cross.empty());
    CHECK(cross != positional);
    CHECK(cross.find(fa_a.filename().string()) != std::string::npos);
    CHECK(cross.find(fa_b.filename().string()) != std::string::npos);
  }

  TEST_CASE("dist -o writes exactly the stdout report" * doctest::skip(!gdiff_available()))
  {
    const auto fa_a = test_util::write_fasta("do_a", 120000, 56);
    const auto fa_b = test_util::write_fasta("do_b", 120000, 57);
    const auto skc_a = run_sketch("do_a", "-i " + fa_a.string());
    const auto skc_b = run_sketch("do_b", "-i " + fa_b.string());
    REQUIRE_FALSE(skc_a.empty());
    REQUIRE_FALSE(skc_b.empty());

    const auto out = test_util::path("dist_o.tsv");
    const std::string cmd =
      GDIFF_BIN + " dist -o " + out.string() + " " + skc_a.string() + " " + skc_b.string() + " >/dev/null 2>&1";
    REQUIRE(std::system(cmd.c_str()) == 0);
    CHECK(test_util::read_file(out) == run_dist(skc_a.string() + " " + skc_b.string()));
  }

} // TEST_SUITE
