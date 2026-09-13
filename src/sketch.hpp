#ifndef _SKETCH_HPP
#define _SKETCH_HPP

#include "CLI11.hpp"
#include "hm.hpp"
#include "lshf.hpp"
#include "msg.hpp"
#include "stils.hpp"
#include "types.hpp"
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>

// Owns one read-only mmap of a container file; defined in sketch.cpp.
class MappedFile;

// GDSK: indexed container of per-genome sketches (header, index, records).
static constexpr uint32_t GDSK_MAGIC = 0x4B534447u; // "GDSK"
static constexpr uint32_t GDSK_VERSION = 1;

// How a record stores its sampled windows.
enum class WinRepr : uint8_t
{
  // One bucket resolve per distinct bucket, no k-mer re-derivation. Freezes -l.
  Pool = 0,
  // 2-bit packed bases plus an N-mask: ~12x smaller, re-derives hashes per query.
  Seq = 1
};

#pragma pack(push, 1)
// Absolute file offsets; a zero section length means the section is absent.
struct record_entry_t
{
  uint64_t offset = 0;
  uint64_t len = 0;
  uint64_t buckets_off = 0;
  uint64_t buckets_len = 0;
  uint64_t windows_off = 0;
  uint64_t windows_len = 0;
};
#pragma pack(pop)

// The per-file configuration every record in a container shares.
struct sketch_config_t
{
  uint64_t seed = 0;
  uint64_t timestamp = 0;
  uint8_t k = 0;
  uint8_t w = 0;
  uint8_t h = 0;
  bool canonical = true;
  uint32_t nrows = 0;
  double frac = 1.0;
  uint64_t tau = 0;         // sampled window length in k-mers (0 = no windows)
  uint64_t sample_size = 0; // windows sampled per genome
  WinRepr win_repr = WinRepr::Pool;
  vec<uint8_t> ppos; // size h
  vec<uint8_t> npos; // size k - h

  // win_repr is excluded: a direction reads only the query side's windows.
  bool compatible_with(const sketch_config_t& other) const;
};

struct sketch_index_t
{
  sketch_config_t cfg;
  vec<record_entry_t> records;
};

// Window metadata; the payload arrays alongside share its index.
struct window_t
{
  str qid;
  uint64_t start = 0;     // 0-based k-mer start in the source sequence
  uint64_t end = 0;       // exclusive k-mer end
  uint32_t nvalid_fw = 0; // valid (non-ambiguous) k-mers in the window
  uint32_t nvalid_rc = 0; // same on the reverse strand; 0 when canonical
};

// Window k-mers sorted by packed key; owned while building, viewed once mapped.
struct hash_pool_t
{
  vec<uint64_t> hashes;
  vec<uint16_t> win_ix;
  const uint64_t* hashes_view = nullptr;
  const uint16_t* win_ix_view = nullptr;
  uint64_t n_view = 0;

  void clear()
  {
    hashes.clear();
    win_ix.clear();
    hashes_view = nullptr;
    win_ix_view = nullptr;
    n_view = 0;
  }

  [[nodiscard]] bool is_view() const noexcept { return hashes_view != nullptr; }
  [[nodiscard]] uint64_t size() const noexcept { return is_view() ? n_view : hashes.size(); }
  [[nodiscard]] const uint64_t* hashes_ptr() const noexcept
  {
    return is_view() ? hashes_view : hashes.data();
  }
  [[nodiscard]] const uint16_t* win_ix_ptr() const noexcept
  {
    return is_view() ? win_ix_view : win_ix.data();
  }
};

// Concatenated 2-bit bases; base_off[i] is window i's bit-pair offset.
struct seq_pack_t
{
  vec<uint64_t> packed;
  vec<uint64_t> nmask;
  vec<uint64_t> base_off; // nwins + 1 entries
  const uint64_t* packed_view = nullptr;
  const uint64_t* nmask_view = nullptr;
  const uint64_t* base_off_view = nullptr;
  uint64_t packed_words_view = 0;
  uint64_t nmask_words_view = 0;

  void clear()
  {
    packed.clear();
    nmask.clear();
    base_off.clear();
    packed_view = nullptr;
    nmask_view = nullptr;
    base_off_view = nullptr;
    packed_words_view = 0;
    nmask_words_view = 0;
  }

  [[nodiscard]] bool is_view() const noexcept { return base_off_view != nullptr; }
  [[nodiscard]] const uint64_t* packed_ptr() const noexcept
  {
    return is_view() ? packed_view : packed.data();
  }
  [[nodiscard]] const uint64_t* nmask_ptr() const noexcept
  {
    return is_view() ? nmask_view : (nmask.empty() ? nullptr : nmask.data());
  }
  [[nodiscard]] const uint64_t* base_off_ptr() const noexcept
  {
    return is_view() ? base_off_view : base_off.data();
  }
  // Expand window wi into cseq as ACGT characters plus 'N' where masked.
  void unpack(uint64_t wi, str& cseq) const;
};

struct window_sample_t
{
  vec<window_t> wins;
  hash_pool_t pool_fw;
  hash_pool_t pool_rc; // empty when canonical
  seq_pack_t packs;

  void clear()
  {
    wins.clear();
    pool_fw.clear();
    pool_rc.clear();
    packs.clear();
  }
};

// Selective load: a reference only needs buckets, a query only needs windows.
enum class SketchPart : uint8_t
{
  Buckets = 1,
  Windows = 2,
  All = 3
};

// A read-only view of one record, keeping the container's mapping alive.
class Sketch
{
public:
  Sketch() = default;
  // Built from a FASTA rather than read from a container.
  Sketch(const sketch_config_t& cfg, str rname, Buckets&& buckets, uint64_t nkmers, double rho);

  [[nodiscard]] const str& get_rname() const noexcept { return rname; }
  [[nodiscard]] double get_rho() const noexcept { return rho; }
  [[nodiscard]] uint64_t get_nkmers() const noexcept { return nkmers; }
  [[nodiscard]] double get_card_est() const noexcept { return card_est; }
  [[nodiscard]] uint64_t get_timestamp() const noexcept { return timestamp; }
  [[nodiscard]] uint64_t get_genome_bp() const noexcept { return genome_bp; }
  [[nodiscard]] uint64_t get_nvalid_bases() const noexcept { return nvalid_bases; }
  [[nodiscard]] const sketch_config_t& get_config() const noexcept { return cfg; }
  [[nodiscard]] bool is_canonical() const noexcept { return cfg.canonical; }
  [[nodiscard]] uint8_t get_k() const noexcept { return cfg.k; }
  [[nodiscard]] uint8_t get_h() const noexcept { return cfg.h; }
  [[nodiscard]] uint8_t get_w() const noexcept { return cfg.w; }
  [[nodiscard]] uint32_t get_nrows() const noexcept { return cfg.nrows; }
  [[nodiscard]] uint64_t get_tau() const noexcept { return cfg.tau; }
  [[nodiscard]] const lshf_sptr_t& get_lshf_sptr() const noexcept { return lshf; }
  [[nodiscard]] const Buckets& get_buckets() const noexcept { return buckets; }
  [[nodiscard]] const window_sample_t& get_windows() const noexcept { return windows; }
  [[nodiscard]] const vec<window_t>& get_wins() const noexcept { return windows.wins; }
  [[nodiscard]] bool has_buckets() const noexcept { return !buckets.empty(); }
  [[nodiscard]] bool has_windows() const noexcept { return loaded_windows; }

  // Rebuild the buckets over canonical k-mers; a no-op when already canonical.
  void canonicalize();

private:
  str rname;
  double rho = 1.0;
  uint64_t nkmers = 0;
  double card_est = 0.0;
  uint64_t timestamp = 0;
  uint64_t genome_bp = 0;
  uint64_t nvalid_bases = 0;
  bool loaded_windows = false;
  sketch_config_t cfg;
  lshf_sptr_t lshf;
  Buckets buckets;
  window_sample_t windows;
  std::shared_ptr<const MappedFile> mapped; // keeps the views above alive

  friend class SketchFile;
};

// A memory-mapped GDSK container; every Sketch it opens shares the mapping.
class SketchFile
{
public:
  explicit SketchFile(std::filesystem::path path);

  [[nodiscard]] const sketch_index_t& get_index() const noexcept { return idx; }
  [[nodiscard]] const sketch_config_t& get_config() const noexcept { return idx.cfg; }
  [[nodiscard]] uint32_t size() const noexcept { return static_cast<uint32_t>(idx.records.size()); }
  [[nodiscard]] const std::filesystem::path& get_path() const noexcept { return path; }

  // Load one record; skipped sections are never faulted in.
  [[nodiscard]] Sketch open(uint32_t rec, SketchPart part = SketchPart::All) const;
  // Advise on a byte range, e.g. MADV_DONTNEED to drop a reference's pages.
  void advise(uint64_t off, uint64_t len, int advice) const noexcept;

private:
  std::filesystem::path path;
  std::shared_ptr<const MappedFile> mapped;
  sketch_index_t idx;
};

// True when `path` starts with the GDSK magic.
bool is_sketch_file(const std::filesystem::path& path);

// Write the header, leaving the stream where the index block starts.
void write_sketch_header(std::ostream& os, const sketch_config_t& cfg, uint64_t nsketches);

// One path per line, '#' lines skipped; "name<TAB>path" names the record.
void read_path_list(const std::filesystem::path& list_path, vec<str>& paths, vec<str>& names);

class BaseLSH
{
public:
  void set_lshf();
  void set_nrows();
  void set_sketch_defaults()
  {
    k = 27;
    w = k + 6;
    h = 11;
    frac = 1.0;
    canonical = true;
    nrows = uint32_t(1) << (2 * h); // recomputed by set_nrows()
  }

protected:
  uint8_t k;
  uint8_t w;
  uint8_t h;
  double frac; // FracMin LSH on top of minimizers: keep if LSH(x) < frac * 2^(2h)
  bool canonical;
  uint32_t nrows;
  lshf_sptr_t lshf = nullptr;
};

class SketchSC : public BaseLSH
{
public:
  explicit SketchSC(CLI::App& sc);
  void process();
  bool validate_configuration();

private:
  // Header plus an index placeholder that process() patches at flush time.
  void write_file_header(std::ostream& os, uint64_t nsketches);
  // Section offsets come out record-relative; process() makes them absolute.
  record_entry_t write_record(str& bytes,
                              const str& rname,
                              uint64_t timestamp,
                              uint64_t genome_bp,
                              uint64_t nvalid_bases,
                              uint64_t nkmers,
                              double card_est,
                              double rho,
                              const Buckets& buckets,
                              const window_sample_t& sample);
  void write_windows(std::ostream& os, const window_sample_t& sample);
  window_sample_t sample_windows(const str& input_path,
                                 uint64_t& genome_bp,
                                 uint64_t& nvalid_bases);
  sketch_config_t make_config(uint64_t timestamp) const;

  vec<str> paths_v;
  vec<str> rnames_v; // optional per-input names from --input-list
  std::filesystem::path input_list_path;
  std::filesystem::path sketch_path;
  uint64_t tau = 500;
  uint64_t sample_size = 1000;
  WinRepr win_repr = WinRepr::Pool;
};

// Build a bucket table from a reference file, with its k-mer count and rho.
struct built_sketch_t
{
  Buckets buckets;
  uint64_t nkmers = 0;
  double card_est = 0.0;
  double rho = 1.0;
};

built_sketch_t build_buckets(const str& input_path,
                             const lshf_sptr_t& lshf,
                             uint8_t w,
                             uint32_t nrows,
                             bool canonical);

#endif
