#ifndef _SKETCH2_HPP
#define _SKETCH2_HPP

// sketch2/dist2 experiment: multi-genome bundles, all-vs-all distances.
#include "CLI11.hpp"
#include "hm.hpp"
#include "lshf.hpp"
#include "sketch.hpp"
#include "types.hpp"
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <utility>

// Binary sketch2 container format ("GSK5"), versioned:
//   [u32] magic "GSK5"  [u32] version  [u64] nsketches
//   shared params block (seed, timestamp, k/w/h, canonical, nrows, frac,
//                        tau, bin_shift, sample_size, ppos, npos)
//   index block: nsketches * sketch2_entry_t (see below) -- one small read
//                so readers never walk the file
//   per record: rname_len/rname, timestamp, genome_length_bp, nvalid_bases,
//                rho, has_sfhm, SFHM ?, windows section (tau/bin_shift/nwins/
//                per-window meta incl nvalid + pool_fw + pool_rc).
// k-mer hashes live in global pools; shared LSH geometry lives in the file header.
static constexpr uint32_t SKETCH2_MAGIC = 0x35304B53u; // "GSK5"
static constexpr uint32_t SKETCH2_VERSION = 1;

#pragma pack(push, 1)
// One record's location + section ranges inside the mmap'd file. Offsets are
// absolute file offsets (except len, relative to the record start).
struct sketch2_entry_t
{
  uint64_t offset = 0;      // start of the per-record preamble
  uint64_t len = 0;         // record size (preamble + SFHM + windows)
  uint64_t buckets_off = 0; // abs offset of SFHM bytes (0 for a record without SF, e.g. old --windows-only files)
  uint64_t buckets_len = 0; // SFHM byte length
  uint64_t windows_off = 0; // abs offset of windows+pools section
  uint64_t windows_len = 0; // windows section byte length
};
#pragma pack(pop)

// File-level header + index, read once by readers. dist2 passes this into
// load_from_offset so each record can rebuild its LSH geometry without
// re-parsing the whole container.
struct sketch2_index_t
{
  uint64_t version = SKETCH2_VERSION;
  uint64_t nsketches = 0;
  uint64_t seed = 0;
  uint64_t timestamp = 0;
  uint8_t k = 0;
  uint8_t w = 0;
  uint8_t h = 0;
  bool canonical = true;
  uint32_t nrows = 0;
  double frac = 1.0;
  uint64_t tau = 0;
  uint64_t bin_shift = 0;
  uint64_t sample_size = 0;
  vec<uint8_t> ppos; // size h
  vec<uint8_t> npos; // size k - h
  vec<sketch2_entry_t> records;
};

// Pack LSH bucket index and residual encoding into one 64-bit hash value.
inline uint64_t pack_win_hash(uint32_t bix, enc_t enc) noexcept { return (uint64_t(bix) << 32) | uint64_t(enc); }

inline uint32_t win_hash_bix(uint64_t h) noexcept { return static_cast<uint32_t>(h >> 32); }

inline enc_t win_hash_enc(uint64_t h) noexcept { return static_cast<enc_t>(h & 0xffffffffu); }

// Window metadata only; k-mer hashes live in the global pools.
struct win2_t
{
  str qid;
  uint64_t start = 0;     // 0-based k-mer start
  uint64_t end = 0;       // exclusive k-mer end
  uint32_t nvalid_fw = 0; // valid k-mers seen in the window (N-excluded)
  uint32_t nvalid_rc = 0; // same for reverse strand (0 in canonical mode)
};

// All retained window k-mers, sorted by packed LSH hash (bucket then enc).
// win_ix[i] is the window that hashes[i] belongs to.
// Owned (hashes/win_ix) or a borrowed view (hashes_v/win_ix_v) into a live
// external buffer (e.g. mmap) that the caller keeps alive.
struct win_hash_pool_t
{
  vec<uint64_t> hashes;
  vec<uint32_t> win_ix;
  const uint64_t* hashes_v = nullptr;
  const uint32_t* win_ix_v = nullptr;
  uint64_t n_ = 0;

  void clear()
  {
    hashes.clear();
    win_ix.clear();
    hashes_v = nullptr;
    win_ix_v = nullptr;
    n_ = 0;
  }

  bool is_view() const noexcept { return hashes_v != nullptr; }
  uint64_t size() const noexcept { return is_view() ? n_ : hashes.size(); }
  const uint64_t* hashes_ptr() const noexcept { return is_view() ? hashes_v : hashes.data(); }
  const uint32_t* win_ix_ptr() const noexcept { return is_view() ? win_ix_v : win_ix.data(); }
};

struct win_sample_t
{
  vec<win2_t> wins_v;
  win_hash_pool_t pool_fw;
  win_hash_pool_t pool_rc; // empty when canonical
};

// Selective load: dist2 only needs query windows + reference buckets.
enum class Sketch2Part : uint8_t
{
  Buckets = 1,
  Windows = 2,
  All = 3
};

class Sketch2
{
public:
  static constexpr uint32_t INVALID_BIX = std::numeric_limits<uint32_t>::max();

  explicit Sketch2(std::filesystem::path sketch_path);
  // mmap-backed selective load. Skipped SFHM pages are not faulted in.
  // Rebuilds LSH geometry from the file-level index (shared params block).
  void load_from_offset(const sketch2_entry_t& e, const sketch2_index_t& idx, Sketch2Part part = Sketch2Part::All);
  // madvise a byte range of the mapped file (e.g. MADV_DONTNEED to evict a
  // reference's bucket pages once its job batch is done). No-op without a map.
  void advise_range(uint64_t off, uint64_t len, int advice) const noexcept;

  uint32_t validate_bucket_ix(uint32_t bix) const noexcept { return bix < nrows ? bix : INVALID_BIX; }
  bool scan_bucket(uint32_t bix, enc_t enc_lr, uint32_t& hdist_min) const noexcept;
  // True iff LSH bucket bix is in-range and contains at least one encoding.
  bool bucket_nonempty(uint32_t bix) const noexcept
  {
    if (bix >= nrows || nonempty_bits.empty()) return false;
    return (nonempty_bits[bix >> 6] >> (bix & 63)) & 1ull;
  }

  str get_rname() const { return rname; }
  double get_rho() const { return rho; }
  bool is_canonical() const { return canonical; }
  uint64_t get_tau() const { return tau; }
  uint64_t get_bin_shift() const { return bin_shift; }
  uint64_t get_nwinmers() const { return nwinmers; }
  uint64_t get_genome_length_bp() const { return genome_length_bp; }
  uint64_t get_nvalid_bases() const { return nvalid_bases; }
  const vec<win2_t>& get_wins() const { return wins_v; }
  const win_hash_pool_t& get_pool_fw() const { return pool_fw; }
  const win_hash_pool_t& get_pool_rc() const { return pool_rc; }
  sfhm_sptr_t get_sfhm_sptr() const { return sfhm; }
  lshf_sptr_t get_lshf_sptr() const { return lshf; }
  uint32_t get_nrows() const { return nrows; }
  uint8_t get_k() const { return k; }
  uint8_t get_h() const { return h; }
  uint8_t get_w() const { return w; }
  bool has_buckets() const { return sfhm != nullptr; }
  bool has_windows() const { return loaded_windows; }

  // True when LSH geometry and window sampling parameters match.
  bool compatible_with(const Sketch2& other) const;

private:
  struct MappedFile;

  static const char* need(const char* p, const char* end, size_t n, const str& what);
  static const char* skip_sfhm_bytes(const char* p, const char* end);
  static const char* skip_windows_bytes(const char* p, const char* end, uint64_t& tau_out, uint64_t& bin_shift_out);
  const char* read_windows_mem(const char* p, const char* end);
  static const char* view_pool_mem(const char* p, const char* end, win_hash_pool_t& pool);
  void build_nonempty_bitmap();
  static void write_pool(std::ostream& stream, const win_hash_pool_t& pool);

  str rname;
  uint8_t k = 0;
  uint8_t w = 0;
  uint8_t h = 0;
  double rho = 1.0;
  uint32_t nrows = 0;
  bool canonical = true;
  uint64_t timestamp = 0;
  uint64_t genome_length_bp = 0;
  uint64_t nvalid_bases = 0;
  uint64_t tau = 0;
  uint64_t bin_shift = 0;
  uint64_t nwinmers = 0;
  bool loaded_windows = false;
  lshf_sptr_t lshf = nullptr;
  sfhm_sptr_t sfhm = nullptr;
  vec<uint64_t> nonempty_bits; // bitset over [0, nrows): 1 => nonempty bucket
  vec<win2_t> wins_v;
  win_hash_pool_t pool_fw;
  win_hash_pool_t pool_rc;
  std::filesystem::path sketch_path;
  std::shared_ptr<MappedFile> mapped;

  friend class Sketch2SC;
  friend sketch2_index_t read_sketch2_index(const std::filesystem::path& sketch_path);
};

// Parse the file header (shared params + index block) in one read. Readers use
// the returned offsets for selective per-record loading and page policy.
sketch2_index_t read_sketch2_index(const std::filesystem::path& sketch_path);

// Read a list file of input paths (one per line; empty lines and lines whose
// first non-space char is '#' skipped). A "name\tpath" line sets an explicit
// record name; a bare path keeps its basename. Shared by sketch2 (--input-list)
// and dist2. Returns the aligned name/path vectors (both length == n lines).
void read_path_list(const std::filesystem::path& list_path, std::vector<str>& paths, std::vector<str>& names);

class Sketch2SC : public BaseLSH
{
public:
  explicit Sketch2SC(CLI::App& sc);
  void process();
  bool validate_configuration();

private:
  // Emit the GSK5 file header: magic/version/nsketches + shared params block +
  // an index-block placeholder that process() patches in place at the end.
  void write_file_header(std::ostream& stream, uint64_t nsketches);
  // Serialize one per-record sketch2 (rname/metadata + SFHM + windows) into
  // `bytes`, returning the entry with section offsets relative to its start.
  sketch2_entry_t write_record(std::string& bytes,
                               const str& rname,
                               uint64_t timestamp,
                               uint64_t genome_length_bp,
                               uint64_t nvalid_bases,
                               double rho,
                               const sfhm_sptr_t& sketch_sfhm,
                               const win_sample_t& sample);
  void write_windows(std::ostream& stream, const win_sample_t& sample);
  win_sample_t sample_windows(const str& input_path, uint64_t& genome_length_bp, uint64_t& nvalid_bases);

  std::vector<str> paths_v;
  std::vector<str> rnames_v; // optional per-input record names (from --input-list "name\tpath")
  std::vector<double> rho_v;
  std::filesystem::path input_list_path;
  std::filesystem::path sketch_path;
  uint64_t tau = 0;
  uint64_t sample_size = 200;
  uint64_t bin_shift = 0;
};

#endif
