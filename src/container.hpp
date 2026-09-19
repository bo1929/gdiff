#ifndef _CONTAINER_HPP
#define _CONTAINER_HPP

#include "types.hpp"
#include <filesystem>
#include <memory>
#include <ostream>

inline constexpr uint32_t container_vgskey = 0x4B534447u;
inline constexpr uint32_t container_version = 1;

#pragma pack(push, 1)
// Absolute file offsets; a zero section length means the section is absent.
struct sketch_entry_t
{
  uint64_t offset = 0;
  uint64_t len = 0;
  uint64_t buckets_off = 0;
  uint64_t buckets_len = 0;
  uint64_t windows_off = 0;
  uint64_t windows_len = 0;
};
#pragma pack(pop)

// The per-file configuration every sketch in a container shares.
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
  // false (default): pre-resolved (bucket, encoding) keys, one resolve per distinct bucket.
  // true: 2-bit packed bases plus an N-mask, ~12x smaller but re-derives hashes per query.
  bool keep_seq = false;
  vec<uint8_t> ppos_v; // size h
  vec<uint8_t> npos_v; // size k - h
};

// Selective load: a reference only needs buckets, a query only needs windows.
enum class SketchLoad : uint8_t
{
  Buckets = 1,
  Windows = 2,
  All = 3
};

// Owns one read-only mmap of a container file; every view into it shares this owner.
class FileMap
{
public:
  FileMap() = default;
  FileMap(const FileMap&) = delete;
  FileMap& operator=(const FileMap&) = delete;
  ~FileMap();

  [[nodiscard]] const char* begin() const noexcept { return static_cast<const char*>(ptr); }
  [[nodiscard]] const char* end() const noexcept { return begin() + len; }
  [[nodiscard]] size_t size() const noexcept { return len; }

  // Advise on a byte range, e.g. MADV_DONTNEED to drop a reference's pages.
  void advise(uint64_t off, uint64_t nbytes, int advice) const noexcept;

  static std::shared_ptr<const FileMap> open(const std::filesystem::path& path);

private:
  int fd = -1;
  void* ptr = nullptr;
  size_t len = 0;
};

class Container
{
public:
  explicit Container(std::filesystem::path path);

  [[nodiscard]] const sketch_config_t& get_config() const noexcept { return cfg; }
  [[nodiscard]] const vec<sketch_entry_t>& get_entries() const noexcept { return entries_v; }
  [[nodiscard]] const sketch_entry_t& get_entry(uint32_t rec) const noexcept { return entries_v[rec]; }
  [[nodiscard]] uint32_t size() const noexcept { return static_cast<uint32_t>(entries_v.size()); }
  [[nodiscard]] const std::filesystem::path& get_path() const noexcept { return path; }

  // Load one sketch; skipped sections are never faulted in.
  [[nodiscard]] Sketch open(uint32_t rec, SketchLoad part = SketchLoad::All) const;
  void advise(uint64_t off, uint64_t len, int advice) const noexcept;

private:
  std::filesystem::path path;
  std::shared_ptr<const FileMap> map;
  sketch_config_t cfg;
  lshf_sptr_t lshf; // built once; every Sketch from this container shares it
  vec<sketch_entry_t> entries_v;
};

// True when `path` starts with the container vgskey.
bool is_container_file(const std::filesystem::path& path);

// Parse a container's header; leaves `p` at the index block.
void read_container_header(const char*& p,
                           const char* end,
                           sketch_config_t& cfg,
                           uint64_t& nsketches,
                           const std::filesystem::path& path);

// Parse the index block into one entry per sketch.
void read_container_index(const char*& p,
                          const char* end,
                          uint64_t nsketches,
                          vec<sketch_entry_t>& entries_v,
                          const std::filesystem::path& path);

// Write the header, leaving the stream where the index block starts.
void write_container_header(std::ostream& os, const sketch_config_t& cfg, uint64_t nsketches);

#endif
