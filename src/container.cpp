#include "container.hpp"

#include <algorithm>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include "lshf.hpp"
#include "serialize.hpp"

FileMap::~FileMap()
{
  if (ptr && ptr != MAP_FAILED && len) munmap(ptr, len);
  if (fd >= 0) close(fd);
}

void FileMap::advise(uint64_t off, uint64_t nbytes, int advice) const noexcept
{
  if (off >= len) return;
  const size_t n = static_cast<size_t>(std::min<uint64_t>(nbytes, len - off));
  if (n == 0) return;
  ::madvise(const_cast<char*>(begin()) + off, n, advice);
}

std::shared_ptr<const FileMap> FileMap::open(const std::filesystem::path& path)
{
  auto m = std::make_shared<FileMap>();
  m->fd = ::open(path.c_str(), O_RDONLY);
  if (m->fd < 0) error_exit("Cannot open container: " + path.string());
  struct stat st
  {};
  if (fstat(m->fd, &st) != 0) error_exit("Cannot fstat container: " + path.string());
  m->len = static_cast<size_t>(st.st_size);
  if (m->len == 0) error_exit("Empty container: " + path.string());
  m->ptr = mmap(nullptr, m->len, PROT_READ, MAP_PRIVATE, m->fd, 0);
  if (m->ptr == MAP_FAILED) error_exit("mmap failed for container: " + path.string());
  return m;
}

bool is_container_file(const std::filesystem::path& path)
{
  std::ifstream in(path, std::ios::binary);
  if (!in) return false;
  uint32_t vgskey = 0;
  in.read(reinterpret_cast<char*>(&vgskey), sizeof(uint32_t));
  return in.good() && vgskey == container_vgskey;
}

void read_container_header(const char*& p,
                           const char* end,
                           sketch_config_t& cfg,
                           uint64_t& nsketches,
                           const std::filesystem::path& path)
{
  const uint32_t vgskey = read_trivial<uint32_t>(p, end, "vgskey");
  const uint32_t version = read_trivial<uint32_t>(p, end, "version");
  if (vgskey != container_vgskey) error_exit("Not a gdiff container: " + path.string());
  if (version != container_version) {
    error_exit(concat_msg(
      "Container format version ", version, " != ", container_version, "; re-run `gdiff sketch`: ", path.string()));
  }
  nsketches = read_trivial<uint64_t>(p, end, "nsketches");

  cfg.seed = read_trivial<uint64_t>(p, end, "seed");
  cfg.timestamp = read_trivial<uint64_t>(p, end, "timestamp");
  cfg.k = read_trivial<uint8_t>(p, end, "k");
  cfg.w = read_trivial<uint8_t>(p, end, "w");
  cfg.h = read_trivial<uint8_t>(p, end, "h");
  cfg.canonical = read_trivial<uint8_t>(p, end, "canonical") != 0;
  cfg.keep_seq = read_trivial<uint8_t>(p, end, "keep_seq") != 0;
  p = need(p, end, 3, "config padding");
  p += 3;
  cfg.nrows = read_trivial<uint32_t>(p, end, "nrows");
  p = need(p, end, 4, "config padding");
  p += 4;
  cfg.frac = read_trivial<double>(p, end, "frac");
  cfg.tau = read_trivial<uint64_t>(p, end, "tau");
  cfg.sample_size = read_trivial<uint64_t>(p, end, "sample_size");

  if (cfg.h > cfg.k || cfg.k == 0) error_exit("Corrupt sketch LSH config: " + path.string());
  p = need(p, end, cfg.k, "LSH positions");
  cfg.ppos_v.assign(p, p + cfg.h);
  p += cfg.h;
  cfg.npos_v.assign(p, p + (cfg.k - cfg.h));
  p += (cfg.k - cfg.h);
  p = need(p, end, round_up_word(cfg.k) - cfg.k, "LSH position padding");
  p += round_up_word(cfg.k) - cfg.k;
}

void read_container_index(const char*& p,
                          const char* end,
                          uint64_t nsketches,
                          vec<scentry>& entries_v,
                          const std::filesystem::path& path)
{
  entries_v.resize(static_cast<size_t>(nsketches));
  for (scentry& r : entries_v) {
    p = need(p, end, sizeof(scentry), "index entry");
    std::memcpy(&r, p, sizeof(scentry));
    p += sizeof(scentry);
  }
}

void write_container_header(std::ostream& os, const sketch_config_t& cfg, uint64_t nsketches)
{
  write_trivial(os, container_vgskey);
  write_trivial(os, container_version);
  write_trivial(os, nsketches);
  write_trivial(os, cfg.seed);
  write_trivial(os, cfg.timestamp);
  write_trivial(os, cfg.k);
  write_trivial(os, cfg.w);
  write_trivial(os, cfg.h);
  write_trivial(os, static_cast<uint8_t>(cfg.canonical ? 1 : 0));
  write_trivial(os, static_cast<uint8_t>(cfg.keep_seq ? 1 : 0));
  const char pad[4] = {};
  os.write(pad, 3);
  write_trivial(os, cfg.nrows);
  os.write(pad, 4);
  write_trivial(os, cfg.frac);
  write_trivial(os, cfg.tau);
  write_trivial(os, cfg.sample_size);
  os.write(reinterpret_cast<const char*>(cfg.ppos_v.data()), static_cast<std::streamsize>(cfg.ppos_v.size()));
  os.write(reinterpret_cast<const char*>(cfg.npos_v.data()), static_cast<std::streamsize>(cfg.npos_v.size()));
  pad_to_word(os);
}

Container::Container(std::filesystem::path path)
  : path(std::move(path))
{
  map = FileMap::open(this->path);
  const char* p = map->begin();
  const char* end = map->end();
  uint64_t nsketches = 0;
  read_container_header(p, end, cfg, nsketches, this->path);
  read_container_index(p, end, nsketches, entries_v, this->path);
  lshf = std::make_shared<LSHF>(cfg.ppos_v, cfg.npos_v);
}

void Container::advise(uint64_t off, uint64_t len, int advice) const noexcept
{
  if (map) map->advise(off, len, advice);
}
