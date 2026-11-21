#include "merge.hpp"

#include "utils.hpp"
#include "container.hpp"
#include "msg.hpp"
#include "sketch.hpp"

#include <filesystem>
#include <fstream>
#include <system_error>

namespace {

  namespace fs = std::filesystem;

  void reject_input_as_output(const fs::path& out_path, const vec<str>& in_path_v)
  {
    std::error_code ec;
    const fs::path out = fs::weakly_canonical(out_path, ec);
    if (ec) return;
    for (const str& in_path : in_path_v) {
      const fs::path in = fs::weakly_canonical(in_path, ec);
      if (!ec && in == out) error_exit("Output container is also an input: " + in_path);
    }
  }

  vec<Container> open_inputs(const vec<str>& in_path_v)
  {
    vec<Container> in_v;
    for (const str& in_path : in_path_v) {
      const Container& c = in_v.emplace_back(in_path);
      const sketch_config_t& first = in_v.front().get_config();
      if (!compatible_configs(c.get_config(), first))
        error_exit("Cannot merge sketches with different configurations: " + in_path);
      if (c.get_config().keep_seq != first.keep_seq)
        error_exit("Cannot merge sketches with different --keep-seq: " + in_path);
    }
    return in_v;
  }

  // A sketch is copied whole, so each of its sections keeps its offset relative to the sketch start.
  scentry relocate(const scentry& e, uint64_t pos)
  {
    scentry r;
    r.offset = pos;
    r.len = e.len;
    r.buckets_offset = pos + (e.buckets_offset - e.offset);
    r.buckets_len = e.buckets_len;
    if (e.windows_len) {
      r.windows_offset = pos + (e.windows_offset - e.offset);
      r.windows_len = e.windows_len;
    }
    return r;
  }

  void write_index(std::ostream& os, const vec<scentry>& index_v)
  {
    os.write(reinterpret_cast<const char*>(index_v.data()), static_cast<std::streamsize>(sizeof(scentry) * index_v.size()));
  }

  // Appends every sketch of `in` to `os`; `pos` is the output offset of the next sketch.
  void append_sketches(std::ostream& os, const Container& in, vec<scentry>& index_v, uint64_t& pos)
  {
    std::ifstream is(in.get_path(), std::ifstream::binary);
    check_fstream(is, "Cannot open container", in.get_path().string());
    vec<char> buf_v;
    for (const scentry& e : in.get_entries()) {
      buf_v.resize(e.len);
      is.seekg(static_cast<std::streamoff>(e.offset));
      is.read(buf_v.data(), static_cast<std::streamsize>(e.len));
      check_fstream(is, "Failed to read a sketch", in.get_path().string());
      os.write(buf_v.data(), static_cast<std::streamsize>(e.len));
      index_v.push_back(relocate(e, pos));
      pos += e.len;
    }
  }

  void write_merged(const fs::path& out_path, const vec<Container>& in_v, uint64_t nsketches)
  {
    std::ofstream os(out_path, std::ofstream::binary);
    check_fstream(os, "Cannot open output file", out_path.string());
    write_container_header(os, in_v.front().get_config(), nsketches);

    const std::streampos index_pos = os.tellp();
    write_index(os, vec<scentry>(nsketches)); // placeholder, rewritten once the offsets are known

    vec<scentry> index_v;
    index_v.reserve(nsketches);
    uint64_t pos = static_cast<uint64_t>(os.tellp());
    for (const Container& in : in_v)
      append_sketches(os, in, index_v, pos);

    os.seekp(index_pos);
    write_index(os, index_v);
    os.seekp(0, std::ios::end);
    check_fstream(os, "Failed to write the merged container", out_path.string());
  }

} // namespace

void MergeSC::merge()
{
  for (const str& in_path : in_path_v) {
    if (!check_input_type(in_path, input_type::container, "merge -i/--sketch-paths")) exit_failure();
  }
  reject_input_as_output(out_path, in_path_v);

  info_msg("Merging ", in_path_v.size(), " container(s)");
  const vec<Container> in_v = open_inputs(in_path_v);
  uint64_t nsketches = 0;
  for (const Container& in : in_v)
    nsketches += in.size();
  if (nsketches == 0) error_exit("Nothing to merge: the inputs hold no sketches");

  write_merged(out_path, in_v, nsketches);
  info_msg("Wrote ", nsketches, " sketch(es) to ", out_path.string());
}

MergeSC::MergeSC(CLI::App& sc)
{
  sc.add_option("-i,--sketch-paths", in_path_v, "Input containers to merge")->required()->check(CLI::ExistingFile);
  sc.add_option("-o,--output-path", out_path, "Path to store the merged container")->required();
}
