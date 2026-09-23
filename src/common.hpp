#ifndef _COMMON_HPP
#define _COMMON_HPP

#include "CLI11.hpp"
#include "msg.hpp"
#include "types.hpp"
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

inline constexpr const char* gdiff_version = "v0.2.1-rc";

inline str invocation;

inline void write_provenance(std::ostream& os)
{
  const std::time_t now = std::time(nullptr);
  str stamp = std::ctime(&now);
  if (!stamp.empty() && stamp.back() == '\n') stamp.pop_back();
  os << "# invocation: " << invocation << '\n';
  os << "# version: gdiff " << gdiff_version << ' ' << stamp << '\n';
}

inline void set_precision(std::ostream& stream, int precision) { stream << std::setprecision(precision); }

inline bool parse_input_entry(const str& line, str& name, std::filesystem::path& path)
{
  const size_t first = line.find_first_not_of(" \t\r");
  if (first == str::npos || line[first] == '#') return false;
  const size_t last = line.find_last_not_of(" \t\r");
  const str entry_line = line.substr(first, last - first + 1);

  name.clear();
  str entry = entry_line;
  const size_t tab = entry_line.find('\t');
  if (tab != str::npos) {
    const str raw_name = entry_line.substr(0, tab);
    const size_t name_end = raw_name.find_last_not_of(" \t");
    name = (name_end == str::npos) ? str{} : raw_name.substr(0, name_end + 1);
    const size_t path_begin = entry_line.find_first_not_of(" \t", tab + 1);
    entry = (path_begin == str::npos) ? str{} : entry_line.substr(path_begin);
  }
  if (entry.empty()) return false;
  path = std::filesystem::path(entry);
  return true;
}

// Open `path` for output and point `stream` at it; warn and no-op when `path` is empty.
inline void open_output(std::ofstream& file, const std::filesystem::path& path, std::ostream*& stream)
{
  if (path.empty()) {
    warn_msg("Empty output path; writing to stdout");
    return;
  }
  file.open(path);
  check_fstream(file, "Cannot open output file", path.string());
  stream = &file;
}

static const std::regex urlexp = std::regex(
  R"(^(?:(?:https?|ftp)://)(?:\S+@)?(?:(?!10(?:\.\d{1,3}){3})(?!127(?:\.\d{1,3}){3})(?!169\.254(?:\.\d{1,3}){2})(?!192\.168(?:\.\d{1,3}){2})(?!172\.(?:1[6-9]|2\d|3[0-1])(?:\.\d{1,3}){2})(?:[1-9]\d?|1\d\d|2[01]\d|22[0-3])(?:\.(?:1?\d{1,2}|2[0-4]\d|25[0-5])){2}(?:\.(?:[1-9]\d?|1\d\d|2[0-4]\d|25[0-4]))|(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+(?:\.(?:[a-z\u00a1-\uffff0-9]+-)*[a-z\u00a1-\uffff0-9]+)*(?:\.(?:[a-z\u00a1-\uffff]{2,})))(?::\d{2,5})?(?:/\S*)?$)");

inline bool match_url(const std::string& input) { return std::regex_match(input, urlexp); }

inline const auto url_validator = CLI::Validator(
  [](std::string& input) {
    if (match_url(input)) {
      return std::string("");
    } else {
      return "Given URL is not valid: " + input;
    }
  },
  "URL",
  "URL validator");

#endif
