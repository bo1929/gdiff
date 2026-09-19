#ifndef _SERIALIZE_HPP
#define _SERIALIZE_HPP

#include "msg.hpp"
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ostream>
#include <type_traits>

inline constexpr uint64_t round_up_word(uint64_t x) noexcept { return (x + 7) & ~uint64_t(7); }

// Pad a stream to the next 8-byte boundary.
inline void pad_to_word(std::ostream& os)
{
  static const char zeros[8] = {};
  const uint64_t pos = static_cast<uint64_t>(os.tellp());
  const uint64_t pad = round_up_word(pos) - pos;
  if (pad) os.write(zeros, static_cast<std::streamsize>(pad));
}

// Advance p by n bytes, failing if the range is short.
inline const char* need(const char* p, const char* end, size_t n, const char* what)
{
  if (p == nullptr || end == nullptr || static_cast<size_t>(end - p) < n) {
    error_exit(concat_msg("Truncated sketch while reading ", what));
  }
  return p;
}

// Read one fixed-size value by memcpy; T must be trivially copyable.
template<typename T>
T read_trivial(const char*& p, const char* end, const char* what)
{
  static_assert(std::is_trivially_copyable_v<T>, "read_trivial memcpy's T");
  p = need(p, end, sizeof(T), what);
  T v{};
  std::memcpy(&v, p, sizeof(T));
  p += sizeof(T);
  return v;
}

// Write one fixed-size value by memcpy; T must be trivially copyable.
template<typename T>
void write_trivial(std::ostream& os, const T& v)
{
  static_assert(std::is_trivially_copyable_v<T>, "write_trivial memcpy's T");
  os.write(reinterpret_cast<const char*>(&v), sizeof(T));
}

// View a length-prefixed array of T; sections are padded to a word.
template<typename T>
const T* view_array(const char*& p, const char* end, uint64_t n, const char* what)
{
  const uint64_t bytes = round_up_word(n * sizeof(T));
  p = need(p, end, bytes, what);
  const T* out = n ? reinterpret_cast<const T*>(p) : nullptr;
  p += bytes;
  return out;
}

template<typename T>
void write_array(std::ostream& os, const T* data, uint64_t n)
{
  if (n) os.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(n * sizeof(T)));
  pad_to_word(os);
}

#endif
