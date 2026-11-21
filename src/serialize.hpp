#ifndef _SERIALIZE_HPP
#define _SERIALIZE_HPP

#include "msg.hpp"
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <ostream>
#include <type_traits>

inline constexpr uint64_t round_up_word(uint64_t x) noexcept { return (x + 7) & ~uint64_t(7); }

inline void pad_to_word(std::ostream& os)
{
  static const char zeros[8] = {};
  const uint64_t pos = static_cast<uint64_t>(os.tellp());
  const uint64_t pad = round_up_word(pos) - pos;
  if (pad) os.write(zeros, static_cast<std::streamsize>(pad));
}

inline void require_bytes(const char* ptr, const char* end, size_t nbytes, const char* field)
{
  if (ptr == nullptr || end == nullptr || static_cast<size_t>(end - ptr) < nbytes) {
    error_exit(concat_msg("Truncated sketch while reading ", field));
  }
}

template<typename T>
T read_trivial(const char*& ptr, const char* end, const char* field)
{
  static_assert(std::is_trivially_copyable_v<T>, "read_trivial memcpy's T");
  require_bytes(ptr, end, sizeof(T), field);
  T value{};
  std::memcpy(&value, ptr, sizeof(T));
  ptr += sizeof(T);
  return value;
}

template<typename T>
void write_trivial(std::ostream& os, const T& value)
{
  static_assert(std::is_trivially_copyable_v<T>, "write_trivial memcpy's T");
  os.write(reinterpret_cast<const char*>(&value), sizeof(T));
}

template<typename T>
const T* view_array(const char*& ptr, const char* end, uint64_t count, const char* field)
{
  const uint64_t nbytes = round_up_word(count * sizeof(T));
  require_bytes(ptr, end, nbytes, field);
  const T* array = count ? reinterpret_cast<const T*>(ptr) : nullptr;
  ptr += nbytes;
  return array;
}

template<typename T>
void write_array(std::ostream& os, const T* data, uint64_t count)
{
  if (count) os.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(count * sizeof(T)));
  pad_to_word(os);
}

#endif
