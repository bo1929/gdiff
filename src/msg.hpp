#ifndef _MSG_HPP
#define _MSG_HPP

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>

#define assertm(exp, msg) assert(((void)(msg), exp))

inline std::string concat_msg() { return {}; }

template<typename... Args>
inline std::string concat_msg(Args&&... args)
{
  std::ostringstream oss;
  oss << std::boolalpha;
  (oss << ... << std::forward<Args>(args));
  return oss.str();
}

[[noreturn]] inline void error_exit(const std::string& msg, int code = EXIT_FAILURE)
{
  std::cerr << "[ERROR] " << msg << std::endl;
  std::exit(code);
}

inline void warn_msg(const std::string& msg)
{
  std::cerr << "[WARNING] " << msg << std::endl;
}

template<typename... Args>
inline void warn_pmsg(const std::string& prefix, Args&&... args)
{
  warn_msg(concat_msg(prefix, ": ", std::forward<Args>(args)...));
}

template<typename... Args>
inline void cerr_msg(Args&&... args)
{
  std::cerr << concat_msg(std::forward<Args>(args)...) << '\n';
}

template<typename stream_t>
inline void check_fstream(const stream_t& stream, const std::string& msg, const std::string& path)
{
  if (!stream.good()) {
    if (!path.empty()) {
      error_exit(msg + ": " + path);
    } else {
      error_exit(msg);
    }
  }
}

#endif
