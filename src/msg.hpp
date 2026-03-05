#ifndef _MSG_HPP
#define _MSG_HPP

#include <ctime>
#include <chrono>
#include <fstream>
#include <iostream>
#include <cassert>

#define assertm(exp, msg) assert(((void)msg, exp))

[[noreturn]] inline void error_exit(const std::string& msg, int code = EXIT_FAILURE)
{
  auto ct = std::chrono::system_clock::now();
  std::time_t tt = std::chrono::system_clock::to_time_t(ct);
  std::cerr << "[ERROR] " << std::ctime(&tt) << ": " << msg << std::endl;
  std::exit(code);
}

inline void warn_msg(const std::string& msg)
{
  auto ct = std::chrono::system_clock::now();
  std::time_t tt = std::chrono::system_clock::to_time_t(ct);
  std::cerr << "[WARNING] " << std::ctime(&tt) << ": " << msg << std::endl;
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

template void check_fstream<std::ifstream>(const std::ifstream&, const std::string&, const std::string&);
template void check_fstream<std::ofstream>(const std::ofstream&, const std::string&, const std::string&);

#endif
