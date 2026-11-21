#ifndef _MSG_HPP
#define _MSG_HPP

#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <utility>

#define assertm(exp, msg) assert(((void)(msg), exp))

// --verbose: report progress even when stderr is not a TTY, and per-reference detail.
extern bool verbose;

inline bool stderr_is_tty() noexcept { return isatty(STDERR_FILENO) != 0; }

// Set while a progress line is on screen, i.e. between progress_continue() and progress_done().
inline bool& progress_line_open() noexcept
{
  static bool open = false;
  return open;
}

inline void close_progress_line()
{
  if (!progress_line_open()) return;
  std::cerr << '\n';
  progress_line_open() = false;
}

// Progress line "<label> <done>/<total>..."; silent unless stderr is a TTY or --verbose.
inline void progress_continue(const std::string& label, uint64_t done, uint64_t total)
{
  if (!stderr_is_tty() && !verbose) return;
  std::cerr << '\r' << label << ' ' << done << '/' << total << "..." << std::flush;
  progress_line_open() = true;
}

// Finish a progress line; pairs with progress_continue().
inline void progress_done()
{
  if (stderr_is_tty() || verbose) std::cerr << '\n';
  progress_line_open() = false;
}

inline std::string concat_msg() { return {}; }

template<typename... Args>
inline std::string concat_msg(Args&&... args)
{
  std::ostringstream oss;
  oss << std::boolalpha;
  (oss << ... << std::forward<Args>(args));
  return oss.str();
}

template<typename... Args>
inline void error_msg(Args&&... args)
{
  close_progress_line();
  std::cerr << "[ERROR] " << concat_msg(std::forward<Args>(args)...) << std::endl;
}

// Unrecoverable error; `main` reports it and exits. Safe to throw on worker threads.
struct fatal_error : std::runtime_error
{
  using std::runtime_error::runtime_error;
};

// For callers that have already reported the error with error_msg().
[[noreturn]] inline void exit_failure(int code = EXIT_FAILURE)
{
  close_progress_line();
  std::exit(code);
}

[[noreturn]] inline void error_exit(const std::string& msg)
{
  throw fatal_error(msg);
}

template<typename... Args>
inline void warn_msg(Args&&... args)
{
  close_progress_line();
  std::cerr << "[WARNING] " << concat_msg(std::forward<Args>(args)...) << std::endl;
}

template<typename... Args>
inline void info_msg(Args&&... args)
{
  close_progress_line();
  std::cerr << concat_msg(std::forward<Args>(args)...) << '\n';
}

template<typename... Args>
inline void verbose_msg(Args&&... args)
{
  if (verbose) info_msg(std::forward<Args>(args)...);
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
