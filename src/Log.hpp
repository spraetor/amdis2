/** \file Log.h */

#pragma once

// std c++ headers: std::ostream, std::ofstream
#include <ostream>

// AMDiS headers: ThreadPrivate
#include "OpenMP.hpp"


#define STATIC_ASSERT(...) \
  static_assert(__VA_ARGS__, #__VA_ARGS__)

namespace AMDiS
{
  extern const char* funcName;

  void waitSec(int seconds);

  void processMemUsage(double& vm_usage, double& resident_set, bool inMegaByte = true);

  /** \brief
   * Manages the output of messages, warnings, errors, ...
   * Used by the macros FUNCNAME, ERROR, ERROR_EXIT, WARNING, TEST, MSG, INFO,
   * PRINT_INFO, WAIT, WAIT_REALLY.
   * Don't use this class directly but only via these macros!
   */
  struct Msg
  {
    /// Prints a formated message to the message stream
    static void print(const char* format, ...);

    /// Prints a formated message to the error stream
    static void print_error(const char* format, ...);

    /// Prints a formated message to the error stream and exits
    static void print_error_exit(const char* format, ...);

    ///
    static void catch_error_exit(const char* /* format */, ...) {}

    /// Prints an error message with funcname, file, and line to the error stream
    static void print_error_funcname(const char* funcname,
                                     const char* file,
                                     int line);

    /// Prints a warning to the message stream
    static void print_warn(const char* format, ...);

    /// Prints a warning with funcname, file, and line to the message stream
    static void print_warn_funcname(const char* funcname,
                                    const char* file,
                                    int line);

    /// Prints the funcname to the message stream
    static void print_funcname(const char* funcname);

    /// Changes the message stream
    // static void change_out(std::ostream*);

    /// Changes the error stream
    // static void change_error_out(std::ofstream *fp);

    /// Creates a filestream and sets the error stream to this filestream
    // static void open_error_file(const char *filename, OPENMODE);

    /// Sets \ref msgInfo
    static void setMsgInfo(int info)
    {
      msgInfo = info;
    }

    /// Returns \ref msgInfo
    static int  getMsgInfo()
    {
      return msgInfo;
    }

    /// Sets \ref msgWait
    static void setMsgWait(bool wait)
    {
      msgWait = wait;
    }

    /// Returns \ref msgWait
    static bool getMsgWait()
    {
      return msgWait;
    }

    /// Waits for enter if w is true
    static void wait(bool w);

    /// Returns \ref out
    static std::ostream* getOutStream()
    {
      return out;
    }

    /// Returns \ref error
    static std::ostream* getErrorStream()
    {
      return error;
    }

#if HAVE_PARALLEL_DOMAIN_AMDIS
    /// In parallel computations, when this variable is true, only the 0 rank will
    /// print messages to the output stream. Error messages and warnings are always
    /// printed from all ranks.
    static bool outputMainRank;
#endif

  protected:
    /// Message stram
    static std::ostream* out;

    /// Error stream
    static std::ostream* error;

    /// Remember funcName to avoid multiple output of funcName within the same
    /// function call
    static ThreadPrivate<const char*> oldFuncName;

    /// Global info level
    static int msgInfo;

    /// Spezifies whether to wait when WAIT is called
    static bool msgWait;
  };

  // ===========================================================================
  // ===== message macros ======================================================
  // ===========================================================================

  /// Should be the first call in every functions. It defines the current
  /// function name nn for message output via MSG, WARNING, ...
#define FUNCNAME(nn) const char *funcName; funcName = nn;
  
#ifdef NDEBUG
  #define FUNCNAME_DBG(nn)
#else
#define FUNCNAME_DBG(nn) const char *funcName; funcName = nn;
#endif

  /// prints an error message
#define ERROR Msg::print_error_funcname(funcName,__FILE__, __LINE__),	\
  Msg::print_error

  /// prints an error message and exits
#define ERROR_EXIT Msg::print_error_funcname(funcName,__FILE__, __LINE__), \
  Msg::print_error_exit

  /// prints a warning
#define WARNING Msg::print_warn_funcname(funcName,__FILE__, __LINE__),	\
  Msg::print_warn

  /// if test is false, an error message is printed
#define TEST(test) if ((test));else ERROR

  /// if test is false, an error message is printed and the program exits
#define TEST_EXIT(test) if ((test));else ERROR_EXIT

  /// In debug mode, it corresponds to ERROR_EXIT, otherwise it is noop.
#ifdef NDEBUG
  #define TEST_EXIT_DBG(test) if (false) Msg::catch_error_exit
  #define DBG_VAR(var)
#else
#define TEST_EXIT_DBG(test) if ((test));else ERROR_EXIT
#define DBG_VAR(var) var
#endif

  /// prints a message
#define MSG Msg::print_funcname(funcName), Msg::print

#ifdef NDEBUG
  #define MSG_DBG
#else
#define MSG_DBG Msg::print_funcname(funcName), Msg::print
#endif

  /// prints a message, if min(Msg::msgInfo, info) >= noinfo
#define INFO(info,noinfo)						\
  if (Msg::getMsgInfo() && (std::min(Msg::getMsgInfo(), (info)) >= (noinfo))) MSG

  /// prints a message, if min(Msg::msgInfo, info) >= noinfo
#define PRINT_INFO(info,noinfo)						\
  if (Msg::getMsgInfo() && (std::min(Msg::getMsgInfo(), (info)) >= (noinfo))) Msg::print


  /// If the value of Msg::wait is not zero the macro will produce the message
  /// 'wait for <enter> ...' and will continue after pressing the return or enter
  /// key. Otherwise the program continues without a message.
#define WAIT Msg::wait(Msg::getMsgWait())

  /// produces the message 'wait for <enter> ...' and will continue after
  /// pressing the return or enter key.
#define WAIT_REALLY Msg::wait(true)

} // end namespace AMDiS
