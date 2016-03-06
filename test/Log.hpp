/** \file Log.h */

#pragma once

// std c++ headers: std::ostream, std::ofstream
#include <iostream>
#include <sstream>

#define STATIC_ASSERT(...) \
  static_assert(__VA_ARGS__, #__VA_ARGS__)

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


  /// prints a message
#define MSG(text) std::cout << text << "\n"

  /// prints an error message
#define ERROR(text) MSG("Error: " << text)

  /// prints an error message and exits
#define ERROR_EXIT(text) do { std::stringstream ss; ss << "Error: " << text; throw std::runtime_error(ss.str()); } while(0)

  /// prints a warning
#define WARNING(text) MSG("Warning: " << text)

  /// if test is false, an error message is printed
#define TEST(test, msg) if (!(test)) { ERROR(msg); }

  /// if test is false, an error message is printed and the program exits
#define TEST_EXIT(test, msg) if (!(test)) { ERROR_EXIT(msg); }

  /// In debug mode, it corresponds to ERROR_EXIT, otherwise it is noop.
#ifdef NDEBUG
  #define TEST_EXIT_DBG(test, msg)
  #define MSG_DBG(text)
  #define DBG_VAR(var)
#else
  #define TEST_EXIT_DBG(test, msg) TEST_EXIT(test, msg)
  #define MSG_DBG(text) MSG(text)
  #define DBG_VAR(var) var
#endif
