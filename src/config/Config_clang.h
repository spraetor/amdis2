/** \file Config_clang.h */

#pragma once

#define CLANG_VERSION (__clang_major__ * 10000 + __clang_minor__ * 100 + __clang_patchlevel__)

#define COMPILER_NAME "clang"
#define COMPILER_VERSION CLANG_VERSION

// alignement specification
// ------------------------
#define ALIGNED(type,name,N)  type name[N] __attribute__ ((aligned(CACHE_LINE)))
#define ASSUME_ALIGNED(var)   __builtin_assume_aligned(var, CACHE_LINE)

typedef double aligned_double   __attribute__ ((aligned(CACHE_LINE)));
typedef float  aligned_float    __attribute__ ((aligned(CACHE_LINE)));
typedef int    aligned_int      __attribute__ ((aligned(CACHE_LINE)));
typedef size_t aligned_size_t   __attribute__ ((aligned(CACHE_LINE)));

// some compiler attributes
// ------------------------
#define NOINLINE                __attribute__ ((noinline))
#define ALWAYS_INLINE           __attribute__ ((always_inline))
#define OPENMODE                std::ios::openmode

// C++11 features
// --------------
#if __cplusplus > 199711L

// __has_feature(cxx_rvalue_references)
#if CLANG_VERSION >= 20900
  #define HAS_VARIADIC_TEMPLATES 1
#endif

#if CLANG_VERSION >= 30000
  #define HAS_ALIAS_TEMPLATES 1
#endif

#if CLANG_VERSION >= 20900
  #define HAS_DECLTYPE 1
#endif

#if CLANG_VERSION >= 30100
  #define HAS_CONSTEXPR 1
#endif

#if CLANG_VERSION >= 30000
  #define HAS_DELEGATING_CONSTRUCTORS 1
#endif

#if CLANG_VERSION >= 30000
  #define HAS_RANGE_BASED_FOR 1
#endif

#if CLANG_VERSION >= 30100
  #define HAS_INITIALIZER_LISTS 1
#endif

#if CLANG_VERSION >= 30000
  #define HAS_OVERRIDE 1
#endif

#if CLANG_VERSION >= 20900
  #define HAS_TYPED_ENUMS 1
#endif

#if CLANG_VERSION >= 20900
  #define HAS_RVALUE_REFERENCES 1
#endif

#endif
