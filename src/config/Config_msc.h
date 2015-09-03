/** \file Config_msc.h */

#pragma once

#define MSC_VERSION _MSC_VER
// MSC_VERSION == 1800 (Visual Studio 2013)
// MSC_VERSION == 1700 (Visual Studio 2012)
// MSC_VERSION == 1600 (Visual Studio 2010)
// MSC_VERSION == 1500 (Visual Studio 2008)

#define COMPILER_NAME "msc"
#define COMPILER_VERSION MSC_VERSION

// alignement specification
// ------------------------
#define ALIGNED(type,name,N)  __declspec(align(CACHE_LINE)) type name[N]
typedef __declspec(align(CACHE_LINE)) double aligned_double;
typedef __declspec(align(CACHE_LINE)) float  aligned_float;
typedef __declspec(align(CACHE_LINE)) int    aligned_int;
typedef __declspec(align(CACHE_LINE)) size_t aligned_size_t;

#include <malloc.h>
#define ALIGNED_ALLOC(type,size) reinterpret_cast<type*>(_aligned_malloc(size*sizeof(type),CACHE_LINE))
#define ALIGNED_FREE(ptr) _aligned_free(ptr);

// some compiler attributes
// ------------------------
#define NOINLINE         __declspec(noinline)
#define ALWAYS_INLINE    __forceinline
#define OPENMODE         std::ios::open_mode

// C++11 features
// --------------
#if __cplusplus > 199711L

#if MSC_VERSION >= 1800
#define HAS_VARIADIC_TEMPLATES 1
#endif

#if MSC_VERSION >= 1800
#define HAS_ALIAS_TEMPLATES 1
#endif

#if MSC_VERSION >= 1600
#define HAS_DECLTYPE 1
#endif

// #if MSC_VERSION >= 2000 (?)
#define HAS_CONSTEXPR 0
// #endif

#if MSC_VERSION >= 1800
#define HAS_DELEGATING_CONSTRUCTORS 1
#endif

#if MSC_VERSION >= 1700
#define HAS_RANGE_BASED_FOR 1
#endif

#if MSC_VERSION >= 1800
#define HAS_INITIALIZER_LISTS 1
#endif

#if MSC_VERSION >= 1700
#define HAS_OVERRIDE 1
#endif

#if MSC_VERSION >= 1700
#define HAS_TYPED_ENUMS 1
#endif

#if MSC_VERSION >= 1600
#define HAS_RVALUE_REFERENCES 1
#endif

#endif
