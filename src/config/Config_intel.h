/** \file Config_intel.h */

#pragma once

// MMmm (M...major, m...minor)
#define INTEL_VERSION __INTEL_COMPILER

#define COMPILER_NAME "icc"
#define COMPILER_VERSION INTEL_VERSION

// alignement specification
// ------------------------
#include <malloc.h>
#define AMDIS_ALIGNED(type,name,N)  __declspec(align(CACHE_LINE)) type name[N]
typedef __declspec(align(CACHE_LINE)) double aligned_double;
typedef __declspec(align(CACHE_LINE)) float  aligned_float;
typedef __declspec(align(CACHE_LINE)) int    aligned_int;
typedef __declspec(align(CACHE_LINE)) size_t aligned_size_t;
#define AMDIS_ASSUME_ALIGNED(var)   var; __assume_aligned(var, CACHE_LINE)

#define AMDIS_ALIGNED_ALLOC(type,size) reinterpret_cast<type*>(_mm_malloc(size*sizeof(type),CACHE_LINE))
#define AMDIS_ALIGNED_FREE(addr) _mm_free(addr)

// some compiler attributes
// ------------------------
#define AMDIS_NOINLINE
#define AMDIS_ALWAYS_INLINE
#define AMDIS_OPENMODE       std::ios::openmode

// C++11 features
// --------------
#if __cplusplus > 199711L

#if INTEL_VERSION >= 1201
#define HAS_VARIADIC_TEMPLATES 1
#endif

#if INTEL_VERSION >= 1201
#define HAS_ALIAS_TEMPLATES 1
#endif

#if INTEL_VERSION >= 1200
#define HAS_DECLTYPE 1
#endif

#if INTEL_VERSION >= 1400
#define HAS_CONSTEXPR 1
#endif

#if INTEL_VERSION >= 1400
#define HAS_DELEGATING_CONSTRUCTORS 1
#endif

#if INTEL_VERSION >= 1400
#define HAS_RANGE_BASED_FOR 1
#endif

#if INTEL_VERSION >= 1400
#define HAS_INITIALIZER_LISTS 1
#endif

#if INTEL_VERSION >= 1400
#define HAS_OVERRIDE 1
#endif

#if INTEL_VERSION >= 1400
#define HAS_TYPED_ENUMS 1
#endif

#if INTEL_VERSION >= 1200
#define HAS_RVALUE_REFERENCES 1
#endif

#endif
