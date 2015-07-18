/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/



/** \file Config_intel.h */

#pragma once

// MMmm (M...major, m...minor)
#define INTEL_VERSION __INTEL_COMPILER

#define COMPILER_NAME "icc"
#define COMPILER_VERSION INTEL_VERSION

// alignement specification
// ------------------------
#define ALIGNED(type,name,N)  __declspec(align(CACHE_LINE)) type name[N]
#define ASSUME_ALIGNED(var)   var; __assume_aligned(var, CACHE_LINE)

typedef __declspec(align(CACHE_LINE)) double aligned_double;
typedef __declspec(align(CACHE_LINE)) float  aligned_float;
typedef __declspec(align(CACHE_LINE)) int    aligned_int;
typedef __declspec(align(CACHE_LINE)) size_t aligned_size_t;

#define ALIGNED_ALLOC(type,size,alignment) (type*)_mm_malloc(size,alignment)
#define ALIGNED_FREE(addr) _mm_free(addr)

// some compiler attributes
// ------------------------
#define NOINLINE
#define ALWAYS_INLINE
#define OPENMODE       std::ios::openmode

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
