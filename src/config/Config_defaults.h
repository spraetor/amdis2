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



/** \file Config_defaults.h */

#pragma once

#ifndef COMPILER_NAME
  #define COMPILER_NAME "Unknown"
#endif
#ifndef COMPILER_VERSION
  #define COMPILER_VERSION 0
#endif

// alignement specification
// ------------------------
#ifndef ALIGNED
#define ALIGNED(type,name,N)  type name[N]
#define ASSUME_ALIGNED(var)   var

typedef double aligned_double;
typedef float  aligned_float;
typedef int    aligned_int;
typedef size_t aligned_size_t;
#endif

#ifndef ALIGNED_ALLOC
  // define aligned_malloc and aligned_free somewhere else, before using the macros
  #define ALIGNED_ALLOC(type,size,alignment) (type*)aligned_malloc(size*sizeof(type), alignment)
  #define ALIGNED_FREE(ptr) aligned_free(ptr)
#endif

// some compiler attributes
// ------------------------
#ifndef NOINLINE
  #define NOINLINE
#endif
#ifndef ALWAYS_INLINE
  #define ALWAYS_INLINE
#endif
#ifndef OPENMODE
  #define OPENMODE std::ios::openmode
#endif

// C++11 features
// --------------
#ifndef HAS_VARIADIC_TEMPLATES
  #define HAS_VARIADIC_TEMPLATES 0
#endif

#ifndef HAS_ALIAS_TEMPLATES
  #define HAS_ALIAS_TEMPLATES 0
#endif

#ifndef HAS_DECLTYPE
  #define HAS_DECLTYPE 0
#endif

#ifndef HAS_CONSTEXPR
  #define HAS_CONSTEXPR 0
#endif

#ifndef HAS_DELEGATING_CONSTRUCTORS
  #define HAS_DELEGATING_CONSTRUCTORS 0
#endif

#ifndef HAS_RANGE_BASED_FOR
  #define HAS_RANGE_BASED_FOR 0
#endif

#ifndef HAS_INITIALIZER_LISTS
  #define HAS_INITIALIZER_LISTS 0
#endif

#ifndef HAS_OVERRIDE
  #define HAS_OVERRIDE 0
  #define override
#endif

#ifndef HAS_TYPED_ENUMS
  #define HAS_TYPED_ENUMS 0
#endif

#ifndef HAS_RVALUE_REFERENCES
  #define HAS_RVALUE_REFERENCES 0
#endif

// #ifdef BOOST_NO_CXX11_NOEXCEPT
//   #define noexcept
// #endif