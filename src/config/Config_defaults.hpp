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
#define ALIGNED_ALLOC(type,size) (type*)aligned_malloc(size*sizeof(type), CACHE_LINE)
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