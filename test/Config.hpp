#pragma once

/** \brief current AMDiS version */
#ifndef AMDIS_VERSION
#define AMDIS_VERSION  "AMDiS: Version 2.0.0"
#endif

//#include <boost/config.hpp>

#define CACHE_LINE 16

// if AMDIS_FIXED_SIZE == 1 use static arrays
#ifndef AMDIS_FIXED_SIZE
#define AMDIS_FIXED_SIZE 0
//  #define DOW 2
//  #define DIM 2
#endif

#if defined(__clang__)					// Clang/LLVM.
#include "config/Config_clang.hpp"

#elif defined(__ICC) || defined(__INTEL_COMPILER)	// Intel ICC/ICPC.
#include "config/Config_intel.hpp"

#elif defined(__GNUC__) || defined(__GNUG__)		// GNU GCC/G++.
#include "config/Config_gcc.hpp"

#elif defined(__HP_cc) || defined(__HP_aCC)
error:
not supported compiler

#elif defined(__IBMC__) || defined(__IBMCPP__)
error:
not supported compiler

#elif defined(_MSC_VER)					// Microsoft Visual Studio.
#include "config/Config_msc.hpp"

#elif defined(__PGI)					// Portland Group PGCC/PGCPP.
error:
not supported compiler
//   #include "Config_pgi.h"

#elif defined(__SUNPRO_C) || defined(__SUNPRO_CC)
error:
not supported compiler
#endif

#include "config/Config_defaults.hpp"

typedef unsigned short small_t;   // only allow small matrices

// some workarounds for mtl (since the namespace has changed)
#define MTL_VEC mtl::vector
#define MTL_MAT mtl::matrix
