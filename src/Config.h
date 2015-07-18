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



/** \file Config.h */

#pragma once

/** \brief current AMDiS version */
#ifndef AMDIS_VERSION
#define AMDIS_VERSION  "AMDiS: Version 0.9.1"
#endif

#include <boost/config.hpp>

#define CACHE_LINE 16

#if defined(__clang__)					// Clang/LLVM.
  #include "config/Config_clang.h"
  
#elif defined(__ICC) || defined(__INTEL_COMPILER)	// Intel ICC/ICPC. 
  #include "config/Config_intel.h"
  
#elif defined(__GNUC__) || defined(__GNUG__)		// GNU GCC/G++.
  #include "config/Config_gcc.h"
  
#elif defined(__HP_cc) || defined(__HP_aCC)
  error: not supported compiler
  
#elif defined(__IBMC__) || defined(__IBMCPP__)
  error: not supported compiler
  
#elif defined(_MSC_VER)					// Microsoft Visual Studio. 
  #include "config/Config_msc.h"
  
#elif defined(__PGI)					// Portland Group PGCC/PGCPP.
  error: not supported compiler
//   #include "Config_pgi.h"  

#elif defined(__SUNPRO_C) || defined(__SUNPRO_CC)
  error: not supported compiler
#endif

#include "config/Config_defaults.h"