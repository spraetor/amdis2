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



/** \file OpenMP.h */

#ifndef AMDIS_OPENMP_H
#define AMDIS_OPENMP_H

#ifdef _OPENMP
#include <omp.h>
#include <vector>
#endif

namespace AMDiS {
  
#ifdef _OPENMP

  template<typename T>
  class ThreadPrivate {
  public:
    ThreadPrivate()
      : data(omp_get_max_threads())
    {}

    ThreadPrivate(T val)
      : data(omp_get_max_threads(), val)
    {}

    inline T& get()
    {
#if (DEBUG != 0)
      if (omp_get_thread_num() >= data.size()) {
	std::cout << "Error in ThreadPrivate::get()!\n";
	exit(0);
      }
#endif
      return data[omp_get_thread_num()];
    }

    inline void set(T& val)
    {
#if (DEBUG != 0)
      if (omp_get_thread_num() >= data.size()) {
	std::cout << "Error in ThreadPrivate::set()!\n";
	exit(0);
      }
#endif
      data[omp_get_thread_num()] = val;
    }

  private:
    std::vector<T> data;
  };

#else

  template<typename T>
  class ThreadPrivate {
  public:
    ThreadPrivate() {}

    ThreadPrivate(T val)
      : data(val)
    {}

    inline T& get()
    {
      return data;
    }

    inline void set(T& val)
    {
      data = val;
    }

  private:
    T data;
  };

#endif

}

#endif
