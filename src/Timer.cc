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
#include "Timer.h"
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

namespace AMDiS {
  
  Timer::Timer() : 
     first_seq(std::chrono::system_clock::now())
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    ,first_mpi(MPI::Wtime())
#endif
  {}

  void Timer::reset()
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    first_mpi = MPI::Wtime();
#else
    first_seq = std::chrono::system_clock::now();
#endif
  }

  double Timer::elapsed() const
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    return MPI::Wtime() - first_mpi;
#else
    using namespace std::chrono;
    auto elapsed_seconds = duration_cast<seconds>(system_clock::now()-first_seq);
    return static_cast<double>(elapsed_seconds.count());
#endif
  }
  
} // end namespace AMDiS
