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
#ifndef AMDIS_TIMER_H
#define AMDIS_TIMER_H

#include <chrono>

namespace AMDiS {
  
  /// Helper class to distinguish between different time measurement methods
  class Timer
  {
  public:
    /// initializes the timer with current time
    Timer();

    /// resets the timer to current time
    void reset();

    /// returns the elapsed time (from construction or last reset) to now in seconds
    double elapsed() const;
    
  private:
    /// begin value for sequentiell measurement
    std::chrono::time_point<std::chrono::system_clock> first_seq;

    /// begin value for parallel measurement
    double first_mpi;

  };
}

#endif // AMDIS_TIMER_H
