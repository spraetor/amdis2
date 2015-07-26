/** \file Traverse.h */

#pragma once

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

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    /// begin value for parallel measurement
    double first_mpi;
#endif
  };
  
}  // end namespace AMDiS
