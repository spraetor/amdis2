/** \file Traverse.h */

#pragma once

#include <chrono>

namespace AMDiS
{
  // TODO: is it necessary to use MPI::Wtime(), or is chrono just fine.

  /// Helper class to distinguish between different time measurement methods
  class Timer
  {
    using Clock     = std::chrono::high_resolution_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    using fsec      = std::chrono::duration<double>;
    
  public:
    /// initializes the timer with current time
    Timer();

    /// resets the timer to current time
    void reset();

    /// returns the elapsed time (from construction or last reset) to now in seconds
    double elapsed() const;

  private:
    /// begin value for sequentiell measurement
    TimePoint t0;

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    /// begin value for parallel measurement
    double t0_mpi;
#endif
  };

}  // end namespace AMDiS
