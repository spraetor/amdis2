#include "Timer.h"
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

namespace AMDiS
{
  Timer::Timer()
    : t0(Clock::now())
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    , t0_mpi(MPI::Wtime())
#endif
  { }

  void Timer::reset()
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    t0_mpi = MPI::Wtime();
#else
    t0 = Clock::now();
#endif
  }

  double Timer::elapsed() const
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    return MPI::Wtime() - t0_mpi;
#else
    auto t1 = Clock::now();
    fsec fs = t1 - t0;
    return fs.count();
#endif
  }

} // end namespace AMDiS
