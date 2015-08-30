#include "Timer.h"
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

namespace AMDiS
{
  Timer::Timer()
    : first_seq(std::chrono::system_clock::now())
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    ,first_mpi(MPI::Wtime())
#endif
  { }

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
    auto elapsed_seconds = duration_cast<microseconds>(system_clock::now()-first_seq);
    return elapsed_seconds.count() / 1.e6;
#endif
  }

} // end namespace AMDiS
