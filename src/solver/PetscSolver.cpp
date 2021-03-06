#ifdef HAVE_SEQ_PETSC

#define AMDIS_NO_EXTERN_PETSC_SOLVER
#include "solver/PetscSolver.hpp"
#undef AMDIS_NO_EXTERN_PETSC_SOLVER

// AMDiS includes
#include "solver/PetscTypes.hpp"

namespace AMDiS
{
  // explicit template instantiation
  template class PetscRunner<PetscMatrix, PetscVector>;
  template class PetscRunner<PetscMatrixNested, PetscVectorNested>;

} // end namespace AMDiS

#endif // HAVE_SEQ_PETSC
