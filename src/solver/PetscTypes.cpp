#include "solver/PetscTypes.hpp"

namespace AMDiS
{

  void PetscMatrix::destroy()
  {
    if (assembled)
    {
      MatDestroy(&matrix);
      matrix = PETSC_NULL;
      for (size_t i = 0; i < nestMat.size(); i++)
      {
        if (nestMat[i] != PETSC_NULL)
          MatDestroy(&(nestMat[i]));
        nestMat[i] = PETSC_NULL;
      }
    }
    assembled = false;
  }


  void PetscVector::destroy()
  {
    if (assembled)
    {
      VecDestroy(&vector);
      vector = PETSC_NULL;
      for (size_t i = 0; i < nestVec.size(); i++)
      {
        VecDestroy(&(nestVec[i]));
        nestVec[i] = PETSC_NULL;
      }
    }
    assembled = false;
  }


  PetscParameters::PetscParameters()
  {
    matSolverPackage["superlu"] = true;
    matSolverPackage["superlu_dist"] = true;
    matSolverPackage["umfpack"] = true;
    matSolverPackage["cholmod"] = true;
    matSolverPackage["essl"] = true;
    matSolverPackage["lusol"] = true;
    matSolverPackage["mumps"] = true;
    matSolverPackage["pastix"] = true;
    matSolverPackage["matlab"] = true;
    matSolverPackage["bas"] = true;
    matSolverPackage["cusparse"] = true;
    matSolverPackage["bstrm"] = true;
    matSolverPackage["sbstrm"] = true;
    matSolverPackage["elemental"] = true;
    matSolverPackage["clique"] = true;
    matSolverPackage["direct"] = true;

    // mtl-name => petsc-name
    // (or petsc-name => petsc-name if solver does not exist for mtl)
    solverMap["cg"] = "cg";
    solverMap["richardson"] = "richardson";
    solverMap["chebyshev"] = "chebyshev";
    solverMap["groppcg"] = "groppcg";
    solverMap["pipecg"] = "pipecg";
    solverMap["cgne"] = "cgne";
    solverMap["nash"] = "nash";
    solverMap["stcg"] = "stcg";
    solverMap["gltr"] = "gltr";
    solverMap["gmres"] = "gmres";
    solverMap["fgmres"] = "fgmres";
    solverMap["lgmres"] = "lgmres";
    solverMap["dgmres"] = "dgmres";
    solverMap["pgmres"] = "pgmres";
    solverMap["tcqmr"] = "tcqmr";
    solverMap["bcgs"] = "bcgs";
    solverMap["bicgstab"] = "bcgs";
    solverMap["ibcgs"] = "ibcgs";
    solverMap["fbcgs"] = "fbcgs";
    solverMap["fbcgsr"] = "fbcgsr";
    solverMap["bcgsl"] = "bcgsl";
    solverMap["bicgstab_ell"] = "bcgsl";
    solverMap["bicgstab2"] = "bcgsl";
    solverMap["cgs"] = "cgs";
    solverMap["tfqmr"] = "tfqmr";
    solverMap["cr"] = "cr";
    solverMap["pipecr"] = "pipecr";
    solverMap["lsqr"] = "lsqr";
    solverMap["preonly"] = "preonly";
    solverMap["qcg"] = "qcg";
    solverMap["bicg"] = "bicg";
    solverMap["minres"] = "minres";
    solverMap["symmlq"] = "symmlq";
    solverMap["lcd"] = "lcd";
    solverMap["python"] = "python";
    solverMap["gcr"] = "gcr";
    solverMap["specest"] = "specest";

    solverMap["superlu"] = "superlu";
    solverMap["superlu_dist"] = "superlu_dist";
    solverMap["umfpack"] = "umfpack";
    solverMap["mumps"] = "mumps";
    solverMap["direct"] = "mumps";

    preconMap["diag"] = "jacobi";
    preconMap["jacobi"] = "jacobi";
    preconMap["ilu"] = "ilu";
    preconMap["ic"] = "icc";

    emptyParam[""] = true;
    emptyParam["0"] = true;
    emptyParam["no"] = true;
    emptyParam["none"] = true;
  }
} // end namespace AMDiS
