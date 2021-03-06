/** \file BlockPreconditioner.h */

#pragma once

#include "solver/ITL_Preconditioner.hpp"
#include "solver/SolverMatrix.hpp"

namespace AMDiS
{
  /// Basis preconditioner structure for block-preconditioners
  template <class MatrixType, class VectorType = MTLTypes::MTLVector>
  struct BlockPreconditioner : public ITL_PreconditionerBase<MatrixType, VectorType>
  {
    using Self        = BlockPreconditioner;
    using Super       = ITL_PreconditionerBase<MatrixType, VectorType>;
    using precon_base = Super;

    BlockPreconditioner()
      : A(NULL), 
	fullMatrix(NULL)
    {}

    /// Implementation of \ref ITL_PreconditionerBase::init
    /// Extract iranges from BlockMatrix to be used to extract sub-vectors and sub-matrices.
    virtual void init(SolverMatrix<Matrix<DOFMatrix*>> const& A_,
                      MatrixType const& fullMatrix_) override
    {
      A = &A_;
      fullMatrix = &fullMatrix_;

      BlockMapper mapper(A_);
      rows.resize(mapper.getNumComponents());
      int start = 0;
      for (int i = 0; i < mapper.getNumComponents(); i++)
      {
        mapper.setRow(i+1);
        int finish = mapper.row(0);
        rows[i].set(start, finish);
        start = finish;
      }
    }

    /// Implementation of \ref ITL_PreconditionerBase::solve
    virtual void solve(VectorType const& b, VectorType& x) const override
    {
      FUNCNAME("BlockPreconditioner::solve()");
      TEST_EXIT(false)("Must be implemented in derived classes!\n");
    }

    /// Implementation of \ref ITL_PreconditionerBase::adjoint_solve
    virtual void adjoint_solve(VectorType const& x, VectorType& y) const override
    {
      FUNCNAME("BlockPreconditioner::adjoint_solve()");
      TEST_EXIT(false)("Must be implemented in derived classes!\n");
    }

    mtl::irange const& getRowRange(size_t i) const
    {
      return rows[i];
    }

    mtl::irange const& getColRange(size_t i) const
    {
      return rows[i];
    }

    DOFMatrix::base_matrix_type const& getSubMatrix(size_t i, size_t j) const
    {
      return A->getSubMatrix(i, j);
    }

    template <class SolverType, class RunnerType>
    static void createSubSolver(std::string param, 
				SolverType*& solver, 
				RunnerType*& runner,
                                std::string solverType = "0", 
				std::string preconType = "no",
                                int max_iter = 100, 
				double tol = 1.e-8)
    {
      // definition of standard-backends
#if defined HAVE_PARALLEL_PETSC
      std::string backend("p_petsc");
#elif defined HAVE_PARALLEL_MTL
      std::string backend("p_mtl");
#elif defined HAVE_PETSC
      std::string backend("petsc");
#else
      std::string backend("mtl");
#endif

      // === read backend-name ===
      std::string initFileStr = param + "->solver";

      // === read solver-name ===
      std::string tmp = "0";
      Parameters::get(initFileStr, tmp);
      if (tmp == "0")
      {
        TEST_EXIT(solverType != "0")("You have to specify the parameter %s in the initfile!\n", initFileStr.c_str());

        Parameters::set(initFileStr, solverType);
        Parameters::set(initFileStr + "->backend", backend);
        Parameters::set(initFileStr + "->left precon", preconType);
        Parameters::set(initFileStr + "->max iteration", max_iter);
        Parameters::set(initFileStr + "->tolerance", tol);
      }
      else
      {
        solverType = tmp;
      }

      Parameters::get(initFileStr + "->backend", backend);

      if (backend != "0" && backend != "no" && backend != "")
        solverType = backend + "_" + solverType;

      LinearSolverCreator* solverCreator =
        dynamic_cast<LinearSolverCreator*>(CreatorMap<LinearSolverInterface>::getCreator(solverType, initFileStr));
      TEST_EXIT(solverCreator)
      ("No valid solver type found in parameter \"%s\"\n", initFileStr.c_str());
      solverCreator->setName(initFileStr);
      solver = dynamic_cast<SolverType*>(solverCreator->create());
      assert(solver != NULL);

      runner = dynamic_cast<RunnerType*>(solver->getRunner());
      assert(runner != NULL);
    }

  protected:
    SolverMatrix<Matrix<DOFMatrix*>> const* A;
    MatrixType const* fullMatrix;

    std::vector<mtl::irange> rows;
  };

} // end namespace AMDiS
