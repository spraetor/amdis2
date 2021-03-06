/** \file KrylovPreconditioner.h */

#pragma once

// AMDiS includes
#include "solver/ITL_Preconditioner.hpp"
#ifdef HAVE_PARALLEL_MTL4
#include "parallel/PITL_Solver.hpp"
#endif

namespace AMDiS
{
  /**
   * \ingroup Solver
   *
   * \brief Implements an inner solver for (flexible) krylov subspace solvers
   *
   * Any AMDiS solver can be used as preconditioner for any other AMDiS solver.
   * Initfile-syntax:<br />
   * [...]->solver: any-krylov-solver<br />
   * [...]->solver->left precon: [krylov|solver]<br />
   * [...]->solver->left precon->solver: any-amdis-solver<br />
   * [...]->solver->left precon->solver->parameter1: ...<br />
   * [...]->solver->left precon->solver->parameter2: ...<br />
   * [...]->solver->left precon->solver->parameter3: ...
   *
   * This Krylov-Preconditioner can also be used as preconditioner for the inner solver, so
   * that a sequence of iner-inner-solver can be loaded.
   **/
  template <class MatrixType, class VectorType>
  struct KrylovPreconditioner : public ITL_PreconditionerBase<MatrixType, VectorType>
  {
    using Self        = KrylovPreconditioner;
    using precon_base = ITL_PreconditionerBase<MatrixType, VectorType>;

    struct Creator : public CreatorInterfaceName<precon_base>
    {
      virtual precon_base* create() override
      {
        return new Self(this->name);
      }
    };

    KrylovPreconditioner(std::string name)
      : fullMatrix(NULL),
        solver(NULL),
        runner(NULL)
    {

#if defined HAVE_PARALLEL_PETSC
      std::string backend("p_petsc");
#elif defined HAVE_PARALLEL_MTL
      std::string backend("p_mtl");
#elif defined HAVE_PETSC || defined HAVE_SEQ_PETSC
      std::string backend("petsc");
#else
      std::string backend("mtl");
#endif

      // === read backend-name ===
      std::string initFileStr = name + "->solver";
      Parameters::get(initFileStr + "->backend", backend);

      // === read solver-name ===
      std::string solverType("0");
      Parameters::get(initFileStr, solverType);

      if (backend != "0" && backend != "no" && backend != "")
        solverType = backend + "_" + solverType;

      LinearSolverCreator* solverCreator =
        dynamic_cast<LinearSolverCreator*>(CreatorMap<LinearSolverInterface>::getCreator(solverType, initFileStr));
      TEST_EXIT(solverCreator)
      ("No valid solver type found in parameter \"%s\"\n", initFileStr.c_str());
      solverCreator->setName(initFileStr);
      solver = solverCreator->create();

      runner = dynamic_cast<RunnerBase<MatrixType, VectorType>*>(solver->getRunner());
    }

    /// Destructor.
    ~KrylovPreconditioner()
    {
      delete solver;
    }

    /// Implementation of \ref ITL_PreconditionerBase::init()
    virtual void init(SolverMatrix<Matrix<DOFMatrix*>> const& A,
                      MatrixType const& fullMatrix_) override
    {
      fullMatrix = &fullMatrix_;
      runner->init(A, fullMatrix_);
    }

    /// Implementation of \ref PreconditionerInterface::init()
    virtual void exit() override
    {
      runner->exit();
    }

    /// Implementation of \ref ITL_PreconditionerBase::solve()
    virtual void solve(VectorType const& b, VectorType& x) const override
    {
      initVector(x);
      runner->solve(*fullMatrix, x, b);
    }

    /// Implementation of \ref ITL_PreconditionerBase::adjoint_solve()
    virtual void adjoint_solve(VectorType const& b, VectorType& x) const override
    {
      initVector(x);
      runner->adjoint_solve(*fullMatrix, x, b);
    }

  protected: // methods

    template <class VectorT>
      Requires_t<mtl::traits::is_distributed<VectorT>>
    initVector(VectorT& x) const
    {
#ifdef HAVE_PARALLEL_MTL4
      x.init_distribution(col_distribution(*fullMatrix), num_cols(*fullMatrix));
#endif
      set_to_zero(x);
    }

    template <class VectorT>
      Requires_t<not_<mtl::traits::is_distributed<VectorT>>>
    initVector(VectorT& x) const
    {
      x.change_dim(num_cols(*fullMatrix));
      set_to_zero(x);
    }

  protected: // member variables
    MatrixType const* fullMatrix;

    LinearSolverInterface* solver;
    RunnerBase<MatrixType, VectorType>* runner;
  };


#ifdef HAVE_PARALLEL_MTL4
  using KrylovPreconditionerParallel = KrylovPreconditioner<MTLTypes::PMTLMatrix, 
							    MTLTypes::PMTLVector>;
#endif
  using KrylovPreconditionerSeq = KrylovPreconditioner<MTLTypes::MTLMatrix, 
						       MTLTypes::MTLVector>;


} // end namespace AMDiS
