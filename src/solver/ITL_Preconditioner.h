/** \file ITL_Preconditioner.h */

#pragma once

#include <solver/LinearSolverInterface.h>
#include <MTL4Types.h>
#include <DOFMatrix.h>
#include <CreatorInterface.h>

#include <boost/numeric/itl/itl.hpp>
#include <boost/numeric/itl/pc/ilu_0.hpp>
#include <boost/numeric/itl/pc/ic_0.hpp>
#include <boost/numeric/mtl/vector/assigner.hpp>

#include "itl/masslumping.hpp"

namespace AMDiS
{
  /**
   * \ingroup Solver
   *
   * \brief Common base class for wrappers to use ITL preconditioners in AMDiS.
   */
  template <class MatrixType, class VectorType>
  struct ITL_PreconditionerBase : public PreconditionerInterface
  {
    virtual void init(SolverMatrix<Matrix<DOFMatrix*>> const& A,
                      MatrixType const& fullMatrix) = 0;
		      
    virtual void exit() {}

    virtual void solve(VectorType const& x, VectorType& y) const = 0;

    virtual void adjoint_solve(VectorType const& x, VectorType& y) const = 0;
  };


  template <class MatrixType, class VectorType>
  itl::pc::solver<ITL_PreconditionerBase<MatrixType, VectorType>, VectorType, false>
  solve(ITL_PreconditionerBase<MatrixType, VectorType> const& P, VectorType const& vin)
  {
    return {P, vin};
  }

  template <class MatrixType, class VectorType>
  itl::pc::solver<ITL_PreconditionerBase<MatrixType, VectorType>, VectorType, true>
  adjoint_solve(ITL_PreconditionerBase<MatrixType, VectorType> const& P, VectorType const& vin)
  {
    return {P, vin};
  }


  /**
   * \ingroup Solver
   *
   * \brief Wrapper for using ITL preconditioners in AMDiS.
   */
  template <class Preconditioner, class MatrixType, class VectorType>
  class ITL_Preconditioner : public ITL_PreconditionerBase<MatrixType, VectorType>
  {
  public:
    using Self        = ITL_Preconditioner;
    using precon_base = ITL_PreconditionerBase<MatrixType, VectorType>;

    /// Creator class
    struct Creator : public CreatorInterfaceName<precon_base>
    {
      virtual precon_base* create() override
      {
        return new Self();
      }
    };

    /// Constructor.
//     ITL_Preconditioner() = default;

    /// Destructor.
    ~ITL_Preconditioner()
    {
      delete precon;
      precon = NULL;
    }

    /// Implementation of \ref ITL_PreconditionerBase::init()
    virtual void init(SolverMatrix<Matrix<DOFMatrix*>> const& A,
                      MatrixType const& fullMatrix) override
    {
      delete precon;
      precon = new Preconditioner(fullMatrix);
    }

    /// Implementation of \ref PreconditionerInterface::exit()
    virtual void exit() override
    {
      delete precon;
      precon = NULL;
    }

    /// Implementation of \ref ITL_PreconditionerBase::solve()
    virtual void solve(VectorType const& vin, VectorType& vout) const override
    {
      TEST_EXIT_DBG(precon)("No preconditioner initialized!\n");
      precon->solve(vin, vout);
    }

    /// Implementation of \ref ITL_PreconditionerBase::adjoint_solve()
    virtual void adjoint_solve(VectorType const& vin, VectorType& vout) const override
    {
      TEST_EXIT_DBG(precon)("No preconditioner initialized!\n");
      precon->adjoint_solve(vin, vout);
    }

  private:
    Preconditioner* precon = NULL;
  };


  /**
   * \ingroup Solver
   * \class AMDiS::DiagonalPreconditioner
   * \brief ITL_Preconditioner implementation of diagonal (jacobi) preconditioner,
   * \implements ITL_Preconditioner
   *
   * Diagonal preconditioner \f$ M^{-1} \f$ for the system \f$ Ax=b \f$ is defined as: \f$ M=diag(A) \f$.
   */
  using DiagonalPreconditioner = 
    ITL_Preconditioner<itl::pc::diagonal<MTLTypes::MTLMatrix>, 
		       MTLTypes::MTLMatrix, MTLTypes::MTLVector>;

  /**
   * \ingroup Solver
   * \class AMDiS::DiagonalPreconditioner
   * \brief ITL_Preconditioner implementation of diagonal (jacobi) preconditioner, 
   * \implements ITL_Preconditioner
   *
   * Diagonal preconditioner \f$ M^{-1} \f$ for the system \f$ Ax=b \f$ is defined as: \f$ M_ii=sum_j(A_ij) \f$.
   */
  using MassLumpingPreconditioner = 
    ITL_Preconditioner<itl::pc::masslumping<MTLTypes::MTLMatrix>, 
		       MTLTypes::MTLMatrix, MTLTypes::MTLVector>;


  /**
   * \ingroup Solver
   * \class AMDiS::IdentityPreconditioner
   * \brief ITL_Preconditioner implementation of identity preconditioner,
   * \implements ITL_Preconditioner
   *
   * Identity preconditioner. Behaves like no preconditioning.
   */
  using IdentityPreconditioner = 
    ITL_Preconditioner<itl::pc::identity<MTLTypes::MTLMatrix>, 
		       MTLTypes::MTLMatrix, MTLTypes::MTLVector>;


  /**
   * \ingroup Solver
   * \class AMDiS::ILUPreconditioner
   * \brief ITL_Preconditioner implementation of ILU (Incomplete LU factorization) 
   * preconditioner. \implements ITL_Preconditioner
   *
   * The preconditioner is used from ITL. It corresponds for instance to
   * "Iterative Methods for Sparce Linear Systems", second edition, Yousef Saad.
   *  The preconditioner is described in chapter 10.3 (algorithm 10.4).
   */
  using ILUPreconditioner = 
    ITL_Preconditioner<itl::pc::ilu_0<MTLTypes::MTLMatrix>, 
		       MTLTypes::MTLMatrix, MTLTypes::MTLVector>;


  /**
   * \ingroup Solver
   * \class AMDiS::ICPreconditioner
   * \brief ITL_Preconditioner implementation of IC (Incomplete Cholesky factorization) 
   * preconditioner. \implements ITL_Preconditioner
   *
   * IC (Incomplete Cholesky factorization) preconditioner.
   */
  using ICPreconditioner = 
    ITL_Preconditioner<itl::pc::ic_0<MTLTypes::MTLMatrix>, 
		       MTLTypes::MTLMatrix, MTLTypes::MTLVector>;


} // namespace AMDiS
