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


/** \file ITL_Preconditioner.h */

#ifndef AMDIS_ITL_PRECONDITIONER_H
#define AMDIS_ITL_PRECONDITIONER_H

#include "solver/LinearSolverInterface.h"
#include "MTL4Types.h"
#include "DOFMatrix.h"
#include "CreatorInterface.h"

#include "itl/masslumping.hpp"

#include <boost/numeric/itl/itl.hpp>
#include <boost/numeric/itl/pc/ilu_0.hpp>
#include <boost/numeric/itl/pc/ic_0.hpp>
#include <boost/numeric/mtl/vector/assigner.hpp>

namespace AMDiS {
  
  /**
   * \ingroup Solver
   * 
   * \brief Common base class for wrappers to use ITL preconditioners in AMDiS.
   */
  template< class MatrixType, class VectorType >
  struct ITL_PreconditionerBase : public PreconditionerInterface
  {
    virtual ~ITL_PreconditionerBase() {}
    
    virtual void init(const SolverMatrix<Matrix<DOFMatrix*> >& A, const MatrixType& fullMatrix) = 0;
    
    virtual void solve(const VectorType& x, VectorType& y) const = 0;
    
    virtual void adjoint_solve(const VectorType& x, VectorType& y) const = 0;
  };
  
  
  template< typename MatrixType, typename VectorType >
  inline itl::pc::solver<ITL_PreconditionerBase< MatrixType, VectorType >, VectorType, false>
  solve(const ITL_PreconditionerBase< MatrixType, VectorType >& P, const VectorType& vin)
  {
    return itl::pc::solver<ITL_PreconditionerBase< MatrixType, VectorType >, VectorType, false>(P, vin);
  }

  template< typename MatrixType, typename VectorType >
  inline itl::pc::solver<ITL_PreconditionerBase< MatrixType, VectorType >, VectorType, true>
  adjoint_solve(const ITL_PreconditionerBase< MatrixType, VectorType >& P, const VectorType& vin)
  {
    return itl::pc::solver<ITL_PreconditionerBase< MatrixType, VectorType >, VectorType, true>(P, vin);
  }

  
  /**
   * \ingroup Solver
   * 
   * \brief Wrapper for using ITL preconditioners in AMDiS.
   */
  template < typename Preconditioner, typename MatrixType, typename VectorType >
  class ITL_Preconditioner : public ITL_PreconditionerBase< MatrixType, VectorType >
  {
  public:
    typedef ITL_PreconditionerBase<MatrixType, VectorType>              precon_base;
    typedef ITL_Preconditioner<Preconditioner, MatrixType, VectorType>  self;
    
    /// Creator class
    struct Creator : CreatorInterfaceName<precon_base>
    {
      virtual ~Creator() {}
      precon_base* create() { return new self(); }
    };
    
    ITL_Preconditioner() 
      : precon(NULL) 
    { }
    
    ~ITL_Preconditioner()
    {
      if (precon) {
	delete precon;
	precon = NULL;
      }
    }
    
    /// Implementation of \ref ITL_PreconditionerBase::init()
    virtual void init(const SolverMatrix<Matrix<DOFMatrix*> >& A, const MatrixType& fullMatrix) override
    {
      if (precon)
	delete precon;
      precon = new Preconditioner(fullMatrix);
    }
    
    /// Implementation of \ref PreconditionerInterface::exit()
    virtual void exit() override
    {
      if (precon) {
	delete precon;
	precon = NULL;
      }
    }
    
    /// Implementation of \ref ITL_PreconditionerBase::solve()
    virtual void solve(const VectorType& vin, VectorType& vout) const override
    {
      assert(precon != NULL);
      precon->solve(vin, vout);
    }
    
    /// Implementation of \ref ITL_PreconditionerBase::adjoint_solve()
    virtual void adjoint_solve(const VectorType& vin, VectorType& vout) const override
    {
      assert(precon != NULL);
      precon->adjoint_solve(vin, vout);
    }
    
  private:
    Preconditioner* precon;
  };
  
  
  /**
   * \ingroup Solver
   * \class AMDiS::DiagonalPreconditioner
   * \brief ITL_Preconditioner implementation of diagonal (jacobi) preconditioner \implements ITL_Preconditioner
   * 
   * Diagonal preconditioner \f$ M^{-1} \f$ for the system \f$ Ax=b \f$ is defined as: \f$ M=diag(A) \f$.
   */
  typedef ITL_Preconditioner<itl::pc::diagonal<MTLTypes::MTLMatrix>, MTLTypes::MTLMatrix, MTLTypes::MTLVector > DiagonalPreconditioner;
  
  /**
   * \ingroup Solver
   * \class AMDiS::DiagonalPreconditioner
   * \brief ITL_Preconditioner implementation of diagonal (jacobi) preconditioner \implements ITL_Preconditioner
   * 
   * Diagonal preconditioner \f$ M^{-1} \f$ for the system \f$ Ax=b \f$ is defined as: \f$ M_ii=sum_j(A_ij) \f$.
   */
  typedef ITL_Preconditioner<itl::pc::masslumping<MTLTypes::MTLMatrix>, MTLTypes::MTLMatrix, MTLTypes::MTLVector > MassLumpingPreconditioner;

  
  /**
   * \ingroup Solver
   * \class AMDiS::IdentityPreconditioner
   * \brief ITL_Preconditioner implementation of identity preconditioner \implements ITL_Preconditioner
   * 
   * Identity preconditioner. Behaves like no preconditioning.
   */
  typedef ITL_Preconditioner<itl::pc::identity<MTLTypes::MTLMatrix>, MTLTypes::MTLMatrix, MTLTypes::MTLVector > IdentityPreconditioner;

  
  /**
   * \ingroup Solver
   * \class AMDiS::ILUPreconditioner
   * \brief ITL_Preconditioner implementation of ILU (Incomplete LU factorization) preconditioner. \implements ITL_Preconditioner
   * 
   * The preconditioner is used from ITL. It corresponds for instance to
   * "Iterative Methods for Sparce Linear Systems", second edition, Yousef Saad.
   *  The preconditioner is described in chapter 10.3 (algorithm 10.4).
   */
  typedef ITL_Preconditioner< itl::pc::ilu_0<MTLTypes::MTLMatrix>, MTLTypes::MTLMatrix, MTLTypes::MTLVector > ILUPreconditioner;

  
  /**
   * \ingroup Solver
   * \class AMDiS::ICPreconditioner
   * \brief ITL_Preconditioner implementation of IC (Incomplete Cholesky factorization) preconditioner. \implements ITL_Preconditioner
   * 
   * IC (Incomplete Cholesky factorization) preconditioner.
   */
  typedef ITL_Preconditioner< itl::pc::ic_0<MTLTypes::MTLMatrix>, MTLTypes::MTLMatrix, MTLTypes::MTLVector > ICPreconditioner;


} // namespace AMDiS

#endif // AMDIS_ITL_PRECONDITIONER_H
