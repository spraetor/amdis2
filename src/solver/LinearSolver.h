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


/** \file LinearSolver.h */

#ifndef AMDIS_LINEAR_SOLVER_BASE_H
#define AMDIS_LINEAR_SOLVER_BASE_H

#include "solver/LinearSolverInterface.h"
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/ParallelSolver.h"
#endif

#include "solver/details/LinearSolver.h"

namespace AMDiS {
  
  template< typename MatrixType, typename VectorType >
  struct RunnerBase : public RunnerInterface
  {
    virtual void init(const SolverMatrix<Matrix<DOFMatrix*> >& A, const MatrixType& fullMatrix) = 0;      
    
    virtual int solve(const MatrixType& A, VectorType& x, const VectorType& b) = 0;  
    
    virtual int adjoint_solve(const MatrixType& A, VectorType& x, const VectorType& b) 
    {
      FUNCNAME("RunnerBase::adjoint_solve()");
      ERROR_EXIT("Must be implemented in derived class!\n");
      return 0;
    }
  };
  
  
  /// Wrapper for template-argument dependent constructors
  template < typename MatrixType, typename Mapper_ = BlockMapper, typename Enable = void >
  struct LinearSolverBase : public LinearSolverInterface
  {
    typedef Mapper_ Mapper;
  
    LinearSolverBase(std::string name) 
      : LinearSolverInterface(name) {}    
    
    MatrixType& getMatrix() 
    {
      return matrix;
    }
    

  protected:
    /// create a sequential BlockMapper
    virtual void initMapper(const SolverMatrix<Matrix<DOFMatrix*> >& A) 
    {
      mapper = new Mapper(A);
    }
    
    virtual void exitMapper()
    {
      delete mapper;
    }

    MatrixType matrix;
    Mapper* mapper;
  };
  
#ifdef HAVE_PARALLEL_MTL4
  template< typename MatrixType >
  struct LinearSolverBase<MatrixType, ParallelMapper, typename enable_if< mtl::traits::is_distributed<MatrixType> > > 
    : public ParallelSolver
  {
    typedef ParallelMapper Mapper;
    
    LinearSolverBase(std::string name) 
      : ParallelSolver(name, false) {}    
    
    MatrixType& getMatrix() 
    {
      return matrix;
    }
    
  protected:
    
    /// create a parallel mapper based on local-to-global mapping
    virtual void initMapper(const SolverMatrix<Matrix<DOFMatrix*> >& A)
    {
      mapper = new ParallelMapper(*ParallelSolver::getDofMapping());
    }
    
    virtual void exitMapper()
    {
      delete mapper;
    }
    MatrixType matrix;
    Mapper* mapper;
  };
#endif

  /** \ingroup Solver
   * 
   * \brief
   * Wrapper class for various MTL4 solvers. These solvers
   * are parametrized by MatrixType and VectorType. The different
   * algorithms, like krylov subspace methods or other external
   * solvers where MTL4 provides an interface, can be assigned
   * by different Runner objects.
   **/
  template< typename MatrixType, typename VectorType, typename Runner, typename Mapper_ = BlockMapper >
  class LinearSolver : public LinearSolverBase<MatrixType, Mapper_>
  {    
  protected:
    typedef LinearSolverBase<MatrixType, Mapper_>                  super;
    typedef LinearSolver<MatrixType, VectorType, Runner, Mapper_>  self;
    typedef typename super::Mapper                                 Mapper;
    
  public:
    /// Creator class used in the LinearSolverInterfaceMap.
    class Creator : public LinearSolverCreator
    {
    public:
      virtual ~Creator() {}
      
      /// Returns a new LinearSolver object.
      LinearSolverInterface* create() 
      { 
	return new self(this->name); 
      }
    };
    
    /// Constructor
    LinearSolver(std::string name)
      : super(name),
	runner(this)
    {}
    
    
    /// Implementation of \ref LinearSolverInterface::getRunner()
    RunnerInterface* getRunner()
    {
      return &runner;
    }
    
    
    /// Implementation of \ref LinearSolverInterface::getLeftPrecon()
    PreconditionerInterface* getLeftPrecon()
    {
      return runner.getLeftPrecon();
    }
    
    
    /// Implementation of \ref LinearSolverInterface::getRightPrecon()
    PreconditionerInterface* getRightPrecon()
    {
      return runner.getRightPrecon();
    }
    
  protected:  

    /// Implementation of \ref LinearSolverInterface::solveLinearSystem()
    virtual int solveLinearSystem(const SolverMatrix<Matrix<DOFMatrix*> >& A,
				  SystemVector& x,
				  SystemVector& b,
				  bool createMatrixData,
				  bool storeMatrixData) override
    {    
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
      MPI::COMM_WORLD.Barrier();
#endif
      this->initMapper(A);
      
      Timer t;
      if (createMatrixData) {
	initMatrix(super::matrix, A, *super::mapper);  
	runner.init(A, super::matrix);
      }

      VectorType mtl_x;
      initVector(mtl_x, x, *super::mapper);

      VectorType mtl_b;
      initVector(mtl_b, b, *super::mapper);
   
      INFO(self::getInfo(), 8)("fill MTL4 matrix needed %.5f seconds\n", t.elapsed());

      int error = runner.solve(super::matrix ,mtl_x, mtl_b);
      
      VecMap<SystemVector, Mapper> xVecMap(x, *super::mapper);
      mtl_x >> xVecMap;

      if (!storeMatrixData)
	runner.exit();

      this->exitMapper();
      return error;
    }
    
    // functions to initialize mtl-matrix and mtl-vector
    // from AMDiS matrix / vectors using mappers
    
    /// initialize a MTL matrix and assign values from an AMDiS matrix
    template< typename Matrix1, typename Matrix2, typename M >
    void initMatrix(Matrix1& target, const Matrix2& source, MapperBase<M>& mapper) 
    {
      dispatch::initMatrix(target, mapper.self());
      dispatch::fillMatrix(target, source, mapper.self());
    }
    
    /// initialize a MTL vector and assign values from an AMDiS vector
    template< typename Vector1, typename Vector2, typename M >
    void initVector(Vector1& target, const Vector2& source, MapperBase<M>& mapper) 
    {
      dispatch::initVector(target, super::matrix);
      dispatch::fillVector(target, source, mapper.self());
    }

    Runner runner; // redirect the implementation to a runner
    
  };
}

#endif // AMDIS_LINEAR_SOLVER_BASE_H
