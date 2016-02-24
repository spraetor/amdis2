/** \file LinearSolver.h */

#pragma once

#include <solver/LinearSolverInterface.h>
#include <solver/Mapper.h>
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include <parallel/ParallelSolver.h>
#endif

#include <solver/details/LinearSolver.h>
#include <traits/basic.hpp>

namespace AMDiS
{
  template <class MatrixType, class VectorType>
  struct RunnerBase : public RunnerInterface
  {
    virtual void init(SolverMatrix<Matrix<DOFMatrix*>> const& A,
		      MatrixType const& fullMatrix) = 0;

    virtual void exit() {}

    virtual int solve(MatrixType const& A, VectorType& x, VectorType const& b) = 0;

    virtual int adjoint_solve(MatrixType const& /*A*/, VectorType& /*x*/,
			      VectorType const& /*b*/)
    {
      FUNCNAME("RunnerBase::adjoint_solve()");
      ERROR_EXIT("Must be implemented in derived class!\n");
      return 0;
    }
  };


  /// Wrapper for template-argument dependent constructors
  template <class MatrixType, class Mapper_ = BlockMapper, class = void>
  struct LinearSolverBase : public LinearSolverInterface
  {
    using Mapper = Mapper_;
    using Super  = LinearSolverInterface;

    LinearSolverBase(std::string name)
      : Super(name)
    {}

    MatrixType& getMatrix()
    {
      return matrix;
    }


  protected:
    /// create a sequential BlockMapper
    virtual void initMapper(SolverMatrix<Matrix<DOFMatrix*>> const& A)
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
  template <class MatrixType>
  struct LinearSolverBase<MatrixType, ParallelMapper,
			  Requires_t<mtl::traits::is_distributed<MatrixType>>>
    : public ParallelSolver
  {
    using Mapper = ParallelMapper;
    using Super  = ParallelSolver;

    LinearSolverBase(std::string name)
      : Super(name, false)
    {}

    MatrixType& getMatrix()
    {
      return matrix;
    }

  protected:

    /// create a parallel mapper based on local-to-global mapping
    virtual void initMapper(SolverMatrix<Matrix<DOFMatrix*>> const& A)
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
  template <class MatrixType, class VectorType, class Runner,
	    class Mapper_ = BlockMapper>
  class LinearSolver : public LinearSolverBase<MatrixType, Mapper_>
  {
  protected:
    using Self   = LinearSolver;
    using Super  = LinearSolverBase<MatrixType, Mapper_>;
    using Mapper = typename Super::Mapper;

  public:
    /// Creator class used in the LinearSolverInterfaceMap.
    struct Creator : public LinearSolverCreator
    {
      /// Returns a new LinearSolver object.
      virtual LinearSolverInterface* create() override
      {
        return new Self(this->name);
      }
    };

    /// Constructor
    LinearSolver(std::string name)
      : Super(name),
        runner(this)
    { }


    /// Implementation of \ref LinearSolverInterface::getRunner()
    virtual RunnerInterface* getRunner() override
    {
      return &runner;
    }


    /// Implementation of \ref LinearSolverInterface::getLeftPrecon()
    virtual PreconditionerInterface* getLeftPrecon() override
    {
      return runner.getLeftPrecon();
    }


    /// Implementation of \ref LinearSolverInterface::getRightPrecon()
    virtual PreconditionerInterface* getRightPrecon() override
    {
      return runner.getRightPrecon();
    }

  protected:
    /// Implementation of \ref LinearSolverInterface::solveSystemImpl()
    virtual int solveSystemImpl(SolverMatrix<Matrix<DOFMatrix*>> const& A,
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
      if (createMatrixData)
      {
        initMatrix(Super::matrix, A, *Super::mapper);
        runner.init(A, Super::matrix);
      }

      VectorType mtl_x;
      initVector(mtl_x, x, *Super::mapper);

      VectorType mtl_b;
      initVector(mtl_b, b, *Super::mapper);

      INFO(Self::getInfo(), 8)("fill MTL4 matrix needed %.5f seconds\n", t.elapsed());

      int error = runner.solve(Super::matrix ,mtl_x, mtl_b);

      VecMap<SystemVector, Mapper> xVecMap{x, *Super::mapper};
      mtl_x >> xVecMap;

      if (!storeMatrixData)
        runner.exit();

      this->exitMapper();
      return error;
    }

    // functions to initialize mtl-matrix and mtl-vector
    // from AMDiS matrix / vectors using mappers

    /// initialize a MTL matrix and assign values from an AMDiS matrix
    template <class Matrix1, class Matrix2, class M>
    void initMatrix(Matrix1& target, Matrix2 const& source, MapperBase<M>& mapper)
    {
      dispatch::initMatrix(target, mapper.self());
      dispatch::fillMatrix(target, source, mapper.self());
    }

    /// initialize a MTL vector and assign values from an AMDiS vector
    template <class Vector1, class Vector2, class M>
    void initVector(Vector1& target, Vector2 const& source, MapperBase<M>& mapper)
    {
      dispatch::initVector(target, Super::matrix);
      dispatch::fillVector(target, source, mapper.self());
    }

    Runner runner; // redirect the implementation to a runner

  };

} // end namespace AMDiS
