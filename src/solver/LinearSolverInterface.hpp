/** \file LinearSolverInterface.h */

/**
 * \defgroup Solver Solver module
 * @{ <img src="solver.png"> @}
 *
 * \brief
 * Contains all classes needed for solving linear and non linear equation
 * systems.
 */

#pragma once

// AMDiS includes
// #include "Global.hpp"
#include "AMDiS_fwd.hpp"
#include "CreatorInterface.hpp"
#include "Log.hpp"

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
// #include "parallel/ParallelMapper.hpp"
#endif

namespace AMDiS
{

  /// Non-templated base-class for Preconditioner types
  struct PreconditionerInterface
  {
    virtual ~PreconditionerInterface() {}
  };

  /// Non-templates base-class for Runner / Worker types
  struct RunnerInterface
  {
    // virtual destructor
    virtual ~RunnerInterface() {}

    virtual PreconditionerInterface* getLeftPrecon()
    {
      return NULL;
    }

    virtual PreconditionerInterface* getRightPrecon()
    {
      return NULL;
    }
  };

  /**
   * \ingroup Solver
   *
   *\brief
   * Solver for linear equation systems.
   */
  class LinearSolverInterface
  {
  public:
    /// The constructor reads needed parameters and sets solvers \ref name.
    /**
     * Reads parameters for a solver with name 'NAME':
     *   NAME->tolerance            \ref tolerance
     *   NAME->relative tolerance   \ref relative
     *   NAME->max iteration        \ref max_iter
     *   NAME->print cycle          \ref print_cycle
     *   NAME->info                 \ref info
     *   NAME->break if tolerance not reached \ref breakTolNotReached
    **/
    explicit LinearSolverInterface(std::string name);

    /// Destructor.
    virtual ~LinearSolverInterface() {}


    /// Public method to call in order to solve a linear system Ax = b.
    /**
     * The method redirects to the specific linear solver and prints statistics
     * and error estimations at the end.
     * 
     * The parameters correspond to
     *  \p A     A block-matrix that represents the system-matrix.
     *  \p x     A block-vector for the unknown components.
     *  \p b     A block-vector for the right-hand side of the linear system.
     *  \p createMatrixData  If true, the matrix will be initialized and the 
     * 	                     corresponding runner of the system receives the 
     *                       matrix in the init() method.
     *  \p storeMatrixData   If true, the exit() method of the runner will be 
     *                       called.
     **/
    int solveSystem(SolverMatrix<Matrix<DOFMatrix*>> const& A,
                    SystemVector& x,
                    SystemVector& b,
                    bool createMatrixData,
                    bool storeMatrixData);

    /** \name getting methods
     * \{
     */

    /// Returns solvers \ref name.
    std::string getName() const
    {
      return name;
    }

    /// Returns \ref tolerance
    double getTolerance() const
    {
      return tolerance;
    }

    /// Returns \ref max_iter
    int getMaxIterations() const
    {
      return max_iter;
    }

    /// Returns number of iterations in last run of an iterative solver
    int getIterations() const
    {
      return iterations;
    }

    /// Returns error code in last run of an iterative solver
    int getErrorCode() const
    {
      return error;
    }

    /// Returns info
    int getInfo() const
    {
      return info;
    }

    /// Returns \ref print_cycle
    int getPrint_cycle() const
    {
      return print_cycle;
    }

    /// Returns \ref residual
    double getResidual() const
    {
      return residual;
    }

    /// Returns \ref relative
    double getRelative() const
    {
      return relative;
    }


    virtual PreconditionerInterface* getLeftPrecon()
    {
      FUNCNAME("LinearSolverInterface::getLeftPrecon()");
      ERROR("no left preconditioner provided!\n");
      return NULL;
    }

    virtual PreconditionerInterface* getRightPrecon()
    {
      FUNCNAME("LinearSolverInterface::getRightPrecon()");
      ERROR("no right preconditioner provided!\n");
      return NULL;
    }

    virtual RunnerInterface* getRunner()
    {
      FUNCNAME("LinearSolverInterface::getRunner()");
      ERROR("no runner provided!\n");
      return NULL;
    }

    /** \} */

    /** \name setting methods
     * \{
     */

    /// Sets \ref tolerance
    void setTolerance(double tol)
    {
      tolerance = tol;
    }

    void setResidual(double r)
    {
      residual = r;
    }

    /// Sets \ref relative
    void setRelative(double rel)
    {
      relative = rel;
    }

    void setRelativeResidual(double r)
    {
      rel_residual = r;
    }

    /// Sets \ref max_iter
    void setMaxIterations(int i)
    {
      max_iter = i;
    }

    void setIterations(int i)
    {
      iterations=i;
    }

    /// set the \ref error
    void setErrorCode(int code)
    {
      error=code;
    }

    /// Sets \ref info
    void setInfo(int i)
    {
      info = i;
    }

    /** \} */

  protected:
    /// main methods that all solvers must implement
    virtual int solveSystemImpl(SolverMatrix<Matrix<DOFMatrix*>> const& A,
				SystemVector& x,
				SystemVector& b,
				bool createMatrixData,
				bool storeMatrixData) = 0;

  protected:
    /// solvers name.
    std::string name;

    /// Solver tolerance |r|. Set in LinearSolverInterface's constructor.
    double tolerance;

    /// Relative solver tolerance |r|/|r0|. Set in LinearSolverInterface's constructor.
    double relative;

    /// maximal number of iterations. Set in LinearSolverInterface's constructor.
    int max_iter;

    /// info level during solving the system. Set in LinearSolverInterface's constructor.
    int info;

    /// current absolute residual norm.
    double residual;

    /// current relative residual norm.
    double rel_residual;

    /// Print cycle, after how many iterations the residuum norm is logged.
    int print_cycle;

    /// How many iterations were performed in last solver (not set by UmfPack)
    int iterations;

    /// Error code  in last solver (not set by UmfPack)
    int error;

    /// break if residual norm > prescribed tolerance
    bool breakTolNotReached;
  };


  using LinearSolverCreator = CreatorInterfaceName<LinearSolverInterface>;

} // end namespace AMDiS
