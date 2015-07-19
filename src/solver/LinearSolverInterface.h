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

#include "Global.h"
#include "AMDiS_fwd.h"
#include "Initfile.h"
#include "solver/SolverMatrix.h"
#include "DOFVector.h"
#include "SystemVector.h"
#include "DOFMatrix.h"
#include "Mapper.h"
#include "Timer.h"

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/ParallelMapper.h"
#endif

namespace AMDiS {
  
  /// Non-templated base-class for Preconditioner types
  struct PreconditionerInterface
  {
    virtual ~PreconditionerInterface() {}
    virtual void exit() {}
  };

  /// Non-templates base-class for Runner / Worker types
  struct RunnerInterface
  {
    // virtual destructor
    virtual ~RunnerInterface() {}
    
    virtual void exit() {}
    
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
    LinearSolverInterface(std::string name_) 
      : name(name_),
      	tolerance(DBL_TOL),
      	relative(0),
      	max_iter(1000),
      	info(0),
      	residual(-1.0),
      	rel_residual(-1.0),
      	print_cycle(100),
      	iterations(-1),
      	error(-1),
      	breakTolNotReached(true)
    {      
      Parameters::get(name + "->tolerance", tolerance);
      Parameters::get(name + "->relative tolerance", relative);
      Parameters::get(name + "->max iteration", max_iter);
      Parameters::get(name + "->print cycle", print_cycle);
      Parameters::get(name + "->info", info);
      Parameters::get(name + "->break if tolerance not reached", breakTolNotReached);
    }

    /// destructor
    virtual ~LinearSolverInterface() { }


    int solveSystem(const SolverMatrix<Matrix<DOFMatrix*> >& A,
            		    SystemVector& x,
            		    SystemVector& b,
            		    bool createMatrixData,
            		    bool storeMatrixData)
    { FUNCNAME("LinearSolverInterface::solveSystem()");
      MSG("LinearSolverInterface::solveSystem()\n");
    
      residual = -1.0;
      rel_residual = -1.0;
      int error_code = solveLinearSystem(A, x, b, createMatrixData, storeMatrixData);
      
      // calculate and print resiual
      if (info > 0) {
      	if (residual >= 0.0 && rel_residual >= 0.0) {
      	  MSG("Residual norm: ||b-Ax|| = %e, ||b-Ax||/||b|| = %e\n", residual, rel_residual);
      	} else if (residual >= 0.0) {
      	  MSG("Residual norm: ||b-Ax|| = %e\n", residual);
      	}
      	
#if DEBUG != 0
      	if (getIterations() > 0) {
      	  MSG("Nr. of iterations needed = %d\n", getIterations());
      	}
      	
      	if (error_code != 0) {
      	  MSG("ERROR-Code = %d\n", error_code);
      	}
      	
      	if (!isNumber(residual) || !isNumber(rel_residual)) {
      	  MSG("Residual or relative residual is NaN/Inf!\n");
      	}
#endif
      
      	// test for absolute tolerance
      	TEST_EXIT((isNumber(residual) && (residual < 0.0  || tolerance < 1.e-30 || residual <= tolerance))
      		  || !breakTolNotReached)
      	  ("Tolerance tol = %e could not be reached!\n Set tolerance by '->solver->tolerance:' \n", tolerance);
      	  
      	// test for relative tolerance
      	TEST_EXIT((isNumber(rel_residual) && (rel_residual < 0.0  || relative < 1.e-30 || rel_residual <= relative))
      		  || (residual < 1.e-30) || !breakTolNotReached)
      	  ("Relative tolerance rtol = %e could not be reached!\n Set tolerance by '->solver->relative tolerance:' \n", relative);
      }
      return error_code;
    }

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
    int getIterations() 
    {
      return iterations;
    }

    /// Returns error code in last run of an iterative solver
    int getErrorCode() 
    {
      return error;
    }

    /// Returns info
    int getInfo()
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
    virtual int solveLinearSystem(const SolverMatrix<Matrix<DOFMatrix*> >& A,
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
  
  
  typedef CreatorInterfaceName<LinearSolverInterface> LinearSolverCreator;
  
} // end namespace AMDiS
