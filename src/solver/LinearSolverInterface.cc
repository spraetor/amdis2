#include "LinearSolverInterface.h"


// AMDiS includes
#include <Initfile.h>

namespace AMDiS
{

  LinearSolverInterface::LinearSolverInterface(std::string name_)
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


  int LinearSolverInterface::solveSystem(SolverMatrix<Matrix<DOFMatrix*>> const& A,
					 SystemVector& x,
					 SystemVector& b,
					 bool createMatrixData,
					 bool storeMatrixData)
  {
    FUNCNAME("LinearSolverInterface::solveSystem()");
    MSG("LinearSolverInterface::solveSystem()\n");

    residual = -1.0;
    rel_residual = -1.0;
    int error_code = solveSystemImpl(A, x, b, createMatrixData, storeMatrixData);

    // calculate and print resiual
    if (info > 0)
    {
      if (residual >= 0.0 && rel_residual >= 0.0)
      {
	MSG("Residual norm: ||b-Ax|| = %e, ||b-Ax||/||b|| = %e\n", 
	    residual, rel_residual);
      }
      else if (residual >= 0.0)
      {
	MSG("Residual norm: ||b-Ax|| = %e\n", residual);
      }

#if DEBUG != 0
      if (getIterations() > 0)
      {
	MSG("Nr. of iterations needed = %d\n", getIterations());
      }

      if (error_code != 0)
      {
	MSG("ERROR-Code = %d\n", error_code);
      }

      if (!isNumber(residual) || !isNumber(rel_residual))
      {
	MSG("Residual or relative residual is NaN/Inf!\n");
      }
#endif

      // test for absolute tolerance
      TEST_EXIT((isNumber(residual) && (tolerance < 1.e-30 || residual <= tolerance))
		|| !breakTolNotReached)
      ("Tolerance tol = %e could not be reached!\n Set tolerance by '->solver->tolerance:' \n", 
	tolerance);

      // test for relative tolerance
      TEST_EXIT((isNumber(rel_residual) && (relative < 1.e-30 || rel_residual <= relative))
		|| (residual < 1.e-30) || !breakTolNotReached)
      ("Relative tolerance rtol = %e could not be reached!\n Set tolerance by '->solver->relative tolerance:' \n", 
	relative);
    }
    return error_code;
  }

} // end namespace AMDiS
