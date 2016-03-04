#include "AMDiS.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

/// Dirichlet boundary function
class G : public AbstractFunction<double, WorldVector<double> >,
	  public TimedObject
{
public:
  
  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    return sin(M_PI * (*timePtr)) * exp(-10.0 * (x * x));
  }
};

/// RHS function
class F : public AbstractFunction<double, WorldVector<double> >,
	  public TimedObject
{
public:
  F(int degree) : AbstractFunction<double, WorldVector<double> >(degree) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    int dim = Global::getGeo(WORLD);
    double r2 = x * x;
    double ux = sin(M_PI * (*timePtr)) * exp(-10.0 * r2);
    double ut = M_PI * cos(M_PI * (*timePtr)) * exp(-10.0 * r2);
    return ut -(400.0 * r2 - 20.0 * dim) * ux;
  }
};

// ===========================================================================
// ===== instationary problem ================================================
// ===========================================================================

/// Instationary problem
class Heat : public ProblemInstat
{
public:
  Heat(ProblemStat &heatSpace)
    : ProblemInstat("heat", heatSpace)
  {
    // init theta scheme
    theta = -1.0;
    Parameters::get(name + "->theta", theta);
    TEST_EXIT(theta >= 0)("theta not set!\n");
    if (theta == 0.0) {
      WARNING("You are using the explicit Euler scheme\n");
      WARNING("Use a sufficiently small time step size!!!\n");
    }
    MSG("theta = %f\n", theta);
    theta1 = theta - 1.0;
  }

  // ===== ProblemInstatBase methods ===================================

  /// set the time in all needed functions!
  void setTime(AdaptInfo *adaptInfo) 
  {
    ProblemInstat::setTime(adaptInfo);
    rhsTime = adaptInfo->getTime() - (1 - theta) * adaptInfo->getTimestep();
  }

  void closeTimestep(AdaptInfo *adaptInfo) 
  {
    ProblemInstat::closeTimestep(adaptInfo);
    WAIT;
  }

  // ===== initial problem methods =====================================

  /// Used by \ref problemInitial to solve the system of the initial problem
  void solve(AdaptInfo *adaptInfo) 
  {
    problemStat->getSolution(0)->interpol(exactSolution);
  }

  /// Used by \ref problemInitial to do error estimation for the initial problem.
  void estimate(AdaptInfo *adaptInfo) 
  {
    double errMax, errSum;

    errSum = Error<double>::L2Err(*exactSolution,
				  *(problemStat->getSolution(0)), 0, &errMax, false);
    adaptInfo->setEstSum(errSum, 0);
    adaptInfo->setEstMax(errMax, 0);
  }

  // ===== setting methods ===============================================

  /// Sets \ref exactSolution;
  void setExactSolution(AbstractFunction<double, WorldVector<double> > *fct) 
  {
    exactSolution = fct;
  } 

  // ===== getting methods ===============================================

  /// Returns pointer to \ref theta.
  double *getThetaPtr() 
  { 
    return &theta; 
  }

  /// Returns pointer to \ref theta1.
  double *getTheta1Ptr() 
  { 
    return &theta1; 
  }

  /// Returns pointer to \ref rhsTime.
  double *getRhsTimePtr() 
  { 
    return &rhsTime; 
  }

private:
  /// Used for theta scheme.
  double theta;

  /// theta - 1
  double theta1;

  /// time for right hand side functions.
  double rhsTime;

  /// Pointer to boundary function. Needed for initial problem.
  AbstractFunction<double, WorldVector<double> > *exactSolution;
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char** argv)
{
  FUNCNAME("main");

  AMDiS::init(argc, argv);

  // ===== create and init stationary problem =====
  ProblemStat heatSpace("heat->space");
  heatSpace.initialize(INIT_ALL);


  // ===== create instationary problem =====
  Heat heat(heatSpace);
  heat.initialize(INIT_ALL);

  // create adapt info
  AdaptInfo adaptInfo("heat->adapt", heatSpace.getNumComponents());

  // create initial adapt info
  AdaptInfo adaptInfoInitial("heat->initial->adapt");

  // create instationary adapt
  AdaptInstationary adaptInstat("heat->adapt",
				heatSpace,
				adaptInfo,
				heat,
				adaptInfoInitial);


  // ===== create rhs functions =====
  int degree = heatSpace.getFeSpace()->getBasisFcts()->getDegree();
  F *rhsFct = new F(degree);
  rhsFct->setTimePtr(heat.getRhsTimePtr());


  // ===== create operators =====
  double one = 1.0;
  double zero = 0.0;

  // create laplace
  Operator A(heatSpace.getFeSpace());
  A.addTerm(new Simple_SOT);
  A.setUhOld(heat.getOldSolution(0));
  if (*(heat.getThetaPtr()) != 0.0)
    heatSpace.addMatrixOperator(A, 0, 0, heat.getThetaPtr(), &one);
  if (*(heat.getTheta1Ptr()) != 0.0)
    heatSpace.addVectorOperator(A, 0, heat.getTheta1Ptr(), &zero);

  // create zero order operator
  Operator C(heatSpace.getFeSpace());
  C.addTerm(new Simple_ZOT);
  C.setUhOld(heat.getOldSolution(0));
  heatSpace.addMatrixOperator(C, 0, 0, heat.getInvTau(), heat.getInvTau());
  heatSpace.addVectorOperator(C, 0, heat.getInvTau(), heat.getInvTau());

  // create RHS operator
  Operator F(heatSpace.getFeSpace());
  F.addTerm(new CoordsAtQP_ZOT(rhsFct));
  heatSpace.addVectorOperator(F, 0);


  // ===== create boundary functions =====
  G *boundaryFct = new G;
  boundaryFct->setTimePtr(heat.getTime());
  heat.setExactSolution(boundaryFct);
  heatSpace.addDirichletBC(1, 0, 0, boundaryFct);


  // ===== start adaption loop =====
  int errorCode = adaptInstat.adapt();

  AMDiS::finalize();

  return errorCode;
}
