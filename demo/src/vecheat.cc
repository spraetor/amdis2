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
    return ut - (400.0 * r2 - 20.0 * dim) * ux;
  };
};

// ===========================================================================
// ===== instationary problem ================================================
// ===========================================================================

/// Instationary problem
class Vecheat : public ProblemInstat
{
public:

  /// Constructor
  Vecheat(ProblemStat &vecheatSpace) 
    : ProblemInstat("vecheat", vecheatSpace)
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
    theta1 = theta - 1;
  }

  // ===== ProblemInstatBase methods ===================================

  /// set the time in all needed functions!
  void setTime(AdaptInfo *adaptInfo) 
  {
    ProblemInstat::setTime(adaptInfo);
    rhsTime = adaptInfo->getTime() - (1-theta) * adaptInfo->getTimestep();
  }

  // ===== initial problem methods =====================================

  /// Used by \ref problemInitial to solve the system of the initial problem
  void solveInitialProblem(AdaptInfo *adaptInfo) 
  {
    int size = problemStat->getNumComponents();
    rhsTime = 0.0;
    for (int i = 0; i < size; i++) {
      problemStat->getMesh(i)->dofCompress();
      problemStat->getSolution()->getDOFVector(i)->interpol(boundaryFct);
    }
  }

  /// Used by \ref problemInitial to do error estimation for the initial problem.
  void estimateInitial(AdaptInfo *adaptInfo) 
  {
    int size = problemStat->getNumComponents();
    double errMax, errSum;

    for (int i = 0; i < size; i++) {
      errSum = Error<double>::L2Err(*boundaryFct,
				    *(problemStat->getSolution()->
				      getDOFVector(i)), 
				    0, &errMax, false);
      adaptInfo->setEstSum(errSum, i);
    }
  }

  /// Used by \ref problemInitial to build before refinement.
  void buildBeforeRefineInitial(Flag) {}

  /// Used by \ref problemInitial to build before coarsening.
  void buildBeforeCoarsenInitial(Flag) {}

  /// Used by \ref problemInitial to build after coarsening.
  void buildAfterCoarsenInitial(Flag) {}

  /// Used by \ref problemInitial to mark elements.
  Flag markElementsInitial() 
  { 
    return 0; 
  }

  /// Used by \ref problemInitial
  Flag refineMeshInitial() 
  { 
    return 0; 
  }

  /// Used by \ref problemInitial
  Flag coarsenMeshInitial() 
  { 
    return 0; 
  }

  /// Used by \ref problemInitial
  void beginIterationInitial(int) {}

  /// Used by \ref problemInitial
  void endIterationInitial(int) {}

  // ===== setting methods ===============================================

  /// Sets \ref boundaryFct;
  void setBoundaryFct(AbstractFunction<double, WorldVector<double> > *fct) 
  {
    boundaryFct = fct;
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

  /// Returns \ref boundaryFct;
  AbstractFunction<double, WorldVector<double> > *getBoundaryFct() 
  {
    return boundaryFct;
  }

protected:
  /// Used for theta scheme.
  double theta;

  /// theta - 1
  double theta1;

  /// time for right hand side functions.
  double rhsTime;

  /// Pointer to boundary function. Needed for initial problem.
  AbstractFunction<double, WorldVector<double> > *boundaryFct;
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char** argv)
{
  FUNCNAME("main");

  AMDiS::init(argc, argv);

  // ===== create and init stationary problem =====
  ProblemStat vecheatSpace("vecheat->space");
  vecheatSpace.initialize(INIT_ALL);

  // ===== create instationary problem =====
  Vecheat vecheat(vecheatSpace);
  vecheat.initialize(INIT_ALL);

  // create adapt info
  AdaptInfo adaptInfo("vecheat->adapt", vecheatSpace.getNumComponents());

  // create initial adapt info
  AdaptInfo adaptInfoInitial("vecheat->initial->adapt", vecheatSpace.getNumComponents());

  // create instationary adapt
  AdaptInstationary adaptInstat("vecheat->adapt",
				vecheatSpace,
				adaptInfo,
				vecheat,
				adaptInfoInitial);

  
  // ===== create rhs functions =====
  F *rhsFct0 = new F(vecheatSpace.getFeSpace(0)->getBasisFcts()->getDegree());
  rhsFct0->setTimePtr(vecheat.getRhsTimePtr());

  F *rhsFct1 = new F(vecheatSpace.getFeSpace(1)->getBasisFcts()->getDegree());
  rhsFct1->setTimePtr(vecheat.getRhsTimePtr());


  // ===== create operators =====
  double one = 1.0;
  double zero = 0.0;

  // create laplace
  Operator A00 (vecheatSpace.getFeSpace(0), vecheatSpace.getFeSpace(0));
  A00.addTerm(new Simple_SOT);
  A00.setUhOld(vecheat.getOldSolution()->getDOFVector(0));

  Operator A11(vecheatSpace.getFeSpace(1), vecheatSpace.getFeSpace(1));
  A11.addTerm(new Simple_SOT);
  A11.setUhOld(vecheat.getOldSolution()->getDOFVector(1));

  if (*(vecheat.getThetaPtr()) != 0.0) {
    vecheatSpace.addMatrixOperator(A00, 0, 0, vecheat.getThetaPtr(), &one);
    vecheatSpace.addMatrixOperator(A11, 1, 1, vecheat.getThetaPtr(), &one);
  }

  if (*(vecheat.getTheta1Ptr()) != 0.0) {
    vecheatSpace.addVectorOperator(A00, 0, vecheat.getTheta1Ptr(), &zero);
    vecheatSpace.addVectorOperator(A11, 1, vecheat.getTheta1Ptr(), &zero);
  }

  // create zero order operator
  Operator C00(vecheatSpace.getFeSpace(0), vecheatSpace.getFeSpace(0));
  C00.addTerm(new Simple_ZOT);
  C00.setUhOld(vecheat.getOldSolution()->getDOFVector(0));
  vecheatSpace.addMatrixOperator(C00, 0, 0, 
				 vecheat.getInvTau(), vecheat.getInvTau());
  vecheatSpace.addVectorOperator(C00, 0, vecheat.getInvTau(), vecheat.getInvTau());
  
  
  Operator C11(vecheatSpace.getFeSpace(1), vecheatSpace.getFeSpace(1));
  C11.addTerm(new Simple_ZOT);
  C11.setUhOld(vecheat.getOldSolution()->getDOFVector(1));
  vecheatSpace.addMatrixOperator(C11, 1, 1, 			    
				 vecheat.getInvTau(), vecheat.getInvTau());
  vecheatSpace.addVectorOperator(C11, 1,			    
				 vecheat.getInvTau(), vecheat.getInvTau());

  // create RHS operator
  Operator F0(vecheatSpace.getFeSpace(0));
  F0.addTerm(new CoordsAtQP_ZOT(rhsFct0));
  vecheatSpace.addVectorOperator(F0, 0);

  Operator F1 = Operator(vecheatSpace.getFeSpace(1));
  F1.addTerm(new CoordsAtQP_ZOT(rhsFct1));
  vecheatSpace.addVectorOperator(F1, 1);


  // ===== create boundary functions =====
  G *boundaryFct = new G;
  boundaryFct->setTimePtr(vecheat.getTime());
  vecheat.setBoundaryFct(boundaryFct);

  vecheatSpace.addDirichletBC(1, 0, 0, boundaryFct);
  vecheatSpace.addDirichletBC(1, 1, 1, boundaryFct);


  // ===== start adaption loop =====
  int errorCode = adaptInstat.adapt();

  AMDiS::finalize();

  return errorCode;
}
