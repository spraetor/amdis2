#include "AMDiS.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

/// Dirichlet boundary function
class G : public AbstractFunction<double, WorldVector<double> >
{
public:
  
  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    return exp(-10.0 * (x * x));
  }
};

/// RHS function
class F : public AbstractFunction<double, WorldVector<double> >
{
public:
  F(int degree) : AbstractFunction<double, WorldVector<double> >(degree) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    int dim = x.getSize();
    double r2 = x * x;
    double ux = exp(-10.0*r2);
    return -(400.0 * r2 - 20.0 * dim) * ux;
  }
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  WorldVector<double> ballCenter;
  ballCenter.set(0.0);
  new BallProject(1, 
		  BOUNDARY_PROJECTION, 
		  ballCenter, 
		  1.0);

  // ===== create and init the scalar problem ===== 
  ProblemStat ball("ball");
  ball.initialize(INIT_ALL);

  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("ball->adapt", ball.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("ball->adapt", &ball, adaptInfo);
  
  // ===== create matrix operator =====
  Operator matrixOperator(ball.getFeSpace());
  matrixOperator.addTerm(new Simple_SOT);
  ball.addMatrixOperator(&matrixOperator, 0, 0);

  // ===== create rhs operator =====
  int degree = ball.getFeSpace()->getBasisFcts()->getDegree();
  Operator rhsOperator(ball.getFeSpace());
  rhsOperator.addTerm(new CoordsAtQP_ZOT(new F(degree)));
  ball.addVectorOperator(&rhsOperator, 0);

  // ===== add boundary conditions =====
  ball.addDirichletBC(1, 0, 0, new G);

  // ===== start adaption loop =====
  adapt->adapt();

  ball.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


