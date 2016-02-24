#include "AMDiS.h"

using namespace std;
using namespace AMDiS;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

/// RHS function
class F : public AbstractFunction<double, WorldVector<double> >
{
public:
  F(int degree) : AbstractFunction<double, WorldVector<double> >(degree) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    return -2.0 * x[0];
  }
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("sphere main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  WorldVector<double> ballCenter;
  ballCenter.set(0.0);
  new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);

  // ===== create and init the scalar problem ===== 
  ProblemStat sphere("sphere");
  sphere.initialize(INIT_ALL);

  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("sphere->adapt", sphere.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("sphere->adapt",
					       &sphere,
					       adaptInfo);
  
  // ===== create matrix operator =====
  Operator matrixOperator(sphere.getFeSpace());
  matrixOperator.addTerm(new Simple_SOT);
  sphere.addMatrixOperator(&matrixOperator, 0, 0);

  // ===== create rhs operator =====
  Operator rhsOperator(sphere.getFeSpace());

  int degree = sphere.getFeSpace()->getBasisFcts()->getDegree();

  rhsOperator.addTerm(new CoordsAtQP_ZOT(new F(degree)));
  sphere.addVectorOperator(&rhsOperator, 0);

  // ===== start adaption loop =====
  adapt->adapt();

  sphere.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


