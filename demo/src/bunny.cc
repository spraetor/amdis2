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
    return -2 * x[0];
  }
};

/// boundary
class G : public AbstractFunction<double, WorldVector<double> >
{
public:

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    return 10000.0;
  }
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("bunny main");

  AMDiS::init(argc, argv);

  // ===== create projection =====
  WorldVector<double> ballCenter;
  ballCenter.set(0.0);
  new BallProject(1, VOLUME_PROJECTION, ballCenter, 1.0);

  // ===== create and init the scalar problem ===== 
  ProblemStat bunny("bunny");
  bunny.initialize(INIT_ALL);

  // === create adapt info ===
  AdaptInfo adaptInfo("bunny->adapt", bunny.getNumComponents());

  // === create adapt ===
  AdaptStationary adapt("bunny->adapt", bunny, adaptInfo);
  
  // ===== create matrix operator =====
  Operator matrixOperator(bunny.getFeSpace());
  matrixOperator.addTerm(new Simple_SOT);
  bunny.addMatrixOperator(&matrixOperator, 0, 0);

  // ===== create rhs operator =====
  Operator rhsOperator(bunny.getFeSpace());

  int degree = bunny.getFeSpace()->getBasisFcts()->getDegree();

  rhsOperator.addTerm(new CoordsAtQP_ZOT(new F(degree)));
  bunny.addVectorOperator(&rhsOperator, 0);

  // ===== start adaption loop =====
  adapt.adapt();

  bunny.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


