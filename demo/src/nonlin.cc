#include "AMDiS.h"

using namespace AMDiS;
using namespace std;

// ===========================================================================
// ===== function definitions ================================================
// ===========================================================================

/// Dirichlet boundary function
class G : public AbstractFunction<double, WorldVector<double> >
{
public:

  G()
    : AbstractFunction<double, WorldVector<double> >()
  {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    return 0.0;
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
    int dow = Global::getGeo(WORLD);
    double r2 = (x * x);
    double ux = exp(-10.0 * r2);
    return -(400.0 * r2 - 20.0 * dow) * ux + pow(ux, 4.0);
  }
};


class NonlinFctLeft : public AbstractFunction<double, double> 
{
public:
  NonlinFctLeft()
    : AbstractFunction<double, double>()
  {}

  double operator()(const double& u) const
  {
    return 4.0 * u * u * u;
  }
};


class NonlinFctRight : public AbstractFunction<double, double> 
{
public:
  NonlinFctRight()
    : AbstractFunction<double, double>()
  {}

  double operator()(const double& u) const
  {
    return -(u * u * u * u);
  }
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("main");

  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemNonLin nonlin("nonlin");
  nonlin.initialize(INIT_ALL);

  // === create adapt info ===
  AdaptInfo adaptInfo("nonlin->adapt");


  // === create adapt ===
  AdaptStationary adapt("nonlin->adapt", nonlin, adaptInfo);

  
  // === create matrix operators ===
  Operator mat01(nonlin.getFeSpace());
  mat01.addTerm(new Simple_SOT());
  nonlin.addMatrixOperator(mat01, 0, 0);

  Operator mat02(nonlin.getFeSpace());
  mat02.addTerm(new VecAtQP_ZOT(nonlin.getSolution(0), new NonlinFctLeft()));
  nonlin.addMatrixOperator(mat02, 0, 0);

  // === create rhs operators ===
  Operator vec01(nonlin.getFeSpace());
  vec01.setUhOld(nonlin.getSolution(0));
  vec01.addTerm(new Simple_SOT(-1.0));
  nonlin.addVectorOperator(vec01, 0);
  
  Operator vec02(nonlin.getFeSpace());
  vec02.addTerm(new VecAtQP_ZOT(nonlin.getSolution(0), new NonlinFctRight()));
  nonlin.addVectorOperator(vec02, 0);
  
  Operator vec03(nonlin.getFeSpace());
  vec03.addTerm(new CoordsAtQP_ZOT(new F(nonlin.getFeSpace()->getBasisFcts()->getDegree())));
  nonlin.addVectorOperator(vec03, 0);

  // ===== add boundary conditions =====
  nonlin.addDirichletBC(1, 0, 0, new G);

  // ===== start adaption loop =====
  adapt.adapt();

  nonlin.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}
