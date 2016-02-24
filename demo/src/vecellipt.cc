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
  F(int degree = 1) : AbstractFunction<double, WorldVector<double> >(degree) {}

  /// Implementation of AbstractFunction::operator().
  double operator()(const WorldVector<double>& x) const 
  {
    int dim = Global::getGeo(WORLD);
    double r2 = (x * x);
    double ux = exp(-10.0 * r2);
    return -(400.0 * r2 - 20.0 * dim) * ux;
  };
};

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

int main(int argc, char* argv[])
{
  FUNCNAME("main");

  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat vecellipt("vecellipt");
  vecellipt.initialize(INIT_ALL);


  // === create adapt info ===
  AdaptInfo adaptInfo("vecellipt->adapt", vecellipt.getNumComponents());


  // === create adapt ===
  AdaptStationary adapt("vecellipt->adapt", vecellipt, adaptInfo);

  
  // ===== create matrix operators =====
  Operator matrixOperator00(vecellipt.getFeSpace(0), vecellipt.getFeSpace(0));
  matrixOperator00.addTerm(new Simple_ZOT);
  vecellipt.addMatrixOperator(matrixOperator00, 0, 0);

  Operator matrixOperator10(vecellipt.getFeSpace(1), vecellipt.getFeSpace(0));
  matrixOperator10.addTerm(new Simple_ZOT(-1.0));
  vecellipt.addMatrixOperator(matrixOperator10, 1, 0);

  Operator matrixOperator11(vecellipt.getFeSpace(0), vecellipt.getFeSpace(0));
  matrixOperator11.addTerm(new Simple_SOT);
  vecellipt.addMatrixOperator(matrixOperator11, 1, 1);

  int degree = vecellipt.getFeSpace(0)->getBasisFcts()->getDegree();
  
  // ===== create rhs operator =====
  Operator rhsOperator0(vecellipt.getFeSpace(0));
  rhsOperator0.addTerm(new CoordsAtQP_ZOT(new F(degree)));
  vecellipt.addVectorOperator(rhsOperator0, 0);

  // ===== add boundary conditions =====
  vecellipt.addDirichletBC(1, 1, 1, new G);


  // ===== start adaption loop =====
  adapt.adapt();

  vecellipt.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


