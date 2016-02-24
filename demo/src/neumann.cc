#include "AMDiS.h"

using namespace AMDiS;
using namespace std;

class N : public AbstractFunction<double, WorldVector<double> >
{
public:
  double operator()(const WorldVector<double>& x) const 
  {
    return 1.0;
  }
};

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
    int dow = x.getSize();
    double r2 = (x*x);
    double ux = exp(-10.0*r2);
    return -(400.0*r2 - 20.0*dow)*ux;
  }
};

// // ===== main program // 
int main(int argc, char* argv[])
{
  FUNCNAME("main");

  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat neumann("neumann");
  neumann.initialize(INIT_ALL);

  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("neumann->adapt", neumann.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("neumann->adapt",
					       &neumann,
					       adaptInfo);
  
  // ===== create matrix operator =====
  Operator matrixOperator(neumann.getFeSpace());
  matrixOperator.addTerm(new Simple_SOT);
  neumann.addMatrixOperator(&matrixOperator, 0, 0);

  // ===== create rhs operator =====
  int degree = neumann.getFeSpace()->getBasisFcts()->getDegree();
  Operator rhsOperator(neumann.getFeSpace());
  rhsOperator.addTerm(new CoordsAtQP_ZOT(new F(degree)));
  neumann.addVectorOperator(&rhsOperator, 0);

  // ===== add boundary conditions =====
  neumann.addNeumannBC(1, 0, 0, new N);
  neumann.addDirichletBC(2, 0, 0, new G);

  // ===== start adaption loop =====
  adapt->adapt();

  neumann.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


