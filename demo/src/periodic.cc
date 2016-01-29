#include "AMDiS.h"

using namespace std;
using namespace AMDiS;

/// Dirichlet boundary function
class G : public AbstractFunction<double, WorldVector<double> >
{
public:
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
    int dim = x.getSize();
    double r2 = (x*x);
    double ux = exp(-10.0*r2);
    return -(400.0*r2 - 20.0*dim)*ux;
  }
};

int main(int argc, char* argv[])
{
  FUNCNAME("main");

  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat periodic("periodic");
  periodic.initialize(INIT_ALL);

  // === create adapt info ===
  AdaptInfo *adaptInfo = new AdaptInfo("periodic->adapt", periodic.getNumComponents());

  // === create adapt ===
  AdaptStationary *adapt = new AdaptStationary("periodic->adapt",
					       &periodic,
					       adaptInfo);
  
  // ===== create matrix operator =====
  Operator matrixOperator(periodic.getFeSpace());
  matrixOperator.addTerm(new Simple_SOT);
  periodic.addMatrixOperator(&matrixOperator, 0, 0);

  // ===== create rhs operator =====
  int degree = periodic.getFeSpace()->getBasisFcts()->getDegree();
  Operator rhsOperator(periodic.getFeSpace());
  rhsOperator.addTerm(new CoordsAtQP_ZOT(new F(degree)));
  periodic.addVectorOperator(&rhsOperator, 0);

  // ===== add boundary conditions =====
  periodic.addDirichletBC(1, 0, 0, new G);

  // ===== add periodic conditions =====
  periodic.addPeriodicBC(-1, 0, 0);

  // ===== start adaption loop =====
  adapt->adapt();

  periodic.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


