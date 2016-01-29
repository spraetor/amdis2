#include "AMDiS.h"

using namespace AMDiS;

// ===========================================================================
// ===== main program ========================================================
// ===========================================================================

struct Dot
{
  template <class V1, class V2>
  auto operator()(V1 const& v1, V2 const& v2) const RETURNS( v1*v2 )
};

int main(int argc, char* argv[])
{
  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat ellipt("ellipt");
  ellipt.initialize(INIT_ALL);

  // === create adapt info ===
  AdaptInfo adaptInfo("ellipt->adapt");

  // === create adapt ===
  AdaptStationary adapt("ellipt->adapt", ellipt, adaptInfo);
  
  // ===== create matrix operator =====
  Operator matrixOperator(ellipt.getFeSpace());
  matrixOperator.addSecondOrderTerm(1.0);
  ellipt.addMatrixOperator(matrixOperator, 0, 0);

  // ===== create rhs operator =====
  Operator rhsOperator(ellipt.getFeSpace());
//   rhsOperator.addZeroOrderTerm(func(Dot(), X(), X()));
  rhsOperator.addZeroOrderTerm(X() * X());
  ellipt.addVectorOperator(rhsOperator, 0);
  
  // ===== add boundary conditions =====
  ellipt.addDirichletBC(1, 0, 0, 0.0);


  // ===== start adaption loop =====
  adapt.adapt();
  
  ellipt.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


