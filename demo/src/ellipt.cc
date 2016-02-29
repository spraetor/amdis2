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

/// write Systemmatrix to file in matlab-format
void writeMatrix(ProblemStat& prob, std::string filename)
{
  mtl::io::matrix_market_ostream out(filename);
  SolverMatrix<Matrix<DOFMatrix*> > solverMatrix;
  solverMatrix.setMatrix(*prob.getSystemMatrix());
  out << solverMatrix.getMatrix();
  out.close();
}

void printSolution(DOFVector<double>& dofvec)
{
  std::cout << "[";
  DOFIterator<double> it(&dofvec, USED_DOFS);
  for (it.reset(); !it.end(); ++it)
    std::cout << *it << " ";
  std::cout << "]\n";
}

int main(int argc, char* argv[])
{
  AMDiS::init(argc, argv);

  // ===== create and init the scalar problem ===== 
  ProblemStat ellipt("ellipt");
  ellipt.initialize(INIT_ALL);
  
  MSG(" NMACROS = %d\n", ellipt.getMesh()->getNumberOfMacros());
  MSG(" NDOFS = %d\n", ellipt.getMesh()->getNumberOfDofs(0));
  MSG(" NALLDOFS = %d\n", ellipt.getMesh()->getNumberOfAllDofs());
  MSG(" NELEMENTS = %d\n", ellipt.getMesh()->getNumberOfElements());
  MSG(" NLEAVES = %d\n", ellipt.getMesh()->getNumberOfLeaves());
  MSG(" NFACES = %d\n", ellipt.getMesh()->getNumberOfFaces());
  MSG(" NEDGES = %d\n", ellipt.getMesh()->getNumberOfEdges());
  MSG(" NVERTICES = %d\n", ellipt.getMesh()->getNumberOfVertices());
  MSG(" NNODES = %d\n", ellipt.getMesh()->getNumberOfNodes());
  
  WorldVector<double> diam = ellipt.getMesh()->getDiameter();
  MSG(" DIAMETER = (%f, %f)\n", diam[0], diam[1]);
  WorldVector<double> diam2; diam2 = ellipt.getMesh()->getDiameter();
  MSG(" DIAMETER2 = (%f, %f)\n", diam2[0], diam2[1]);
  
  MSG(" USED_SIZE = %d\n", ellipt.getSolution(0)->getUsedSize());

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
  ellipt.addDirichletBC(4, 0, 0, 0.0);


  // ===== start adaption loop =====
  adapt.adapt();
  writeMatrix(ellipt, "matrix.mtx");
  printSolution(*ellipt.getSolution(0));
  printSolution(*ellipt.getRhsVector(0));

  
  ellipt.writeFiles(adaptInfo, true);

  AMDiS::finalize();
}


