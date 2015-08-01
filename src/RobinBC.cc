#include "RobinBC.h"
#include "Assembler.h"
#include "DOFVector.h"
#include "DOFMatrix.h"
#include "SurfaceOperator.h"
#include "est/Estimator.h"

namespace AMDiS 
{

  void RobinBC::fillBoundaryCondition(DOFVectorBase<double>* vector, 
                        				      ElInfo* elInfo,
                        				      const DegreeOfFreedom* dofIndices,
                        				      const BoundaryType* localBound,
                        				      int nBasFcts)
  {
    FUNCNAME_DBG("RobinBC::fillBoundaryCondition()");
    TEST_EXIT_DBG(vector->getFeSpace() == rowFeSpace)("invalid row fe space\n");

    int dim = elInfo->getMesh()->getDim();

    if (neumannOperators) {
      for (int i = 0; i < dim + 1; i++)
      	if (elInfo->getBoundary(i) == boundaryType)
      	  vector->assemble(1.0, elInfo, localBound, (*neumannOperators)[i]);
    }
  }


  void RobinBC::fillBoundaryCondition(DOFMatrix* matrix,
                        				      ElInfo* elInfo,
                        				      const DegreeOfFreedom* dofIndices,
                        				      const BoundaryType* localBound,
                        				      int nBasFcts) 
  {
    if (robinOperators) {
      int dim = elInfo->getMesh()->getDim();

      for (int i = 0; i < dim + 1; i++)
      	if (elInfo->getBoundary(i) == boundaryType)
      	  matrix->assemble(1.0, elInfo, localBound, (*robinOperators)[i]);
    }
  }
  

  double RobinBC::boundResidual(ElInfo *elInfo,
                        				DOFMatrix *matrix,
                        				const DOFVectorBase<double> *dv)
  {
    FUNCNAME("RobinBC::fillBoundaryCondition()");
    TEST_EXIT(matrix->getRowFeSpace() == rowFeSpace)("invalid row fe space\n");
    TEST_EXIT(matrix->getColFeSpace() == colFeSpace)("invalid col fe space\n");

    int dim = elInfo->getMesh()->getDim();
    DimVec<double>  lambda(dim);
    double n_A_grdUh, val = 0.0;
    WorldVector<double> normal;
    const DimVec<WorldVector<double> > &Lambda = elInfo->getGrdLambda();
    double det = elInfo->getDet();
    bool neumannQuad = false;
    const BasisFunction *basFcts = dv->getFeSpace()->getBasisFcts();

    TEST_EXIT(basFcts == rowFeSpace->getBasisFcts())("invalid basFcts\n");

    ElementVector uhEl(basFcts->getNumber());
    dv->getLocalVector(elInfo->getElement(), uhEl);

    TEST_EXIT(neumannOperators || robinOperators)
      ("neither neumann nor robin operators set!\n");

    if (!robinOperators)
      neumannQuad = true;
    else {
      if (neumannOperators) {
      	if ((*neumannOperators)[0]->getAssembler()->
      	    getZeroOrderAssembler()->getQuadrature()->getNumPoints() > 
      	    (*robinOperators)[0]->getAssembler()->
      	    getZeroOrderAssembler()->getQuadrature()->getNumPoints()) 
      	  neumannQuad = true;
      }
    }

    std::vector<Operator*>::iterator op;
    for (op = matrix->getOperatorsBegin(); op != matrix->getOperatorsEnd(); ++op)
      (*op)->getAssembler()->initElement(elInfo);        

    for (int face = 0; face < dim + 1; face++) {
      elInfo->getNormal(face, normal);

      Quadrature *quadrature = neumannQuad ? 
      	(*neumannOperators)[face]->getAssembler()->
      	   getZeroOrderAssembler()->getQuadrature() :
      	(*robinOperators)[face]->getAssembler()->
      	   getZeroOrderAssembler()->getQuadrature();
      
      if (elInfo->getBoundary(face) == boundaryType) {
      	(*neumannOperators)[face]->getAssembler()->
      	  getZeroOrderAssembler()->initElement(elInfo);
      
      	int nPoints = quadrature->getNumPoints();
      	DenseVector<double> uhAtQp(nPoints);
      	DenseVector<WorldVector<double> > grdUhAtQp;
      	DenseVector<WorldMatrix<double> > D2UhAtQp;
      	dv->getVecAtQPs(elInfo, quadrature, NULL, uhAtQp);
        
      	ElementVector f(nPoints);
      	set_to_zero(f);
      
      	if (robinOperators)
      	  (*robinOperators)[face]->evalZeroOrder(nPoints, uhAtQp, grdUhAtQp,  D2UhAtQp, f, -1.0);
      	
      	std::vector<WorldVector<double> > grdUh(nPoints);
      	std::vector<WorldVector<double> > A_grdUh(nPoints);
      
      	for (int iq = 0; iq < nPoints; iq++) {
      	  A_grdUh[iq].set(0.0);	
      	  lambda = quadrature->getLambda(iq);
      	  basFcts->evalGrdUh(lambda, Lambda, uhEl, grdUh[iq]);
      	}
      	
      	for (op = matrix->getOperatorsBegin(); op != matrix->getOperatorsEnd(); ++op)
      	  (*op)->weakEvalSecondOrder(grdUh, A_grdUh);		
      
      	if (neumannOperators)
      	  (*neumannOperators)[face]->getC(elInfo, nPoints, f);
      
      	val = 0.0;
      	for (int iq = 0; iq < nPoints; iq++) {
      	  n_A_grdUh = (normal * A_grdUh[iq]) - f[iq]; 
      	  val += quadrature->getWeight(iq) * math::sqr(n_A_grdUh);
      	}
      }
    }

    return det * val;
  }

} // end namespace AMDiS
