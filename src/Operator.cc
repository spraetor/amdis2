#include "Operator.h"
#include "ElInfo.h"
#include "Assembler.h"
#include "FixVec.h"
#include "DOFVector.h"
#include "Quadrature.h"

namespace AMDiS 
{
  Operator::Operator(const FiniteElemSpace *row, const FiniteElemSpace *col)
    : rowFeSpace(row), 
      colFeSpace(col ? col : row),
      fillFlag(Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | 
	       Mesh::FILL_DET | Mesh::FILL_GRD_LAMBDA),
      needDualTraverse(false),
      assembler(NULL),
      uhOld(NULL),
      optimized(true)
  {
    secondOrder.resize(0);
    firstOrderGrdPsi.resize(0);
    firstOrderGrdPhi.resize(0);
    zeroOrder.resize(0);

    nRow = rowFeSpace->getBasisFcts()->getNumber();
    nCol = colFeSpace->getBasisFcts()->getNumber();
  }


  void Operator::setUhOld(const DOFVectorBase<double> *vec)
  {
    uhOld = vec;
    auxFeSpaces.insert(vec->getFeSpace());
  }


  void Operator::getElementMatrix(const ElInfo *elInfo, 
				  ElementMatrix& userMat, 
				  double factor)
  {
    if (!assembler.get())
      initAssembler(NULL, NULL, NULL, NULL);

    assembler.get()->calculateElementMatrix(elInfo, userMat, factor);
  }


  void Operator::getElementMatrix(const ElInfo *rowElInfo, const ElInfo *colElInfo,
				  const ElInfo *smallElInfo, const ElInfo *largeElInfo,
				  bool rowColFeSpaceEqual,
				  ElementMatrix& userMat, 
				  double factor)
  {
    if (!assembler.get())
      initAssembler(NULL, NULL, NULL, NULL);

    assembler.get()->calculateElementMatrix(rowElInfo, colElInfo, 
					    smallElInfo, largeElInfo,
					    rowColFeSpaceEqual,
					    userMat, factor);
  }


  void Operator::getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor)
  {
    if (!assembler.get())
      initAssembler(NULL, NULL, NULL, NULL);

    assembler.get()->calculateElementVector(elInfo, userVec, factor);
  }


  void Operator::getElementVector(const ElInfo *mainElInfo, const ElInfo *auxElInfo,
				  const ElInfo *smallElInfo, const ElInfo *largeElInfo,
				  ElementVector& userVec,
				  double factor)
  {
    if (!assembler.get())
      initAssembler(NULL, NULL, NULL, NULL);

    assembler.get()->calculateElementVector(mainElInfo, auxElInfo, 
					    smallElInfo, largeElInfo,
					    userVec, factor);
  }


  void Operator::initAssembler(Quadrature *quad2,
			       Quadrature *quad1GrdPsi,
			       Quadrature *quad1GrdPhi,
			       Quadrature *quad0) 
  {    
    if (optimized) {
      Assembler *aptr = 
	new OptimizedAssembler(this, quad2, quad1GrdPsi, quad1GrdPhi, 
			       quad0, rowFeSpace, colFeSpace);
      assembler.set(aptr);
    } else {
      Assembler *aptr = 
	new StandardAssembler(this, quad2, quad1GrdPsi, quad1GrdPhi, 
			      quad0, rowFeSpace, colFeSpace);
      assembler.set(aptr);
    }
  }


  int Operator::getQuadratureDegree(int order, FirstOrderType firstOrderType) 
  {
    std::vector<OperatorTerm*>* terms = NULL;

    switch(order) {
    case 0:
      terms = &zeroOrder;
      break;
    case 1:
      if (firstOrderType == GRD_PHI)
	terms = &firstOrderGrdPhi;
      else 
	terms = &firstOrderGrdPsi;
      break;
    case 2:
      terms = &secondOrder;
      break;
    }

    const BasisFunction *psi = rowFeSpace->getBasisFcts();
    const BasisFunction *phi = colFeSpace->getBasisFcts();

    int psiDegree = psi->getDegree();
    int phiDegree = phi->getDegree();
    int maxTermDegree = 0;

    for (OperatorTerm* term : terms)
      maxTermDegree = std::max(maxTermDegree, term->degree);
   
    return psiDegree + phiDegree - order + maxTermDegree;
  }


  void Operator::finishAssembling()
  {
    if (assembler.get())
      assembler.get()->finishAssembling();
  }


  Assembler* Operator::getAssembler() 
  {
    if (!assembler.get())
      initAssembler(NULL, NULL, NULL, NULL);    

    return assembler.get();
  }


  void Operator::weakEvalSecondOrder(const std::vector<WorldVector<double> > &grdUhAtQP,
				     std::vector<WorldVector<double> > &result) const
  {
    for (std::vector<OperatorTerm*>::const_iterator termIt = secondOrder.begin(); 
	 termIt != secondOrder.end(); ++termIt)
      static_cast<SecondOrderTerm*>(*termIt)->weakEval(grdUhAtQP, result);
  }
  
} // end namespace AMDiS
