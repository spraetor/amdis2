#include "SubAssembler.h"
#include "Assembler.h"
#include "FiniteElemSpace.h"
#include "Operator.h"
#include "BasisFunction.h"
#include "Mesh.h"
#include "Quadrature.h"
#include "DOFVector.h"

using namespace std;

namespace AMDiS 
{
  SubAssembler::SubAssembler(Operator *op,
                  			     Assembler *assembler,
                  			     Quadrature *quadrat,
                  			     int order, 
                  			     bool optimized,
                  			     FirstOrderType type) 
    : rowFeSpace(assembler->rowFeSpace),
      colFeSpace(assembler->colFeSpace),
      nRow(0),
      nCol(0),
      coordsNumAllocated(0),
      quadrature(quadrat),
      psiFast(NULL),
      phiFast(NULL),
      symmetric(true),
      opt(optimized),
      firstCall(true),
      name("")
  {
    FUNCNAME("SubAssembler::SubAssembler()");

    TEST_EXIT(rowFeSpace)("No row FE space defined!\n");
    TEST_EXIT(colFeSpace)("No col FE space defined!\n");

    const BasisFunction *psi = rowFeSpace->getBasisFcts();
    const BasisFunction *phi = colFeSpace->getBasisFcts();

    nRow = psi->getNumber();
    nCol = phi->getNumber();

    switch (order) {
    case 0:
      terms = op->zeroOrder;
      break;
    case 1:
      if (type == GRD_PHI)
        terms = op->firstOrderGrdPhi;
      else 
        terms = op->firstOrderGrdPsi;
      break;
    case 2:
      terms = op->secondOrder;
      break;
    }

    // If the number of basis functions in row and col fe space are equal,
    // the element matrix may be symmetric if all operator terms are symmetric.
    // If the row and col fe space are different, the element matrix is neither
    // quadratic, and therefore cannot be symmetric.
    if (nRow == nCol) {
      symmetric = true;
      for (size_t i = 0; i < terms.size(); i++) {
      	if (!(terms[i])->isSymmetric()) {
      	  symmetric = false;
      	  break;
      	}
      }  
    } else {
      symmetric = false;
    }

    dim = rowFeSpace->getMesh()->getDim();
  }


  FastQuadrature *SubAssembler::updateFastQuadrature(FastQuadrature *quadFast,
						     const BasisFunction *psi,
						     Flag updateFlag)
  {
    if (!quadFast) {
      quadFast = FastQuadrature::provideFastQuadrature(psi, *quadrature, updateFlag);
    } else {
// #pragma omp critical 
      {
      	if (!quadFast->initialized(updateFlag))
      	  quadFast->init(updateFlag);
      }
    }

    return quadFast;
  }


  void SubAssembler::initImpl(const ElInfo* elInfo,
                  			      const ElInfo*,
                  			      Quadrature *quad)
  {
    // set corrdsAtQPs invalid
    coordsValid = false;

    // set values at QPs invalid
    for (auto& itAny : cachedValuesAtQPs)
      itAny.second->valid = false;

    // set gradients at QPs invalid
    for (auto& itAnyGrd : cachedGradientsAtQPs)
      itAnyGrd.second->valid = false;
    
    // calls initElement of each term
    for (OperatorTerm* term : terms)
      term->initElement(elInfo, this, quad);
  }


  void SubAssembler::getCoordsAtQPs(const ElInfo* elInfo, 
                        				    Quadrature *quad,
                        				    DenseVector<WorldVector<double> > &coordsAtQPs)
  {
    // already calculated for this element ?
    if (coordsValid) {
      coordsAtQPs = cacheCoordsAtQPs;
      return;
    }

    Quadrature *localQuad = quad ? quad : quadrature; 
    const int nPoints = localQuad->getNumPoints();
  
    coordsAtQPs.change_dim(nPoints);
    for (int i = 0; i < nPoints; i++)
      elInfo->coordToWorld(localQuad->getLambda(i), coordsAtQPs[i]);

    // mark values as valid
    coordsValid = true;
    cacheCoordsAtQPs.change_dim(num_rows(coordsAtQPs));
    cacheCoordsAtQPs = coordsAtQPs;
  }

} // end namespace AMDiS
