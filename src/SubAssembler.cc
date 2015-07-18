/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/


#include "SubAssembler.h"
#include "Assembler.h"
#include "FiniteElemSpace.h"
#include "Operator.h"
#include "BasisFunction.h"
#include "Mesh.h"
#include "Quadrature.h"
#include "DOFVector.h"

using namespace std;

namespace AMDiS {

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


  void SubAssembler::initImpl(const ElInfo* smallElInfo,
			      const ElInfo *largeElInfo,
			      Quadrature *quad)
  {
    // set corrdsAtQPs invalid
    coordsValid = false;

    // set values at QPs invalid
    map<const void*, ValuesAtQPs*>::iterator itAny;
    for (itAny = cachedValuesAtQPs.begin(); itAny != cachedValuesAtQPs.end(); ++itAny)
      (*itAny).second->valid = false;

    // set gradients at QPs invalid
    map<const void*, ValuesAtQPs*>::iterator itAnyGrd;
    for (itAnyGrd = cachedGradientsAtQPs.begin(); itAnyGrd != cachedGradientsAtQPs.end(); ++itAnyGrd)
      (*itAnyGrd).second->valid = false;
    
    // calls initElement of each term
    for (vector<OperatorTerm*>::iterator it = terms.begin(); 
	 it != terms.end(); ++it) {
      if (largeElInfo == NULL || smallElInfo == largeElInfo)
	(*it)->initElement(smallElInfo, this, quad);
      else
	(*it)->initElement(smallElInfo, largeElInfo, this, quad);      
    }
  }


  void SubAssembler::getCoordsAtQPs(const ElInfo* elInfo, 
				    Quadrature *quad,
				    mtl::dense_vector<WorldVector<double> > &coordsAtQPs)
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

}
