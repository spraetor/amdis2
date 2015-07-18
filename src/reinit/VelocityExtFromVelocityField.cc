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


#include "VelocityExtFromVelocityField.h"

namespace reinit {

// double NormEps::eps = 0.0;

void 
VelocityExtFromVelocityField::calcVelocityBoundary(DegreeOfFreedom *locInd, 
						   const int indexV)
{
  // ===== Calculate normal velocity in element vertices. =====

  // Get gradient of lSFct on element.
  DimVec<double> lambda(dim, DEFAULT_VALUE, 0.0);
  lambda[0] = 1.0;
  WorldVector<double> elGrd;

  for (int i=0; i<=dim; ++i)
    lSFctVal[i] = (*lSFct)[locInd[i]];

  basFcts->evalGrdUh(lambda, elInfo->getGrdLambda(), lSFctVal, elGrd);

  // Calculate normal velocity.
  double fac = 1.0 / NormEps::calcNormEps(elGrd);
  double sP;
  for (int i=0; i<=dim; ++i) {
    sP = 0.0;
    for (int j=0; j<dim; ++j) {
      sP += (*velField[j])[locInd[i]] * elGrd[j];
    }
    elNormalVel[i] = fac * sP;
  }
    

  // ===== Extend normal velocity from interface to element vertices. =====
  double tempV=0.0;
  for(int i=0; i<=dim; i++)
    {
      tempV += elNormalVel[i]*lamVec[indexV][i];
    }
  (*(velDOF[0]))[locInd[indexV]]=tempV;
}

}
