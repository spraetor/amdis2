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


#include "SubElInfo.h"
#include "ElInfo.h"
#include "Mesh.h"

namespace AMDiS {

  SubElInfo::SubElInfo(VectorOfFixVecs<DimVec<double> > *lambda_, 
		       const ElInfo *elInfo_) 
    : elInfo(elInfo_)
  {
    FUNCNAME("SubElInfo::SubElInfo");

    int nPoints = lambda_->getSize();
    int dim = elInfo_->getMesh()->getDim();

    TEST_EXIT(nPoints == dim + 1)
      ("invalid number of vertices of subelement\n");

    FixVec<WorldVector<double>, VERTEX> worldCoords(dim, NO_INIT);

    lambda = new VectorOfFixVecs<DimVec<double> >(dim, dim + 1, NO_INIT);
    for (int i = 0; i <= dim; i++) {
      (*lambda)[i] = (*lambda_)[i];
    }

    /** 
     * Get worldcoordinates of the vertices of subelement in order to
     * calculate the corresponding determinant.
     */  
    for (int i = 0; i < nPoints; i++) {
      elInfo->coordToWorld((const DimVec<double>)((*lambda)[i]), 
			   worldCoords[i]);
    }  

    det = elInfo->calcDet(worldCoords);
  }

}
