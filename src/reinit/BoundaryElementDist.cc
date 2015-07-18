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


#include "BoundaryElementDist.h"

namespace reinit {

using namespace AMDiS;

void 
BoundaryElementDist::calcNormal(
	             const FixVec<WorldVector<double>, DIMEN> &vecs, 
		     WorldVector<double> &normalVec)
{
  FUNCNAME("BoundaryElementDist::calcNormal");
  
  switch(dim) {
  case 2: calcNormal_2d(vecs, normalVec);
    break;
  case 3: calcNormal_3d(vecs, normalVec);
    break;
  default: ERROR_EXIT("illegal dimension !\n");
  }
}

void 
BoundaryElementDist::calcNormal_2d(
		     const FixVec<WorldVector<double>, DIMEN> &vecs, 
		     WorldVector<double> &normalVec)
{
  FUNCNAME("BoundaryElementDist::calcNormal_2d");
  
  TEST_EXIT(vecs.getSize() == dim)("illegal number of world vectors in vecs !\n");
  
  double norm2 = 0.0;
  double val;
  for (int i=0; i<dim; ++i){
    val = vecs[0][i] - vecs[1][i];
    norm2 += val*val;
    normalVec[dim-1-i] = val;
  } 
  normalVec[0] *= -1;
  double norm = sqrt(norm2);
  for(int i=0; i<dim; ++i) {
    normalVec[i] = 1/norm * normalVec[i];
  }
}

void 
BoundaryElementDist::calcNormal_3d(
		     const FixVec<WorldVector<double>, DIMEN> &vecs, 
		     WorldVector<double> &normalVec)
{
  FUNCNAME("BoundaryElementDist::calcNormal_3d");
  
  TEST_EXIT(vecs.getSize() == dim)("illegal number of world vectors in vecs !\n");
  
  WorldVector<double> A = vecs[1]-vecs[0];
  WorldVector<double> B = vecs[2]-vecs[0];
  vectorProduct(A, B, normalVec);
  
  double norm = sqrt(normalVec * normalVec);
  for(int i=0; i<dim; ++i) {
    normalVec[i] = 1/norm * normalVec[i];
  }
}

}
