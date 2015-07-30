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


#include "ElementLevelSet.h"
#include "ElInfo.h"

namespace reinit {

using namespace AMDiS;
int
ElementLevelSet::createElementLevelSet(const ElInfo *elInfo_,
				       const bool doCalcIntersecPts_)
{
  Element *el = elInfo_->getElement();

  if (elInfo == NULL  ||  el != lastEl) {
    /**
     * Element has changed. New calculation.
     */

    // Set new element.
    lastEl = el;

    // Set information for element and reset data.
    setElement(elInfo_);

    // Calculate level set function values for each vertex of element.
    calculateElementLevelSetVal();

    // Check whether level set function values are not too small.
    checkElementLevelSetVal();

    // Calculate status of element and element vertices.
    calculateElementStatus();

    // Calculate intersection points with zero level set if necessary.
    if (doCalcIntersecPts_ && elStatus == LEVEL_SET_BOUNDARY) {
      calculateIntersecPoints();

      // dim == 3: sort intersection points if necessary
      if (dim == 3 && numIntersecPoints == 4) 
	sortIntersecPoints_4IP3D();	
    }
  }

//   else {
    /**
     * LevelSet-Status for element has already been calculated.
     */
//   }

  return elStatus;
}

void
ElementLevelSet::calculateElementLevelSetVal()
{
  DimVec<double> lambda(dim, 0.0);

  for (int i = 0; i <= dim; i++) {

    lambda[i] = 1.0;
    lSFct->setElInfo(elInfo);
    elVertexLevelSetVec[i] = (*lSFct)((const DimVec<double>) lambda);
    lambda[i] = 0.0;
  }
}

int
ElementLevelSet::calculateElementStatus()
{
  for (int i=0; i<=dim; ++i) {
    if (elVertexLevelSetVec[i] < 0) {
      elVertexStatusVec[i] = LEVEL_SET_INTERIOR;
      numElVertexInterior++;
    }    
    else {
      elVertexStatusVec[i] = LEVEL_SET_EXTERIOR;
      numElVertexExterior++;
    }
  }

  /**
   * Calculate level set status of element.
   */
  if (numElVertexInterior == dim+1) {
    elStatus = LEVEL_SET_INTERIOR;
    return LEVEL_SET_INTERIOR;
  }
  if (numElVertexExterior == dim+1) {
    elStatus = LEVEL_SET_EXTERIOR;
    return LEVEL_SET_EXTERIOR;
  }
  elStatus = LEVEL_SET_BOUNDARY;
  return LEVEL_SET_BOUNDARY;
}

void
ElementLevelSet::calculateIntersecPoints()
{
  //**************************************************************************
  // The level set function is linearly approximated by a hyperplane through 
  // the points of the graph of the level set function in the vertices 
  // of the element.
  // This routine calculates the intersection points of the hyperplane 
  // with the edges of the element.
  //**************************************************************************
  DimVec<double> tempPoint(dim);
  DimVec<double> zeroVec(dim, 0.0);
  int i, j;

  /**
   * Get intersection points (in barycentric coordinates).
   *
   * An edge is intersected if its vertices have level sets with
   * different sign.
   */
  for (i = 0; i <= dim; i++) {
    for (j = i+1; j <= dim; j++) {

      if (elVertexStatusVec[i] * elVertexStatusVec[j] < 0.0) {

	  tempPoint = zeroVec;
	  tempPoint[j] = elVertexLevelSetVec[i] / 
	    (elVertexLevelSetVec[i] - elVertexLevelSetVec[j]);
	  checkIntersecBary(tempPoint[j]);
	  tempPoint[i] = 1 - tempPoint[j];

	  (*elIntersecPoints)[numIntersecPoints] = tempPoint;
	  ++numIntersecPoints;
      }

    } // for(j ...
  } // for(i ...
}

int
ElementLevelSet::checkElementLevelSetVal()
{
  int changed = 0;
  double abs_grd_val = 0.0;
  double abs_min_val = 0.0;
  double comp = 0.0;

  for (int i=0; i<dim; ++i) {
    comp = elVertexLevelSetVec[i]-elVertexLevelSetVec[dim];
    abs_grd_val += comp*comp;
  }
  abs_grd_val = sqrt(abs_grd_val);

  abs_min_val = LS_VAL_TOL*abs_grd_val;
  abs_min_val = (abs_min_val > LS_VAL_MIN) ? abs_min_val : LS_VAL_MIN;

  for (int i=0; i<=dim; ++i) {
    if (fabs(elVertexLevelSetVec[i]) < abs_min_val) {
//       elVertexLevelSetVec[i] = abs_min_val;
      elVertexLevelSetVec[i] = (elVertexLevelSetVec[i] < 0) ?
	-abs_min_val : abs_min_val;
      ++changed;
    }
  }

  return changed;
}

bool
ElementLevelSet::checkIntersecBary(double &bary)
{
  if (bary < SP_BARY_TOL) {
    bary = SP_BARY_TOL;
    return true;
  }
  if (bary > 1-SP_BARY_TOL) {
    bary = 1-SP_BARY_TOL;
    return true;
  }

  return false;
}

void
ElementLevelSet::sortIntersecPoints_4IP3D()
{
  int indexFace1 = 0; 
  int indexFace2 = 0; 
  int indexOpVert = 0;
  DimVec<double> tempPoint(dim);
  int i,j;

  /**
   * Consider the 4 intersection points S0, S1, S2 and S3. If the components 
   * of S0 indexFace1 and indexFace2 are not zero, i.e. S0 lies in 
   * the element faces indexFace1 and indexFace2, the intersection 
   * point with zero components indexFace1 and indexFace2 is
   * opposite to S0 in the intersection plane. Move this vertex to the end
   * of the intersection point list elIntersecPoints.
   */ 

  // Get indexFace1 and indexFace2.
  for (i = 0; i < numIntersecPoints; i++) {
    if (fabs((*elIntersecPoints)[0][i]) > 1.e-15) {
      indexFace1 = i;
      break;
    }
  }
  for (j = i+1; j < numIntersecPoints; j++) {
    if (fabs((*elIntersecPoints)[0][j]) > 1.e-15) {
      indexFace2 = j;
      break;
    }
  }

  // Get index of vertex opposite to S0.
  for (i = 1; i < numIntersecPoints; i++) {
    if (fabs((*elIntersecPoints)[i][indexFace1]) <= 1.e-15  && 
	fabs((*elIntersecPoints)[i][indexFace2]) <= 1.e-15) {
      indexOpVert = i;
      break;
    }
  }

  // Move vertex to the end of \ref elIntersecPoints.
  if (indexOpVert != numIntersecPoints-1) {
    tempPoint = (*elIntersecPoints)[indexOpVert];
    (*elIntersecPoints)[indexOpVert] = (*elIntersecPoints)[numIntersecPoints-1];
    (*elIntersecPoints)[numIntersecPoints-1] = tempPoint;
  }

  return;
}

void 
ElementLevelSet::calcIntersecNormal(WorldVector<double> &normalVec)
{
  FUNCNAME("ElementLevelSet::calcIntersecNormal");
  
  switch(dim) {
  case 2: calcIntersecNormal_2d(normalVec);
    break;
  case 3: calcIntersecNormal_3d(normalVec);
    break;
  default: ERROR_EXIT("illegal dimension !\n");
  }
}

void 
ElementLevelSet::calcIntersecNormal_2d(WorldVector<double> &normalVec)
{
  FUNCNAME("ElementLevelSet::calcIntersecNormal_2d");

  TEST_EXIT(numIntersecPoints == 2)("illegal number of intersection points !\n");

  // ===== Get world coordinates of intersection points. =====
  WorldVector<double> sP1;
  WorldVector<double> sP2;
  elInfo->coordToWorld((*elIntersecPoints)[0], sP1);
  elInfo->coordToWorld((*elIntersecPoints)[1], sP2);

  // ===== Calculate normal vector. =====
  double norm2 = 0.0;
  double val;
  for (int i=0; i<dim; ++i){
    val = sP1[i] - sP2[i];
    norm2 += val*val;
    normalVec[dim-1-i] = val;
  } 
  normalVec[0] *= -1;
  double norm = sqrt(norm2);
  for(int i=0; i<dim; ++i) {
    normalVec[i] = 1/norm * normalVec[i];
  }

  // ===== Set correct orientation (exterior normal vector). =====

  // Calculate barycenter.
  WorldVector<double> baryc;
  for (int i=0; i<dim; ++i) {
    baryc[i] = 0.5*sP1[i] + 0.5*sP2[i];
  }
  double d = sqrt((sP1-sP2)*(sP1-sP2));

  // Barycenter + factor*normalVec.
  double elSize = sqrt(fabs(elInfo->getDet()));
  double factor = 0.01*d/elSize;
  WorldVector<double> tmpPoint;
  int cntr = 0;
  DimVec<double> lambda(dim, NO_INIT);
  
  while (1) {
    ++cntr;

    for(int i=0; i<dim; ++i) {
      tmpPoint[i] = baryc[i] + factor*normalVec[i];
    }

    elInfo->worldToCoord(tmpPoint, &lambda);
    for (int i=0; i<=dim; ++i) {
      if (lambda[i] < 0) {
	factor *= 0.1;

	if (cntr == 10) {
	  WARNING("inefficient normal vector calculation !\n");
	}
	if (cntr == 100) {
	  ERROR_EXIT("infinite loop !\n");
	}

	continue;
      }
    }

    break;
  }

  if (ElementLevelSet::calcLevelSetFct(lambda) < 0) {
    for (int i=0; i<dim; ++i) {
      normalVec[i] *= -1;
    }
  }
}


void 
ElementLevelSet::calcIntersecNormal_3d(WorldVector<double> &normalVec)
{
  FUNCNAME("ElementLevelSet::calcIntersecNormal_3d");

  TEST_EXIT(numIntersecPoints == 3 || numIntersecPoints == 4)("illegal number of intersection points !\n");

  // ===== Get world coordinates of intersection points. =====
  WorldVector<double> sP1;
  WorldVector<double> sP2;
  WorldVector<double> sP3;
  elInfo->coordToWorld((*elIntersecPoints)[0], sP1);
  elInfo->coordToWorld((*elIntersecPoints)[1], sP2);
  elInfo->coordToWorld((*elIntersecPoints)[2], sP3);

  // ===== Calculate normal vector. =====
  WorldVector<double> A = sP2-sP1;
  WorldVector<double> B = sP3-sP1;
  vectorProduct(A, B, normalVec);
  
  double norm = sqrt(normalVec * normalVec);
  for(int i=0; i<dim; ++i) {
    normalVec[i] = 1/(norm+1.e-4) * normalVec[i];
  }

  // ===== Set correct orientation (exterior normal vector). =====

  // Calculate barycenter.
  WorldVector<double> baryc;
  for (int i=0; i<dim; ++i) {
    baryc[i] = 1.0/3*sP1[i] + 1.0/3*sP2[i] + 1.0/3*sP3[i];
  }
  double d = sqrt((sP1-sP2)*(sP1-sP2));

  // Barycenter + factor*normalVec.
  double elSize = pow(fabs(elInfo->getDet()), 1.0/3);
  double factor = 0.01*d/elSize;
  WorldVector<double> tmpPoint;
  int cntr = 0;
  DimVec<double> lambda(dim, NO_INIT);
  
  while (1) {
    ++cntr;

    for(int i=0; i<dim; ++i) {
      tmpPoint[i] = baryc[i] + factor*normalVec[i];
    }

    elInfo->worldToCoord(tmpPoint, &lambda);
    for (int i=0; i<=dim; ++i) {
      if (lambda[i] < 0) {
	factor *= 0.1;

	if (cntr == 10) {
	  WARNING("inefficient normal vector calculation !\n");
	}
	if (cntr == 100) {
	  ERROR_EXIT("infinite loop !\n");
	}

	continue;
      }
    }

    break;
  }

  if (ElementLevelSet::calcLevelSetFct(lambda) < 0) {
    for (int i=0; i<dim; ++i) {
      normalVec[i] *= -1;
    }
  }
}


int 
ElementLevelSet::getVertexPos(const DimVec<double> barCoords)
{
  double vertex_val= 0.0;

  for (int i=0; i<=dim; ++i) {
    if (barCoords[i] > 1-1.e-15) {
      vertex_val = elVertexLevelSetVec[i];
      break;
    }
  }

  if (vertex_val > 0) {
    return LEVEL_SET_EXTERIOR;
  }
  else {
    return LEVEL_SET_INTERIOR;
  }
}

int
ElementLevelSet::getElPos(const DimVec<double> barCoords)
{
  double ls_val = calcLevelSetFct(barCoords);
  if (ls_val > 0) {
    return LEVEL_SET_EXTERIOR;
  }
  else {
    return LEVEL_SET_INTERIOR;
  }
}

}
