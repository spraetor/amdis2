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


#include "PenaltyOperator.h"
#include "SurfaceOperator.h"

namespace compositeFEM {

using namespace AMDiS;
using namespace std;

double 
PenaltyOperator::getPenaltyCoeff(const ElInfo *elInfo)
{
  if (penaltyCoeffFlag) {
    eps = pow(fabs(elInfo->getDet()), degree*1.0/dim) * factor;
    return 1.0/eps;
  }
  else {
    eps = 1.0;
    return 1.0;
  }
}

void 
PenaltyOperator::getElementMatrix(const ElInfo *elInfo, 
				  ElementMatrix& userMat, 
				  double factor)
{
  VectorOfFixVecs<DimVec<double> > *intersecPoints = NULL;
  double penaltyCoeff = getPenaltyCoeff(elInfo);

  /**
   * Get element status. Does element lie completely inside the integration 
   * domain, completely outside of the integration domain or is it 
   * intersected by the boundary ?
   */
  elStatus = elLS->createElementLevelSet(elInfo);

  /**
   * element status == completely inside or outside  
   *                                         --->  penalty term does not
   *                                               depend on this element;
   *                                               no integration has to
   *                                               be done
   *                                                       
   * element status == lies on boundary  ---> surface integration on boundary
   */
  if (elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY) {

    /**
     * Get intersection points.
     */
    intersecPoints = elLS->getElIntersecPoints();

    /**
     * Create SurfaceOperator and calculate element matrix.
     */
    if (dim == 3  &&  elLS->getNumElIntersecPoints() == 4) {

      /**
       * Intersection plane must be divided into two simplices.
       *
       * Note: The intersection points S0, S1, S2, S3 are supposed to be 
       *       alligned in such a manner that a line through S1 and S2 
       *       divides the intersection plane.
       */
      
      // Treat first simplex.
      (*tempCoords)[0] = (*intersecPoints)[0];
      (*tempCoords)[1] = (*intersecPoints)[1];
      (*tempCoords)[2] = (*intersecPoints)[2];
      if (!surfaceOp) {
	surfaceOp = new SurfaceOperator(this, (*tempCoords));
      }
      else {
	surfaceOp->adaptSurfaceOperator((*tempCoords));
      }
      surfaceOp->getElementMatrix(elInfo, userMat, factor*penaltyCoeff);
      
      // Treat second simplex.
      (*tempCoords)[0] = (*intersecPoints)[3];
      surfaceOp->adaptSurfaceOperator((*tempCoords));
      surfaceOp->getElementMatrix(elInfo, userMat, factor*penaltyCoeff);
    }

    else {
      
      /**
       * Intersection points form a simplex.
       */      

      // *intersecPoints always has MAX_INTERSECTION_POINTS entries
      (*tempCoords)[0] = (*intersecPoints)[0];
      if (dim > 1)
	(*tempCoords)[1] = (*intersecPoints)[1];
      if (dim == 3) {
	(*tempCoords)[2] = (*intersecPoints)[2];
      }

      if (!surfaceOp) {
	surfaceOp = new SurfaceOperator(this, *tempCoords);
      }
      else {
	surfaceOp->adaptSurfaceOperator(*tempCoords);
      }
      surfaceOp->getElementMatrix(elInfo, userMat, factor*penaltyCoeff);
    }
    
  } // if (elstatus == LEVEL_SET_BOUNDARY)

//   else if (dim == 1 && ElementLevelSet::getNumVertIntPoints() != 0) {

//     // ===== intersection points are element vertices =====
//     DimVec<double> lambda(dim, DEFAULT_VALUE, 0.0);
//     const int *statusVec = ElementLevelSet::getElVertStatusVec();

//     for (int i=0; i<dim; ++i) {

//       if (statusVec[i] == ElementLevelSet::LEVEL_SET_BOUNDARY) {
// 	lambda[i] = 1.0;
// 	(*tempCoords)[0] = lambda;

// 	if (!surfaceOp) {
// 	  surfaceOp = new SurfaceOperator(this, *tempCoords);
// 	}
// 	else {
// 	  surfaceOp->adaptSurfaceOperator(*tempCoords);
// 	}
// 	surfaceOp->getElementMatrix(elInfo, userMat, factor*penaltyCoeff);

// 	lambda[i] = 0.0;
//       }
//     }

//   }

  return;
}

void 
PenaltyOperator::getElementVector(const ElInfo *elInfo, 
				  ElementVector& userVec, 
				  double factor)
{
  VectorOfFixVecs<DimVec<double> > *intersecPoints = NULL;
  double penaltyCoeff = getPenaltyCoeff(elInfo);

  /**
   * Get element status. Does element lie completely inside the integration 
   * domain, completely outside of the integration domain or is it 
   * intersected by the boundary ?
   */
  elStatus = elLS->createElementLevelSet(elInfo);

  /**
   * element status == completely inside or outside  
   *                                        --->  penalty term does not
   *                                              depend on this element;
   *                                              no integration has to
   *                                              be done
   *                                                       
   * element status == lies on boundary  ---> surface integration on boundary
   */
  if (elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY) {

    /**
     * Get intersection points.
     */
    intersecPoints = elLS->getElIntersecPoints();

    /**
     * Create SurfaceOperator and calculate element vector.
     */
    if (dim == 3  &&  elLS->getNumElIntersecPoints() == 4) {

      /**
       * Intersection plane must be divided into two simplices.
       *
       * Note: The intersection points S0, S1, S2, S3 are supposed to be 
       *       alligned in such a manner that a line through S1 and S2 
       *       divides the intersection plane.
       */
      
      // Treat first simplex.
      (*tempCoords)[0] = (*intersecPoints)[0];
      (*tempCoords)[1] = (*intersecPoints)[1];
      (*tempCoords)[2] = (*intersecPoints)[2];
      if (!surfaceOp) {
	surfaceOp = new SurfaceOperator(this, (*tempCoords));
      }
      else {
	surfaceOp->adaptSurfaceOperator((*tempCoords));
      }
      surfaceOp->getElementVector(elInfo, userVec, factor * penaltyCoeff);
      
      // Treat second simplex.
      (*tempCoords)[0] = (*intersecPoints)[3];
      surfaceOp->adaptSurfaceOperator((*tempCoords));
      surfaceOp->getElementVector(elInfo, userVec, factor * penaltyCoeff);
    }

    else {
      
      /**
       * Intersection points form a simplex.
       */

      // *intersecPoints always has MAX_INTERSECTION_POINTS entries
      (*tempCoords)[0] = (*intersecPoints)[0];
      if (dim > 1)
	(*tempCoords)[1] = (*intersecPoints)[1];
      if (dim == 3) {
	(*tempCoords)[2] = (*intersecPoints)[2];
      }

      if (!surfaceOp) {
	surfaceOp = new SurfaceOperator(this, *tempCoords);
      }
      else {
	surfaceOp->adaptSurfaceOperator(*tempCoords);
      }
      surfaceOp->getElementVector(elInfo, userVec, factor*penaltyCoeff);
    }

  } // if (elstatus == LEVEL_SET_BOUNDARY)

//   else if (dim == 1 && ElementLevelSet::getNumVertIntPoints() != 0) {

//     // ===== intersection points are element vertices =====
//     DimVec<double> lambda(dim, DEFAULT_VALUE, 0.0);
//     const int *statusVec = ElementLevelSet::getElVertStatusVec();

//     for (int i=0; i<dim; ++i) {

//       if (statusVec[i] == ElementLevelSet::LEVEL_SET_BOUNDARY) {
// 	lambda[i] = 1.0;
// 	(*tempCoords)[0] = lambda;

// 	if (!surfaceOp) {
// 	  surfaceOp = new SurfaceOperator(this, *tempCoords);
// 	}
// 	else {
// 	  surfaceOp->adaptSurfaceOperator(*tempCoords);
// 	}
// 	surfaceOp->getElementVector(elInfo, userVec, factor*penaltyCoeff);

// 	lambda[i] = 0.0;
//       }
//     }

//   }

  return;
}

}
