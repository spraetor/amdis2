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


#include "HL_SignedDistTraverse.h"
#include "VelocityExtFromVelocityField.h"

namespace reinit {

void HL_SignedDistTraverse::initializeBoundary() 
{
  FUNCNAME("HL_SignedDistTraverse::initializeBoundary()");
 
  // ===== All non-boundary vertices are initialized with "infinity". =====
  sD_DOF->set(inftyValue);
  
  // ===== Traverse mesh and initialize boundary elements. =====
  TraverseStack stack;
  FixVec<double, VERTEX> distVec(dim, NO_INIT);
  int elStatus;
  
  const int nBasFcts = feSpace->getBasisFcts()->getNumber();
  if (static_cast<int>(locInd.size()) < nBasFcts)
    locInd.resize(nBasFcts);
  
  ElInfo *elInfo;
  if (velExt && velExtType.isSet(VEL_EXT_FROM_VEL_FIELD)) {
    elInfo = stack.traverseFirst(feSpace->getMesh(),
				 -1, 
				 Mesh::CALL_LEAF_EL | 
				 Mesh::FILL_BOUND |
				 Mesh::FILL_COORDS |
				 Mesh::FILL_GRD_LAMBDA);
  } else {
    elInfo = stack.traverseFirst(feSpace->getMesh(),
				 -1, 
				 Mesh::CALL_LEAF_EL | 
				 Mesh::FILL_BOUND |
				 Mesh::FILL_COORDS);
  }

  while(elInfo) {

    // Set elInfo in case velocity extension from velocity field is used.
    if (velExt && velExtType.isSet(VEL_EXT_FROM_VEL_FIELD))
      ((VelocityExtFromVelocityField *)velExt)->setElInfo(elInfo);    
    
    // Get local indices of vertices.
    feSpace->getBasisFcts()->getLocalIndices(const_cast<Element*>(elInfo->getElement()),
					     const_cast<DOFAdmin*>(feSpace->getAdmin()),
					     locInd); 
    
    // Get element status.
    elStatus = elLS->createElementLevelSet(elInfo);
    
    // Is element cut by the interface ?
    if (elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY) {
      
      // Reset element distance vector.
      for (int i=0; i <= dim; i++)
	distVec[i] = inftyValue;      
      
      // Mark all vertices as boundary vertices.
      for (int i = 0; i <= dim; i++)
	(*bound_DOF)[locInd[i]] = 1.0;      
      
      // Calculate distance for all vertices.
      if (bndElDist->calcDistOnBoundaryElement(elInfo, distVec) !=
	  ElementLevelSet::LEVEL_SET_BOUNDARY) {
	ERROR_EXIT("error in distance calculation !\n");
      } else {
	
	// If distance is smaller, correct to new distance.
	for (int i = 0; i <= dim; i++) {
	  
	  // --> for test purposes:
	  if (distVec[i] > 1000) {
	    MSG("Element %d  Knoten %d  keine Randwertinitialisierung!\n",
		elInfo->getElement()->getIndex(), i);
	  }
		  
	  // --> end: for test purposes	  
	  if ((*sD_DOF)[locInd[i]] > distVec[i]) {
	    (*sD_DOF)[locInd[i]] = distVec[i];
	    //If Distance is corrected, calculate new velocity.
	    if (velExt != NULL)	      
	      velExt->calcVelocityBoundary(&(locInd[0]), i);	      
	  }
	}
      }
    }  // end of: elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY
    
    elInfo = stack.traverseNext(elInfo);
  }  // end of: mesh traverse 
}


void HL_SignedDistTraverse::HL_updateIteration()
{
  FUNCNAME("HL_SignedDistTraverse::HL_updateIteration()");

  // ===== Create DOF vector for the last iteration step. =====
  if (sDOld_DOF)
    delete sDOld_DOF;

  sDOld_DOF = new DOFVector<double>(feSpace, "sDOld_DOF");
  sDOld_DOF->copy(const_cast<DOFVector<double> &>(*sD_DOF));
  
  // ===== Gauss-Seidel or Jacobi iteration ? =====
  if (GaussSeidelFlag)
    update_DOF = sD_DOF;  
  else
    update_DOF = sDOld_DOF;  
  
  // ===== Iteration loop: proceed until tolerance is reached. =====
  TraverseStack stack;
  ElInfo *elInfo;
  tol_reached = false;
  int itCntr = 0;
  while (!tol_reached && itCntr != maxIt) {   
    ++itCntr;
    tol_reached = true;
    
    // ===== Traverse mesh: perform Hopf-Lax element update on each element. =====
    elInfo = 
      stack.traverseFirst(feSpace->getMesh(), -1, 
			  Mesh::CALL_LEAF_EL | Mesh::FILL_BOUND | Mesh::FILL_COORDS);
    while (elInfo) {
      HL_elementUpdate(elInfo);
      elInfo = stack.traverseNext(elInfo);
    }
    
    // ===== Is tolerance reached ? =====
    tol_reached = checkTol();
    
    sDOld_DOF->copy(const_cast<DOFVector<double> &>(*sD_DOF));
  }
  

  MSG("Calculation of signed distance function via mesh traverse iteration:\n");
  if (GaussSeidelFlag)
    MSG("\tGauss-Seidel iteration\n");
  else
    MSG("\tJacobi iteration\n");
  
  MSG("\tnumber of iterations needed: %d\n", itCntr);
}


void HL_SignedDistTraverse::HL_elementUpdate(ElInfo *elInfo) 
{
  // ===== Get global indices of vertices of element. =====
  const int nBasFcts = feSpace->getBasisFcts()->getNumber();
  if (static_cast<int>(locInd.size()) < nBasFcts)
    locInd.resize(nBasFcts);
  
  feSpace->getBasisFcts()->getLocalIndices(const_cast<Element *>(elInfo->getElement()),
					   const_cast<DOFAdmin *>(feSpace->getAdmin()),
					   locInd);
  
  // ===== Hopf-Lax element update for each vertex of element. =====
  for (int i = 0; i <= dim; i++) {
    
    // ===== Calculate update for non-boundary vertex. =====
    if ((*bound_DOF)[locInd[i]] < 1.e-15) {
      //save permutation of vertexes for calculation of the velocity
      if (velExt != NULL)
	velExt->setPermutation(i, 1);
	
      double update = calcElementUpdate(elInfo, i, &(locInd[0]));
      // ---> for test purposes: count number of calculated updates
      ++calcUpdate_Cntr;
      // ---> end: for test purposes
      
      // Calculates minimum of all element updates on elements
      // containing vertex i.
      if (update < (*sD_DOF)[locInd[i]]) {
	(*sD_DOF)[locInd[i]] = update;

	// If Distance is corrected, calculate new velocity.
	if(velExt != NULL)
	  velExt->calcVelocity(&(locInd[0]), i);
	
	// ---> for test purposes: count number of calculated updates
	++setUpdate_Cntr;
	// ---> end: for test purposes
      }
    }
  } 
}


double HL_SignedDistTraverse::calcElementUpdate(ElInfo *elInfo, 
						int nXh,
						const DegreeOfFreedom *locInd) 
{
  // ===== Get local indices of element vertices (except xh). =====
  int nWh = 0;  
  int nYh = (nXh + 1) % (dim+1);
  int nZh = (nXh + 2) % (dim+1);
  if (dim == 3) 
    nWh = (nXh + 3) % (dim+1);
  
  // ===== Get world coordinates of vertices of element and their values
  //       of uh. 
  //       The coordinates of the vertex the update is calculated for
  //       are stored at the end of the vector elVert. =====
  switch (dim) {
  case 2: 
    elVert[0] = &(elInfo->getCoord(nYh));
    elVert[1] = &(elInfo->getCoord(nZh));
    elVert[2] = &(elInfo->getCoord(nXh));
    
    uhVal[0] = (*update_DOF)[locInd[nYh]];
    uhVal[1] = (*update_DOF)[locInd[nZh]];
    uhVal[2] = (*update_DOF)[locInd[nXh]];
    
    break;
  case 3: 
    elVert[0] = &(elInfo->getCoord(nYh));
    elVert[1] = &(elInfo->getCoord(nZh));
    elVert[2] = &(elInfo->getCoord(nWh));
    elVert[3] = &(elInfo->getCoord(nXh));

    uhVal[0] = (*update_DOF)[locInd[nYh]];
    uhVal[1] = (*update_DOF)[locInd[nZh]];
    uhVal[2] = (*update_DOF)[locInd[nWh]];
    uhVal[3] = (*update_DOF)[locInd[nXh]];
    
    break;
  default: ERROR_EXIT("illegal dimension !\n");
    break;
  }
  
  // ===== Calculate Hopf-Lax element update for vertex. =====
  return elUpdate->calcElementUpdate(elVert, uhVal);
}


bool HL_SignedDistTraverse::checkTol()
{
  DOFVector<double>::Iterator it_sD(sD_DOF, USED_DOFS);
  DOFVector<double>::Iterator it_sDOld(sDOld_DOF, USED_DOFS);
  
  for (it_sD.reset(), it_sDOld.reset(); !it_sD.end(); ++it_sD, ++it_sDOld)
    if ((*it_sDOld) - (*it_sD) > tol || (*it_sDOld) - (*it_sD) < 0)
      return false;
  
  return true;
}

}
