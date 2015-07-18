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


#include "Traverse.h"

#include "ElementLevelSet.h"
#include "LevelSetAdaptMesh.h"

namespace compositeFEM {

using namespace AMDiS;

void 
LevelSetAdaptMesh::adaptMesh(AdaptInfo *adaptInfo)
{
  TraverseStack *stack;
  ElInfo *elInfo;

  doRefine = true;
  doCoarsen = true; 

  // === Process refinement/coarsening until prescribed element sizes are
  //     reached.
  while (doRefine || doCoarsen) {
    doRefine = false;
    doCoarsen = false;
    refinementMarkIsSet = false;
    coarseningMarkIsSet = false;
    
    // === Mark elements for refinement/coarsening.
    stack = new TraverseStack;
    elInfo = stack->traverseFirst(mesh, -1, 
				  Mesh::CALL_LEAF_EL | 
				  Mesh::FILL_NEIGH | 
				  Mesh::FILL_BOUND | 
				  Mesh::FILL_DET |
				  Mesh::FILL_COORDS);

    while (elInfo) {

      if (bandAtInterf)
	markElement_Band(elInfo);
      else
	markElement_noBand(elInfo);

      elInfo = stack->traverseNext(elInfo);
    }    
    delete stack;
    
    // === Refine/coarsen mesh.
    if(refinementMarkIsSet) {
      if (prob->refineMesh(adaptInfo) != MESH_REFINED) {
	doRefine = false;
      }
    }
    if (coarseningMarkIsSet) {
      if(prob->coarsenMesh(adaptInfo) != MESH_COARSENED) {
	  doCoarsen = false;
      }
    }
  } // end: while (doRefine || doCoarsen)
}

void 
LevelSetAdaptMesh::markElement_noBand(ElInfo *elInfo)
{
  // === Get position of element with respect to the zero level set.
  int elStatus = elLS->createElementLevelSet(elInfo, false);
      
  double elSize = fabs(elInfo->getDet());

  if (elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY) { 
    markElement(elInfo, elSize, sizeInterf);
    return;
  }

  if (!onlyInterfRef) {
    if (elStatus == ElementLevelSet::LEVEL_SET_INTERIOR) {
      markElement(elInfo, elSize, sizeNegLs);
      return;
    }

    if (elStatus == ElementLevelSet::LEVEL_SET_EXTERIOR) { 
      markElement(elInfo, elSize, sizePosLs);
      return;
    }      
  }
}

void 
LevelSetAdaptMesh::markElement_Band(ElInfo *elInfo)
{
  // === Get position of element with respect to the band around the 
  //     zero level set.
  int elStatusBand = calculateElementStatusBand(elInfo);
      
  double elSize = fabs(elInfo->getDet());

  if (elStatusBand == LS_ADAPT_MESH_IN_BAND) { 
    markElement(elInfo, elSize, sizeInterf);
    return;
  }

  if (!onlyInterfRef) {
    if (elStatusBand == LS_ADAPT_MESH_INTERIOR) {
      markElement(elInfo, elSize, sizeNegLs);
      return;
    }

    if (elStatusBand == LS_ADAPT_MESH_EXTERIOR) { 
      markElement(elInfo, elSize, sizePosLs);
      return;
    }   
  }
}

void 
LevelSetAdaptMesh::markElement(ElInfo *elInfo,
			       double elSize,
			       double sizeWanted)
{
  // === After mesh adaption: sizeWanted is an upper bound for element size,
  //     i.e. if sizeWanted < 2*elSize the element is not coarsened.

  // Mark for refinement.
  if (elSize > sizeWanted) {
    elInfo->getElement()->setMark(1);
    refinementMarkIsSet = true;
    
    // Is further refinement needed?
    if (elSize > 2*sizeWanted) {
      doRefine = true;
    }
  }

  // Mark for coarsening.
  if (2*elSize <= sizeWanted) {
    elInfo->getElement()->setMark(-1);
    coarseningMarkIsSet = true;
    
    // Is further coarsening needed?
    if (4*elSize <= sizeWanted) {
      doCoarsen = true;
    }
  }
}

int
LevelSetAdaptMesh::calculateElementStatusBand(ElInfo *elInfo)
{
  bool inBand = false;
  const double *elVertexLevelSetVec = elLS->getElVertLevelSetVec();

  // === Get position of element with respect to band around the zero 
  //     level set.
  int elStatus = elLS->createElementLevelSet(elInfo, false);

  for (int i=0; i<=dim; ++i) {
    if (fabs(elVertexLevelSetVec[i]) < bandSize) {
      inBand = true;
      break;
    }
  }

  if (inBand || elStatus == ElementLevelSet::LEVEL_SET_BOUNDARY) {
    return LS_ADAPT_MESH_IN_BAND;
  }
  else if (elStatus == ElementLevelSet::LEVEL_SET_INTERIOR) {
    return LS_ADAPT_MESH_INTERIOR;
  }
  else {
    return LS_ADAPT_MESH_EXTERIOR;
  }
}

void 
LevelSetAdaptMesh::getElementSizesFromInit(const std::string probName,
					   double &sizeInterf_,
					   double &sizeNegLs_,
					   double &sizePosLs_,
					   double &bandSize_)
{
  FUNCNAME("LevelSetAdaptMesh::getElementSizesFromInit()");

  // === Read parameter from file. ===
  int dim;
  Parameters::get(probName + "->dim", dim);

  std::string meshName;
  Parameters::get(probName + "->mesh", meshName);

  int globalRef = 0;
  Parameters::get(meshName + "->global refinements", globalRef);

  int refInterf = 0;
  Parameters::get("LevelSetAdaptMesh->additional refinements at interface", 
		  refInterf);

  int refNegLs = 0;
  Parameters::get("LevelSetAdaptMesh->additional refinements in negative domain", 
		  refNegLs);

  int refPosLs = 0;
  Parameters::get("LevelSetAdaptMesh->additional refinements in positive domain", 
		  refPosLs);

  int numElBand = 0;
  Parameters::get("LevelSetAdaptMesh->number of elements for band", 
		  numElBand);

  double macroElSize = -1.0;
  Parameters::get("LevelSetAdaptMesh->macro element size", 
		  macroElSize);

  // === Calculate element sizes. ===
  double fac = pow(0.5, 1.0/dim);

  sizeInterf_ = macroElSize * pow(fac, globalRef + refInterf);
  sizeNegLs_ = macroElSize * pow(fac, globalRef + refNegLs);
  sizePosLs_ = macroElSize * pow(fac, globalRef + refPosLs);

  bandSize_ = numElBand * sizeInterf_;

  sizeInterf_ = pow(sizeInterf_, dim);
  sizeNegLs_ = pow(sizeNegLs_, dim);
  sizePosLs_ = pow(sizePosLs_, dim);


  TEST_EXIT(sizeInterf_ >= 0)("illegal sizeInterf !\n");
  TEST_EXIT(sizeNegLs_ >= 0)("illegal sizeNegLs !\n");
  TEST_EXIT(sizePosLs_ >= 0)("illegal sizePosLs !\n");
  TEST_EXIT(bandSize_ >= 0)("illegal bandSize !\n");
}

}
