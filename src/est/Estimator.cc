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


#include "Estimator.h"
#include "Traverse.h"
#include "Initfile.h"
#include "DualTraverse.h"

namespace AMDiS {

  Estimator::Estimator(std::string name_, int r) 
    : name(name_),
      norm(NO_NORM),
      row(r),
      mesh(NULL),
      auxMesh(NULL),
      traverseInfo(0)
  {
    Parameters::get(name + "->error norm", norm);
  }


  double Estimator::estimate(double ts)
  {
//     FUNCNAME("Estimator::estimate()");
    
    bool dualTraverse = false;

/*    
    for (unsigned int i = 0; i < matrix.size(); i++) {
      TEST_EXIT(traverseInfo.getStatus(row, i) != SingleComponentInfo::DIF_SPACES_WITH_DIF_AUX)
	("Not yet implemented!\n");

      if (traverseInfo.getStatus(row, i) == SingleComponentInfo::EQ_SPACES_WITH_DIF_AUX ||
	  traverseInfo.getStatus(row, i) == SingleComponentInfo::DIF_SPACES_NO_AUX ||
	  traverseInfo.getStatus(row, i) == SingleComponentInfo::DIF_SPACES_WITH_AUX)
	dualTraverse = true;
#ifndef NDEBUG
      MSG("traverseInfo = %d, dualTraverse = %d\n", traverseInfo.getStatus(row, i), int(dualTraverse));
#endif
    }

    if (!dualTraverse) {
      mesh = uh[row == -1 ? 0 : row]->getFeSpace()->getMesh();
      auxMesh = NULL;
    } else {
      const FiniteElemSpace *mainFeSpace = traverseInfo.getRowFeSpace(row);
      const FiniteElemSpace *auxFeSpace = traverseInfo.getNonRowFeSpace(row);

      TEST_EXIT(mainFeSpace)("No main FE space!\n");
      TEST_EXIT(auxFeSpace)("No aux FE space!\n"); 

      mesh = mainFeSpace->getMesh();
      auxMesh = auxFeSpace->getMesh();

      TEST_EXIT_DBG(mainFeSpace->getBasisFcts()->getDegree() ==
		    auxFeSpace->getBasisFcts()->getDegree())
	("Mh, do you really want to do this? Think about it ...\n");
    }
*/    

     mesh = uh[row == -1 ? 0 : row]->getFeSpace()->getMesh();
     auxMesh = NULL;

    init(ts);

    if (!dualTraverse)
      singleMeshTraverse();
    else
      dualMeshTraverse();

    exit();

    return est_sum;
  }


  void Estimator::singleMeshTraverse()
  {
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1, traverseFlag);
    while (elInfo) {
      estimateElement(elInfo);
      elInfo = stack.traverseNext(elInfo);
    }  
  }


  void Estimator::dualMeshTraverse()
  {
    DualTraverse dualTraverse;
    DualElInfo dualElInfo;

    bool cont = dualTraverse.traverseFirst(mesh, auxMesh, -1, -1, 
					   traverseFlag, traverseFlag,
					   dualElInfo);
    while (cont) {
      estimateElement(dualElInfo.rowElInfo, &dualElInfo);      
      cont = dualTraverse.traverseNext(dualElInfo);
    }
  }
}
