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
#include "DOFVector.h"

namespace AMDiS
{

  Estimator::Estimator(std::string name_, int r)
    : name(name_),
      norm(NO_NORM),
      row(r),
      mesh(NULL),
      auxMesh(NULL),
      traverseInfo(0)
  {
    int norm_(norm);
    Parameters::get(name + "->error norm", norm_);
    norm = (Norm)norm_;
  }


  double Estimator::estimate(double ts)
  {
    //     FUNCNAME("Estimator::estimate()");

    bool dualTraverse = false;

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
    ElInfo* elInfo = stack.traverseFirst(mesh, -1, traverseFlag);
    while (elInfo)
    {
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
    while (cont)
    {
      estimateElement(dualElInfo.rowElInfo, &dualElInfo);
      cont = dualTraverse.traverseNext(dualElInfo);
    }
  }
}
