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


#include "ProblemInterpol.h"
#include "Mesh.h"
#include "Error.h"
#include "SystemVector.h"
#include "AdaptInfo.h"

// TODO: move to deprecate
namespace AMDiS {

  using namespace std;

  ProblemInterpol::ProblemInterpol(const char *nameStr,
				   ProblemStatSeq *spaceProblem,
				   vector<AbstractFunction<double, WorldVector<double> >*> *fct,
				   vector<AbstractFunction<WorldVector<double>, WorldVector<double> >*> *grdFct)
    : ProblemStatSeq(nameStr),
      interpolFct(fct),
      grdInterpolFct(grdFct)
  {
    Flag adoptFlag = INIT_SYSTEM | INIT_MESH | INIT_FE_SPACE;
    Flag initFlag = INIT_ALL & ~adoptFlag & ~INIT_SOLVER & ~INIT_ESTIMATOR;
  
    initialize(initFlag, spaceProblem, adoptFlag);
  }


  void ProblemInterpol::solve(AdaptInfo *adaptInfo, bool, bool) 
  {
    int size = static_cast<int>(meshes.size());
    for (int i = 0; i < size; i++)
      meshes[i]->dofCompress();
    
    solution->interpol(interpolFct);
  }


  void ProblemInterpol::estimate(AdaptInfo *adaptInfo) 
  {
    FUNCNAME("ProblemIterpolVec::estimate()");

    double errMax = 0.0, errSum = 0.0;
    int errorNorm = 0;
    int size = 
      static_cast<int>(interpolFct ? interpolFct->size() : grdInterpolFct->size());
    int relErr = 0;
    Parameters::get(name + "->rel error", relErr);

    if (grdInterpolFct) 
      errorNorm = 1;
    else
      if (interpolFct) 
	errorNorm = 2;

    switch (errorNorm) {
    case 1:
      for (int i = 0; i < size; i++) {
	errSum = Error<double>::H1Err((*(*grdInterpolFct)[i]), 
				      *(solution->getDOFVector(i)), 
				      relErr, &errMax, true, i);
	adaptInfo->setEstSum(errSum, i);
	adaptInfo->setEstMax(errMax, i);
      }
      break;
    case 2:
      for (int i = 0; i < size; i++) {
	errSum = Error<double>::L2Err((*(*interpolFct)[i]), 
				      *(solution->getDOFVector(i)), 
				      relErr, &errMax, true, i);
	adaptInfo->setEstSum(errSum, i);
	adaptInfo->setEstMax(errMax, i);
      }
      break;
    default: 
      ERROR_EXIT("invalid error norm\n");
    }

    MSG("estimate: %e\n", errSum);
  }

}
