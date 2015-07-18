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


#include "Marker.h"

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/MpiHelper.h"
#endif

namespace AMDiS {

  Marker *Marker::createMarker(std::string name, int row) 
  {
    int strategy = 0;
    Parameters::get(name + "->strategy", strategy);
  
    Marker *marker = NULL;

    switch (strategy) {
    case 0: 
      break;
    case 1:
      marker = new GRMarker(name, row);
      break;
    case 2:
      marker = new MSMarker(name, row);
      break;
    case 3:
      marker = new ESMarker(name, row);
      break;
    case 4:
      marker = new GERSMarker(name, row);
      break;
    default: 
      ERROR_EXIT("invalid strategy\n");
    }
    
    return marker;
  }

  
  void Marker::initMarking(AdaptInfo *adaptInfo, Mesh *mesh)
  {
    FUNCNAME_DBG("Marker::initMarking()");

    TEST_EXIT_DBG(adaptInfo)("No AdaptInfo object!\n");

    elMarkRefine = 0;
    elMarkCoarsen = 0;
    estSum = pow(adaptInfo->getEstSum(row == -1 ? 0 : row), p);
    estMax = adaptInfo->getEstMax(row == -1 ? 0 : row);
  }


  void Marker::finishMarking(AdaptInfo *adaptInfo) 
  {
    FUNCNAME("Marker::finishMarking()");
    
    INFO(info, 4)("%d elements marked for refinement\n", elMarkRefine);
    INFO(info, 4)("%d elements marked for coarsening\n", elMarkCoarsen);
  }


  void Marker::markElement(AdaptInfo *adaptInfo, ElInfo *elInfo) 
  {
    Element *el = elInfo->getElement();
    double lError = el->getEstimation(row);

    if (adaptInfo->isRefinementAllowed(row == -1 ? 0 : row) && lError > markRLimit) {
      if (maxRefineLevel == -1 || elInfo->getLevel() < maxRefineLevel)
	setMark(el, adaptInfo->getRefineBisections(row == -1 ? 0 : row));
    } else {
      if (adaptInfo->isCoarseningAllowed(row == -1 ? 0 : row) && lError <= markCLimit) {
	if (minRefineLevel == -1 || elInfo->getLevel() > minRefineLevel) {
 	  if (!el->getElementData()->getElementData(COARSENABLE) || 
 	      lError + el->getCoarseningEstimation(row) <= markCLimit)
 	    setMark(el, -adaptInfo->getCoarseBisections(row == -1 ? 0 : row));	
	}
      }
    }
  }


  Flag Marker::markMesh(AdaptInfo *adaptInfo, Mesh *mesh) 
  {
    FUNCNAME_DBG("Marker::markMesh()");

    TEST_EXIT_DBG(mesh)("No mesh!\n");

    initMarking(adaptInfo, mesh);
   
    if (!adaptInfo->isCoarseningAllowed(row == -1 ? 0 : row) && 
	!adaptInfo->isRefinementAllowed(row == -1 ? 0 : row))
      return 0;    
    
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
    while (elInfo) {
      markElement(adaptInfo, elInfo);
      elInfo = stack.traverseNext(elInfo);
    }
    
    finishMarking(adaptInfo);
    
    Flag markFlag;
    if (elMarkRefine) 
      markFlag = 1;
    if (elMarkCoarsen) 
      markFlag |= 2;

    return markFlag; 
  } 


  void MSMarker::initMarking(AdaptInfo *adaptInfo, Mesh *mesh)
  {
    FUNCNAME("MSMarker::initMarking()");

    Marker::initMarking(adaptInfo, mesh);
    
    double MSGammaP = pow(MSGamma, p);
    double MSGammaCP = pow(MSGammaC, p);
    
    markRLimit = MSGammaP * adaptInfo->getEstMax(row == -1 ? 0 : row);
    markCLimit = MSGammaCP * adaptInfo->getEstMax(row == -1 ? 0 : row);

    MSG("start max_est: %.3le mark_limits: %.3le %.3le\n", 
	adaptInfo->getEstMax(row == -1 ? 0 : row), markRLimit, markCLimit);
  }


  void ESMarker::initMarking(AdaptInfo *adaptInfo, Mesh *mesh) 
  {
    FUNCNAME("ESMarker::initMarking()");
    
    Marker::initMarking(adaptInfo, mesh);
   
    double ESThetaP = pow(ESTheta, p);
    double ESThetaCP = pow(ESThetaC, p);
    double epsP = pow(adaptInfo->getSpaceTolerance(row == -1 ? 0 : row), p);

    int nLeaves = mesh->getNumberOfLeaves();
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(nLeaves);
#endif
    
    markRLimit = ESThetaP * epsP / nLeaves;
    markCLimit = ESThetaCP * epsP / nLeaves;

    INFO(info, 2)("start mark_limits: %.3le %.3le   nt = %d\n", 
		  markRLimit, markCLimit, nLeaves);
  }


  Flag GERSMarker::markMesh(AdaptInfo *adaptInfo, Mesh *mesh) 
  {
    FUNCNAME("GERSMarker::markMesh()");

    initMarking(adaptInfo, mesh);

    if (!adaptInfo->isCoarseningAllowed(row == -1 ? 0 : row) && 
	!adaptInfo->isRefinementAllowed(row == -1 ? 0 : row))
      return 0;    

    GERSSum = 0.0;

    double LTheta = pow(1.0 - GERSThetaStar, p);
    double epsP = pow(adaptInfo->getSpaceTolerance(row == -1 ? 0 : row), p);

    if (estSum < oldErrSum) {
      double improv = estSum / oldErrSum;
      double wanted = 0.8 * epsP / estSum;
      double redfac = std::min((1.0 - wanted) / (1.0 - improv), 1.0);
      redfac = std::max(redfac, 0.0);
      
      if (redfac < 1.0) {
	LTheta *= redfac;
	INFO(info, 1)("GERS: use extrapolated theta_star = %lf\n",
		      pow(LTheta, 1.0 / p));
      }
    }

    oldErrSum = estSum;
    double GERSGamma = 1.0;

    if (adaptInfo->isRefinementAllowed(row == -1 ? 0 : row)) {
      if (LTheta > 0) {
	do {
	  GERSSum = 0.0;	  
	  GERSGamma -= GERSNu;
	  markRLimit = GERSGamma * estMax;

	  TraverseStack stack;
	  ElInfo *elInfo = NULL;
	  elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
	  while (elInfo) {
	    markElementForRefinement(adaptInfo, elInfo);
	    elInfo = stack.traverseNext(elInfo);
	  }
	  
	} while((GERSGamma > 0) && (GERSSum < LTheta * estSum));
      }

      INFO(info, 2)("GERS refinement with gamma = %.3lf\n", GERSGamma);
    }

    if (adaptInfo->isCoarseningAllowed(row == -1 ? 0 : row)) {
      GERSGamma = 0.3;
      LTheta = GERSThetaC * epsP;
      
      do {
	GERSSum = 0.0;
	GERSGamma -= GERSNu;
	markCLimit = GERSGamma * estMax;
	
	TraverseStack stack;
	ElInfo *elInfo = NULL;
	elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
	while (elInfo) {
	  markElementForCoarsening(adaptInfo, elInfo);
	  elInfo = stack.traverseNext(elInfo);
	}

	INFO(info, 6)("coarse loop: gamma = %.3e, sum = %.3e, limit = %.3e\n",
		      GERSGamma, GERSSum, LTheta);
      } while(GERSSum > LTheta);
      
      INFO(info, 2)("GERS coarsening with gamma = %.3lf\n", GERSGamma);      
    }
    
    finishMarking(adaptInfo);
    
    Flag markFlag;
    if (elMarkRefine) 
      markFlag = 1;
    if (elMarkCoarsen) 
      markFlag |= 2;
      
    return(markFlag); 
  }


  void GERSMarker::markElementForRefinement(AdaptInfo *adaptInfo, ElInfo *elInfo) 
  {
    Element *el = elInfo->getElement();
    double lError = el->getEstimation(row);

    if (lError > markRLimit) {
      GERSSum += lError;
      setMark(el, adaptInfo->getRefineBisections(row == -1 ? 0 : row));
    }
  }


  void GERSMarker::markElementForCoarsening(AdaptInfo *adaptInfo, ElInfo *elInfo) 
  {
    Element *el = elInfo->getElement();
    double lError = el->getEstimation(row);
    
    if (el->getMark() <= 0) {
      if (el->getElementData()->getElementData(COARSENABLE))
	lError += el->getCoarseningEstimation(row);      
      
      if (lError <= markCLimit) {
	GERSSum += lError;
	setMark(el, -adaptInfo->getCoarseBisections(row == -1 ? 0 : row));
      } else {
	setMark(el, 0);
      }
    }
  }


}
