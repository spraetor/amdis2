#include "Marker.h"

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/MpiHelper.h"
#endif

namespace AMDiS
{
  
  Marker::Marker(std::string name_, int row_)
    : name(name_),
      row(row_),
      maximumMarking(false),
      p(2),
      info(10),
      maxRefineLevel(-1),
      minRefineLevel(-1)
  {
    Parameters::get(name + "->p", p);
    Parameters::get(name + "->info", info);
    Parameters::get(name + "->max refinement level", maxRefineLevel);
    Parameters::get(name + "->min refinement level", minRefineLevel);
  }
  
  Marker* Marker::createMarker(std::string name, int row)
  {
    int strategy = 0;
    Parameters::get(name + "->strategy", strategy);

    Marker* marker = NULL;

    switch (strategy)
    {
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

  
  void Marker::setMark(Element* el, char newMark)
  {
    char oldMark = el->getMark();

    if (!maximumMarking || (newMark > oldMark))
    {
      el->setMark(newMark);

      if (oldMark > 0)
      {
	elMarkRefine--;
      }
      else
      {
	if (oldMark < 0)
	  elMarkCoarsen--;
      }

      if (newMark > 0)
      {
	elMarkRefine++;
      }
      else
      {
	if (newMark < 0)
	  elMarkCoarsen++;
      }
    }
  }

  
  void Marker::initMarking(AdaptInfo& adaptInfo, Mesh* mesh)
  {
    FUNCNAME_DBG("Marker::initMarking()");

    TEST_EXIT_DBG(adaptInfo)("No AdaptInfo object!\n");

    elMarkRefine = 0;
    elMarkCoarsen = 0;
    estSum = pow(adaptInfo.getEstSum(row == -1 ? 0 : row), p);
    estMax = adaptInfo.getEstMax(row == -1 ? 0 : row);
  }


  void Marker::finishMarking(AdaptInfo& adaptInfo)
  {
    FUNCNAME("Marker::finishMarking()");

    INFO(info, 4)("%d elements marked for refinement\n", elMarkRefine);
    INFO(info, 4)("%d elements marked for coarsening\n", elMarkCoarsen);
  }


  void Marker::markElement(AdaptInfo& adaptInfo, ElInfo* elInfo)
  {
    Element* el = elInfo->getElement();
    double lError = el->getEstimation(row);

    if (adaptInfo.isRefinementAllowed(row == -1 ? 0 : row) && lError > markRLimit)
    {
      if (maxRefineLevel == -1 || elInfo->getLevel() < maxRefineLevel)
        setMark(el, adaptInfo.getRefineBisections(row == -1 ? 0 : row));
    }
    else
    {
      if (adaptInfo.isCoarseningAllowed(row == -1 ? 0 : row) && lError <= markCLimit)
      {
        if (minRefineLevel == -1 || elInfo->getLevel() > minRefineLevel)
        {
          if (!el->getElementData()->getElementData(COARSENABLE) ||
              lError + el->getCoarseningEstimation(row) <= markCLimit)
            setMark(el, -adaptInfo.getCoarseBisections(row == -1 ? 0 : row));
        }
      }
    }
  }


  Flag Marker::markMesh(AdaptInfo& adaptInfo, Mesh* mesh)
  {
    FUNCNAME_DBG("Marker::markMesh()");

    TEST_EXIT_DBG(mesh)("No mesh!\n");

    initMarking(adaptInfo, mesh);

    if (!adaptInfo.isCoarseningAllowed(row == -1 ? 0 : row) &&
        !adaptInfo.isRefinementAllowed(row == -1 ? 0 : row))
      return 0;

    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
    while (elInfo)
    {
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
  
  
  // ---------------------------------------------------------------------------
  
  
  GRMarker::GRMarker(std::string name_, int row_)
    : Marker(name_, row_)
  {}

  
  void GRMarker::markElement(AdaptInfo& adaptInfo, ElInfo* elInfo)
  {
    Element* el = elInfo->getElement();
    if (adaptInfo.isRefinementAllowed(row == -1 ? 0 : row))
      setMark(el, adaptInfo.getRefineBisections(row == -1 ? 0 : row));
  }
  
  
  // ---------------------------------------------------------------------------
  
  
  MSMarker::MSMarker(std::string name_, int row_)
    : Marker(name_, row_),
      MSGamma(0.5),
      MSGammaC(0.1)
  {
    Parameters::get(name + "->MSGamma", MSGamma);
    Parameters::get(name + "->MSGammaC", MSGammaC);
  }


  void MSMarker::initMarking(AdaptInfo& adaptInfo, Mesh* mesh)
  {
    FUNCNAME("MSMarker::initMarking()");

    Marker::initMarking(adaptInfo, mesh);

    double MSGammaP = pow(MSGamma, p);
    double MSGammaCP = pow(MSGammaC, p);

    markRLimit = MSGammaP * adaptInfo.getEstMax(row == -1 ? 0 : row);
    markCLimit = MSGammaCP * adaptInfo.getEstMax(row == -1 ? 0 : row);

    MSG("start max_est: %.3le mark_limits: %.3le %.3le\n",
        adaptInfo.getEstMax(row == -1 ? 0 : row), markRLimit, markCLimit);
  }
  
  
  // ---------------------------------------------------------------------------

  
  ESMarker::ESMarker(std::string name_, int row_)
    : Marker(name_, row_),
      ESTheta(0.9),
      ESThetaC(0.2)
  {
    Parameters::get(name + "->ESTheta", ESTheta);
    Parameters::get(name + "->ESThetaC", ESThetaC);
  }
  

  void ESMarker::initMarking(AdaptInfo& adaptInfo, Mesh* mesh)
  {
    FUNCNAME("ESMarker::initMarking()");

    Marker::initMarking(adaptInfo, mesh);

    double ESThetaP = pow(ESTheta, p);
    double ESThetaCP = pow(ESThetaC, p);
    double epsP = pow(adaptInfo.getSpaceTolerance(row == -1 ? 0 : row), p);

    int nLeaves = mesh->getNumberOfLeaves();
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalAdd(nLeaves);
#endif

    markRLimit = ESThetaP * epsP / nLeaves;
    markCLimit = ESThetaCP * epsP / nLeaves;

    INFO(info, 2)("start mark_limits: %.3le %.3le   nt = %d\n",
                  markRLimit, markCLimit, nLeaves);
  }
  
  
  // ---------------------------------------------------------------------------

  
  GERSMarker::GERSMarker(std::string name_, int row_)
    : Marker(name_, row_),
      oldErrSum(0.0),
      GERSThetaStar(0.6),
      GERSNu(0.1),
      GERSThetaC(0.1)
  {
    Parameters::get(name + "->GERSThetaStar", GERSThetaStar);
    Parameters::get(name + "->GERSNu", GERSNu);
    Parameters::get(name + "->GERSThetaC", GERSThetaC);
  }
  

  Flag GERSMarker::markMesh(AdaptInfo& adaptInfo, Mesh* mesh)
  {
    FUNCNAME("GERSMarker::markMesh()");

    initMarking(adaptInfo, mesh);

    if (!adaptInfo.isCoarseningAllowed(row == -1 ? 0 : row) &&
        !adaptInfo.isRefinementAllowed(row == -1 ? 0 : row))
      return 0;

    GERSSum = 0.0;

    double LTheta = pow(1.0 - GERSThetaStar, p);
    double epsP = pow(adaptInfo.getSpaceTolerance(row == -1 ? 0 : row), p);

    if (estSum < oldErrSum)
    {
      double improv = estSum / oldErrSum;
      double wanted = 0.8 * epsP / estSum;
      double redfac = std::min((1.0 - wanted) / (1.0 - improv), 1.0);
      redfac = std::max(redfac, 0.0);

      if (redfac < 1.0)
      {
        LTheta *= redfac;
        INFO(info, 1)("GERS: use extrapolated theta_star = %lf\n",
                      pow(LTheta, 1.0 / p));
      }
    }

    oldErrSum = estSum;
    double GERSGamma = 1.0;

    if (adaptInfo.isRefinementAllowed(row == -1 ? 0 : row))
    {
      if (LTheta > 0)
      {
        do
        {
          GERSSum = 0.0;
          GERSGamma -= GERSNu;
          markRLimit = GERSGamma * estMax;

          TraverseStack stack;
          ElInfo* elInfo = NULL;
          elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
          while (elInfo)
          {
            markElementForRefinement(adaptInfo, elInfo);
            elInfo = stack.traverseNext(elInfo);
          }

        }
        while((GERSGamma > 0) && (GERSSum < LTheta * estSum));
      }

      INFO(info, 2)("GERS refinement with gamma = %.3lf\n", GERSGamma);
    }

    if (adaptInfo.isCoarseningAllowed(row == -1 ? 0 : row))
    {
      GERSGamma = 0.3;
      LTheta = GERSThetaC * epsP;

      do
      {
        GERSSum = 0.0;
        GERSGamma -= GERSNu;
        markCLimit = GERSGamma * estMax;

        TraverseStack stack;
        ElInfo* elInfo = NULL;
        elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
        while (elInfo)
        {
          markElementForCoarsening(adaptInfo, elInfo);
          elInfo = stack.traverseNext(elInfo);
        }

        INFO(info, 6)("coarse loop: gamma = %.3e, sum = %.3e, limit = %.3e\n",
                      GERSGamma, GERSSum, LTheta);
      }
      while(GERSSum > LTheta);

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


  void GERSMarker::markElementForRefinement(AdaptInfo& adaptInfo, ElInfo* elInfo)
  {
    Element* el = elInfo->getElement();
    double lError = el->getEstimation(row);

    if (lError > markRLimit)
    {
      GERSSum += lError;
      setMark(el, adaptInfo.getRefineBisections(row == -1 ? 0 : row));
    }
  }


  void GERSMarker::markElementForCoarsening(AdaptInfo& adaptInfo, ElInfo* elInfo)
  {
    Element* el = elInfo->getElement();
    double lError = el->getEstimation(row);

    if (el->getMark() <= 0)
    {
      if (el->getElementData()->getElementData(COARSENABLE))
        lError += el->getCoarseningEstimation(row);

      if (lError <= markCLimit)
      {
        GERSSum += lError;
        setMark(el, -adaptInfo.getCoarseBisections(row == -1 ? 0 : row));
      }
      else
      {
        setMark(el, 0);
      }
    }
  }

} // end namespace AMDiS
