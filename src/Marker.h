/** \file Marker.h */

#pragma once

#include "Element.h"
#include "ElInfo.h"
#include "Flag.h"
#include "Mesh.h"
#include "AdaptInfo.h"
#include "Traverse.h"
#include "MatrixVector.h"
#include "Initfile.h"

namespace AMDiS 
{

  /**
   * \ingroup Adaption
   * 
   * \brief
   * Base class for all scalar markers.
   */
  class Marker
  {
  public:
    /// Constructor.
    Marker() {}

    /// Constructor.
    Marker(std::string name_, int row_)
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

    /// destructor
    virtual ~Marker() {}

    /// Marks element with newMark. If \ref maximumMarking is set, the element
    /// is marked only if the new mark is bigger than the old one. The return
    /// value specifies whether the element has been marked, or not.
    void setMark(Element *el, char newMark) 
    {
      // TODO: move implementation ouside of class
      char oldMark = el->getMark();

      if (!maximumMarking || (newMark > oldMark)) {
      	el->setMark(newMark);
      
      	if (oldMark > 0) {
      	  elMarkRefine--; 
      	} else {
      	  if (oldMark < 0) 
      	    elMarkCoarsen--;
      	}
      
      	if (newMark > 0) {
      	  elMarkRefine++; 
      	} else { 
      	  if (newMark < 0) 
      	    elMarkCoarsen++;
      	}
      }
    }

    /// Can be used by sub classes. Called before traversal.
    virtual void initMarking(AdaptInfo *adaptInfo, Mesh *mesh);

    /// Can be used by sub classes. Called after traversal.
    virtual void finishMarking(AdaptInfo *adaptInfo);

    /// Marks one element.
    virtual void markElement(AdaptInfo *adaptInfo, ElInfo *elInfo);
  
    /// Marking of the mesh.
    virtual Flag markMesh(AdaptInfo *adaptInfo, Mesh *mesh);

    /// Sets \ref maximumMarking.
    void setMaximumMarking(bool maxMark) 
    {
      maximumMarking = maxMark;
    }

    int getElMarkRefine() 
    { 
      return elMarkRefine; 
    }

    int getElMarkCoarsen() 
    { 
      return elMarkCoarsen; 
    }

    /// Returns \ref name of the Marker
    std::string getName() const 
    { 
      return name; 
    }

    /// Creates a scalr marker depending on the strategy set in parameters.
    static Marker *createMarker(std::string name, int row_);

  protected:
    /// Name of the scalar marker.
    std::string name;

    /// Equal to -1 for scalar problems. Component number if marker is
    /// part of a vector valued marker.
    int row;

    /// estimation sum
    double estSum;

    /// estmation maximum
    double estMax;

    /// If true, elements are marked only if the new value is bigger than
    /// the current marking.
    bool maximumMarking;

    /// Lower limit for error estimation, from which an element is marked for
    /// refinement
    double markRLimit;

    /// Upper limit for error estimation, from which an element is marked for 
    /// coarsening
    double markCLimit;

    /// power in estimator norm
    double p;

    /// Info level.
    int info;

    /// Counter for elements marked for refinement
    int elMarkRefine;

    /// Counter for elements marked for coarsening
    int elMarkCoarsen;

    /// Maximal level of all elements.
    int maxRefineLevel;

    /// Minimal level of all elements.
    int minRefineLevel;
  };


  /**
   * \ingroup Adaption
   * 
   * \brief
   * Global refinement.
   */
  class GRMarker : public Marker
  {
  public:
    /// Constructor.
    GRMarker(std::string name_, int row_) 
      : Marker(name_, row_) 
    {}

    /// Implementation of Marker::markElement().
    virtual void markElement(AdaptInfo *adaptInfo, ElInfo *elInfo) 
    {
      Element *el = elInfo->getElement();
      if (adaptInfo->isRefinementAllowed(row == -1 ? 0 : row))
        setMark(el, adaptInfo->getRefineBisections(row == -1 ? 0 : row));
    }
  };


  /**
   * \ingroup Adaption
   * 
   * \brief
   * Maximum strategy.
   */
  class MSMarker : public Marker
  {
  public:
    /// Constructor.
    MSMarker(std::string name_, int row_) 
      : Marker(name_, row_),
      	MSGamma(0.5),
      	MSGammaC(0.1)
    {
      Parameters::get(name + "->MSGamma", MSGamma);
      Parameters::get(name + "->MSGammaC", MSGammaC);
    }

    /// Implementation of MarkScal::initMarking().
    void initMarking(AdaptInfo *adaptInfo, Mesh *mesh);
  
  protected:
    /// Marking parameter.
    double MSGamma;

    /// Marking parameter.
    double MSGammaC;
  };


  /**
   * \ingroup Adaption
   * 
   * \brief
   * Equidistribution strategy
   */
  class ESMarker : public Marker
  {
  public:
    /// Constructor.
    ESMarker(std::string name_, int row_) 
      : Marker(name_, row_),
      	ESTheta(0.9),
      	ESThetaC(0.2)
    {
      Parameters::get(name + "->ESTheta", ESTheta);
      Parameters::get(name + "->ESThetaC", ESThetaC);
    }

    /// Implementation of MarkScal::initMarking().
    virtual void initMarking(AdaptInfo *adaptInfo, Mesh *mesh);
  
  protected:
    /// Marking parameter.
    double ESTheta;

    /// Marking parameter.
    double ESThetaC;
  };


  /**
   * \ingroup Adaption
   * 
   * \brief
   * Guaranteed error reduction strategy 
   */
  class GERSMarker : public Marker
  {
  public:
    /// Constructor.
    GERSMarker(std::string name_, int row_) 
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

    /// Implementation of Marker::markMesh().
    virtual Flag markMesh(AdaptInfo *adaptInfo, Mesh *mesh);

  protected:
    /// Refinement marking function.
    void markElementForRefinement(AdaptInfo *adaptInfo, ElInfo *elInfo);

    /// Coarsening marking function.
    void markElementForCoarsening(AdaptInfo *adaptInfo, ElInfo *elInfo);

  protected:
    /// Marking parameter.
    double GERSSum;

    /// Marking parameter.
    double oldErrSum;

    /// Marking parameter.
    double GERSThetaStar;

    /// Marking parameter.
    double GERSNu;

    /// Marking parameter.
    double GERSThetaC;
  };

} // end namespace AMDiS
