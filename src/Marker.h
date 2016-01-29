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
    /**
     * Parameters read by Marker, with name 'MARKER'
     *   MARKER->p:                    Power in estimator norm, \ref p
     *   MARKER->info:                 Info level, \ref info
     *   MARKER->max refinement level: \ref maxRefineLevel
     *   MARKER->min refinement level: \ref minRefineLevel
     **/
    Marker(std::string name, int row);

    /// destructor
    virtual ~Marker() {}

    /// Marks element with newMark. If \ref maximumMarking is set, the element
    /// is marked only if the new mark is bigger than the old one. The return
    /// value specifies whether the element has been marked, or not.
    void setMark(Element* el, char newMark);

    /// Can be used by sub classes. Called before traversal.
    virtual void initMarking(AdaptInfo& adaptInfo, Mesh* mesh);

    /// Can be used by sub classes. Called after traversal.
    virtual void finishMarking(AdaptInfo& adaptInfo);

    /// Marks one element.
    virtual void markElement(AdaptInfo& adaptInfo, ElInfo* elInfo);

    /// Marking of the mesh.
    virtual Flag markMesh(AdaptInfo& adaptInfo, Mesh* mesh);

    /// Sets \ref maximumMarking.
    void setMaximumMarking(bool maxMark)
    {
      maximumMarking = maxMark;
    }

    int getElMarkRefine() const
    {
      return elMarkRefine;
    }

    int getElMarkCoarsen() const
    {
      return elMarkCoarsen;
    }

    /// Returns \ref name of the Marker
    std::string getName() const
    {
      return name;
    }

    /// Creates a scalr marker depending on the strategy set in parameters.
    static Marker* createMarker(std::string name, int row_);

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
    /**
     * No extra parameters read by GRMarker.
     **/
    GRMarker(std::string name, int row);

    /// Implementation of Marker::markElement().
    virtual void markElement(AdaptInfo& adaptInfo, ElInfo* elInfo) override;
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
    /**
     * Parameters read by MSMarker, with name 'MARKER'
     *   MARKER->MSGamma:   \ref MSGamma
     *   MARKER->MSGammaC:  \ref MSGammaC
     **/
    MSMarker(std::string name, int row);

    /// Implementation of Marker::initMarking().
    virtual void initMarking(AdaptInfo& adaptInfo, Mesh* mesh) override;

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
    /**
     * Parameters read by ESMarker, with name 'MARKER'
     *   MARKER->ESTheta:   \ref ESTheta
     *   MARKER->ESThetaC:  \ref ESThetaC
     **/
    ESMarker(std::string name, int row);

    /// Implementation of Marker::initMarking().
    virtual void initMarking(AdaptInfo& adaptInfo, Mesh* mesh) override;

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
    /**
     * Parameters read by GERSMarker, with name 'MARKER'
     *   MARKER->GERSThetaStar:   \ref GERSThetaStar
     *   MARKER->GERSNu:          \ref GERSNu
     *   MARKER->GERSThetaC:      \ref GERSThetaC
     **/
    GERSMarker(std::string name, int row);

    /// Implementation of Marker::markMesh().
    virtual Flag markMesh(AdaptInfo& adaptInfo, Mesh* mesh);

  protected:
    /// Refinement marking function.
    void markElementForRefinement(AdaptInfo& adaptInfo, ElInfo* elInfo);

    /// Coarsening marking function.
    void markElementForCoarsening(AdaptInfo& adaptInfo, ElInfo* elInfo);

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
