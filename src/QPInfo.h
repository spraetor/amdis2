/** \file QPInfo.h */

#pragma once

#include <map>

#include "AMDiS_fwd.h"
#include "FixVec.h"

namespace AMDiS 
{
  /// Stores informations at quadrature points of the current element.
  class QPInfo 
  {
  public:
    /// Sets \ref currentElInfo_ to elInfo and all valid flags to false.
    void initElement(const ElInfo *elInfo);

    /// Returns coordinates at quadrature points.
    WorldVector<double> *getCoordsAtQPs(int numPoints);

    /** \brief
     * Returns vector values at quadrature points. If quadFast is set it will be 
     * used for a more efficient evaluation.
     */
    DenseVector<double>& getVecAtQPs(const DOFVector<double>*, 
                        					   int numPoints,
                        					   const FastQuadrature *quadFast = NULL);

    /** \brief
     * Returns gradient values at quadrature points. If quadFast is set it will be 
     * used for a more efficient evaluation.
     */
    DenseVector<WorldVector<double> >& getGrdAtQPs(const DOFVector<double>*,
                                      				     int numPoints,
                                      				     const FastQuadrature *quadFast = NULL);

    /** \brief
     * Returns D2 values at quadrature points. If quadFast is set it will be 
     * used for a more efficient evaluation.
     */
    DenseVector<WorldMatrix<double> >& getD2AtQPs(const DOFVector<double>*,
                                      				    int numPoints,
                                      				    const FastQuadrature *quadFast = NULL);


    /// Returns element normals at quadrature points.
    WorldVector<double> **getElementNormalAtQPs(int numPoints);

    ///
    DimVec<WorldVector<double> > **getGrdLambdaAtQPs(int numPoints);

    /// Returns a QPInfo instance for the given quadrature.
    static QPInfo *provideQPInfo(const Quadrature*, const FastQuadrature*);

    /// Deletes the QPInfo instance for the given quadrature.
    static void clearQPInfo(const Quadrature*, const FastQuadrature*);

    /// Deletes all QPInfo instances.
    static void clearAllQPInfos();

  protected:
    /// Constructor. Called by \ref provideQPInfo().
    QPInfo(const Quadrature*);

    /// Destructor. Called by \ref clearQPInfo() and \ref clearAllQPInfos().
    ~QPInfo();

  protected:
    /// Structure which stores infos about one DOFVector.
    class VecQPInfo 
    {
    public:
      /// Constructor.
      VecQPInfo() 
      	: valNumPointsValid_(0),
      	  grdNumPointsValid_(0),
      	  D2NumPointsValid_(0)
      {}

      /// Values at quadrature points.
      DenseVector<double> valAtQPs;

      /// Gradients at quadrature points
      DenseVector<WorldVector<double> > grdAtQPs_;
    
      /// D2 at quadrature points
      DenseVector<WorldMatrix<double> > D2AtQPs_;

      /// valid flag for values
      int valNumPointsValid_;

      /// valid flag for gradients
      int grdNumPointsValid_;

      /// valid flag for D2
      bool D2NumPointsValid_;
    };

    /// Quadrature of this QPInfo
    const Quadrature *quadrature_;

    /// Set to \ref quadrature_->getNumPoints().
    int numPoints_;

    /// ElInfo of the current element
    const ElInfo *currentElInfo_;

    /// Coords at quadrature points
    WorldVector<double> *coordsAtQPs_;

    /// Valid flag for coords
    int coordsNumPointsValid_;

    /// Map of all vector infos
    std::map<const DOFVector<double>*, VecQPInfo*> vecQPInfos_;

    /// element normal at quadrature points (array of pointers)
    WorldVector<double> **elementNormalAtQPs_;

    /// for constant values at all QPs (all entries point to same memory)
    WorldVector<double> **elementNormalConst_;
  
    /// valid flag for element normals
    int elementNormalNumPointsValid_;

    /// gradient of barycentric coordinates at QPs (array of pointers)
    DimVec<WorldVector<double> > **grdLambdaAtQPs_;

    /// for constant values at all QPs (all entries point to same memory)
    DimVec<WorldVector<double> > **grdLambdaConst_;

    /// number of valid points of grdLambdaAtQPs_
    int grdLambdaNumPointsValid_;

    /// Static map of all QPInfos. Used by \ref provideQPInfo().
    static std::map<const Quadrature*, QPInfo*> qpInfos_;
  };

} // end namespace AMDiS
