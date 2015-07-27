#include "SurfaceAssembler.h"

#include "Mesh.h"
#include "FiniteElemSpace.h"
#include "Operator.h"
#include "SubQuadrature.h"
#include "QPInfo.h"

namespace AMDiS 
{
  /// Creates a SurfaceAssembler conforming to operate for the given \ref coords.
  SurfaceAssembler::SurfaceAssembler(Operator *operat,
                        				     const FiniteElemSpace *rowFeSpace,
                        				     const FiniteElemSpace *colFeSpace,
                        				     VectorOfFixVecs<DimVec<double> > &coords) 
    : Assembler(operat, rowFeSpace, colFeSpace, -1), 
      coords_(coords)
  {
    TEST_EXIT(rowDim_ == colDim_)("rowDim_ != colDim_\n");

    SubQuadrature *subQuadrature;

    if (rowQuad11_) {
      subQuadrature = new SubQuadrature(rowQuad11_, rowDim_);
      rowQuad11_ = colQuad11_ = subQuadrature;
      subQuadrature->scaleQuadrature(coords_);
      rowQPInfo11_ = colQPInfo11_ = QPInfo::provideQPInfo(rowQuad11_, NULL);
    }
    if (rowQuad10_) {
      subQuadrature = new SubQuadrature(rowQuad10_, rowDim_);
      rowQuad10_ = colQuad10_ = subQuadrature;
      subQuadrature->scaleQuadrature(coords_);
      rowQPInfo10_ = colQPInfo10_ = QPInfo::provideQPInfo(rowQuad10_, NULL);
    }
    if (rowQuad01_) {
      subQuadrature = new SubQuadrature(rowQuad01_, rowDim_);
      rowQuad01_ = colQuad01_ = subQuadrature;
      subQuadrature->scaleQuadrature(coords_);
      rowQPInfo01_ = colQPInfo01_ = QPInfo::provideQPInfo(rowQuad01_, NULL);
    }
    if (rowQuad00_) {
      subQuadrature = new SubQuadrature(rowQuad00_, rowDim_);
      rowQuad00_ = colQuad00_ = subQuadrature;
      subQuadrature->scaleQuadrature(coords_);
      rowQPInfo00_ = colQPInfo00_ = QPInfo::provideQPInfo(rowQuad00_, NULL);
    }
  }

  /// Destructor
  SurfaceAssembler::~SurfaceAssembler()
  {
    if (rowQuad11_) 
      delete rowQuad11_;
    if (rowQuad10_) 
      delete rowQuad10_;
    if (rowQuad01_) 
      delete rowQuad01_;
    if (rowQuad00_) 
      delete rowQuad00_;
  }

  
  void SurfaceAssembler::adaptSurfaceAssembler(VectorOfFixVecs<DimVec<double> > &coords)
  {
    coords_ = coords;
    if (rowQuad11_) 
      dynamic_cast<SubQuadrature*>(rowQuad11_)->scaleQuadrature(coords_);
    if (rowQuad10_) 
      dynamic_cast<SubQuadrature*>(rowQuad10_)->scaleQuadrature(coords_);
    if (rowQuad01_) 
      dynamic_cast<SubQuadrature*>(rowQuad01_)->scaleQuadrature(coords_);
    if (rowQuad00_) 
      dynamic_cast<SubQuadrature*>(rowQuad00_)->scaleQuadrature(coords_);
  }

  
  bool SurfaceAssembler::initElementVector(const ElInfo *elInfo)
  {
    if (Assembler::initElementVector(elInfo)) {
      FixVec<WorldVector<double>, VERTEX> worldCoords(rowDim_ - 1, NO_INIT);

      // transform barycentric coords to world coords
      for (int i = 0; i < rowDim_; i++)
        elInfo->coordToWorld(coords_[i], &worldCoords[i]);

      // set determinant for world coords of the side
      det_ = elInfo->calcDet(worldCoords);

      return true;
    } else {
      return false;
    }
  }
  
} // end namespace AMDiS
