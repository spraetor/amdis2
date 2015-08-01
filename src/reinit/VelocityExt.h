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




#ifndef VELOCITYEXT_H
#define VELOCITYEXT_H

#include "AdaptInfo.h"
#include "DOFVector.h"
#include "Element.h"
#include "io/FileWriter.h"
#include "FixVec.h"

namespace reinit {

using namespace AMDiS;

class VelocityExt
{
public:
  VelocityExt(int dim_)
    : dim(dim_),
      nVelDOFs(0),
      lamVec(dim_, dim_ + 1, NO_INIT),
      permutation(dim_)
  {
    indexFace = -1;
  }

  virtual ~VelocityExt()
  {}

  /// Print velDOF to file.
  void printVelDOF(AdaptInfo *adaptInfo, int i = 0)
  {
    FUNCNAME("VelocityExt::printVelDOF()");

    TEST_EXIT(i < static_cast<int>(velDOF.size()))("Illegal index!\n");

    FileWriter *fileWriter = new FileWriter(
			     "VelocityExt->velocity output", 
			     (velDOF[i])->getFeSpace()->getMesh(),
			     const_cast<DOFVector<double> *>(velDOF[i]));
    fileWriter->writeFiles(adaptInfo, false);

    delete fileWriter;
  }

  /// Print origVelDOF to file.
  void printOrigVelDOF(AdaptInfo *adaptInfo, int i = 0)
  {
    FUNCNAME("VelocityExt::printOrigVelDOF()");

    TEST_EXIT(i < static_cast<int>(origVelDOF.size()))("Illegal index!\n");

    FileWriter *fileWriter = new FileWriter(
			     "VelocityExt->interface velocity output", 
			     (origVelDOF[i])->getFeSpace()->getMesh(),
			     const_cast<DOFVector<double> *>(origVelDOF[i]));
    fileWriter->writeFiles(adaptInfo, false);

    delete fileWriter;
  }

  /// Set velocity (one velocity vector).
  void setVelocity(DOFVector<double> *origVelDOF_,
		   DOFVector<double> *velDOF_)
  {
    FUNCNAME("VelocityExt::setVelocity()");

    nVelDOFs = 1;

    TEST_EXIT(origVelDOF_)("illegal velocity vector origVelDOF !\n");
    TEST_EXIT(velDOF_)("illegal velocity vector velDOF !\n");
    TEST_EXIT(origVelDOF_->getFeSpace() == velDOF_->getFeSpace())
      ("different fe spaces !\n");

    origVelDOF.clear();
    velDOF.clear();

    origVelDOF.push_back(origVelDOF_);
    velDOF.push_back(velDOF_);
  }

  /// Set velocity (multiple velocity vectors).
  void setVelocity(std::vector<DOFVector<double> *> &origVelDOF_,
		   std::vector<DOFVector<double> *> &velDOF_)
  {
    FUNCNAME("VelocityExt::setVelocity()");

    nVelDOFs = origVelDOF_.size();

    TEST_EXIT(nVelDOFs > 0)("Illegal number of velocity vectors!\n");
    TEST_EXIT(nVelDOFs == static_cast<int>(velDOF_.size()))("Different sizes!\n");

    for (int i=0; i<nVelDOFs; ++i) {
      TEST_EXIT(origVelDOF_[i])("illegal velocity vector origVelDOF !\n");
      TEST_EXIT(velDOF_[i])("illegal velocity vector velDOF !\n");
      TEST_EXIT((origVelDOF_[i])->getFeSpace() == (velDOF_[i])->getFeSpace())
	("different fe spaces !\n");
    }
      
    origVelDOF = origVelDOF_;
    velDOF = velDOF_;
  };

  /**
   * Calculates the velocity for a vertex on a boundary element.
   * lokInd gives the lokal index of the element vertices.
   * indexV is the index of the vertex(with respect to the numeration on the element),
   * for that the velocity shall be calculated.
   */
  virtual void calcVelocityBoundary(DegreeOfFreedom *lokInd, 
				    const int indexV);

  /**
   * Calculates the velocity for a vertex on a non-boundary element.
   * lokInd gives the lokal index of the element vertices.
   * indexV is the index of the vertex(with respect to the numeration on the element),
   * for that the velocity shall be calculated.
   */
  void calcVelocity(DegreeOfFreedom *lokInd, const int indexV);

  /**
   * Sets barycentric coordinates for 2D for a vertex on a boundary element.
   * l_0. l_1, l_2 are the coordinates.
   * indexV is the index of the vertex(with respect to the numeration on the element),
   * for that the coordinates shall be stored.
   */
  void setBarycentricCoords_2D_boundary(const double &l_0, 
					const double &l_1, 
					const double &l_2, 
					const int indexV);

  /**
   * Sets barycentric coordinates for 2D for a vertex on a non-boundary element.
   * l_0. l_1, l_2 are the coordinates.
   */
  void setBarycentricCoords_2D(const double &l_0, 
			       const double &l_1, 
			       const double &l_2);

  /**
   * Calculates the barycentric coordinates for 2D for a vertex on a 
   * non-boundary element.
   */
  void calcBarycentricCoords_2D(const double &c_delta, 
				const double &c_alpha, 
				const double &norm_zhminusyh, 
				const double &norm_xhminusyh);

  /**
   * Sets barycentric coordinates for 3D for a vertex on a boundary element
   * (in case the coordinates are already known).
   * l_0, l_1, l-2, l_3 are the coordinates.
   * indexV is the index of the vertex(with respect to the numeration on the element),
   * for that the coordinates shall be stored.
   */
  void setBarycentricCoords_3D_boundary(const double &l_0, 
					const double &l_1, 
					const double &l_2, 
					const double &l_3, 
					const int indexV);

  /// Calculates the barycentric coordinates for 3D for a boundary vertex.
  void calcBarycentricCoords_3D_boundary(const DimVec<double> sp1, 
					 const DimVec<double> sp2, 
					 const double lambda, 
					 int i);

  /**
   * Expands 3 coordinates from a face update to 4 coordinates and stores them in 
   * lamVec[index]. vertNum is the index of the element face the coordinates 
   * belong to (0, 1 or 2).
   */
  void copyAndExpandFaceCoords_3D(int vertNum, int index);

  /**
   * Sets barycentric coordinates for 3D for a vertex on a non-boundary element.
   * l_0, l_1, l-2, l_3 are the coordinates.
   */
  void setBarycentricCoords_3D(const double &l_0, 
			       const double &l_1, 
			       const double &l_2, 
			       const double &l_3);

  /**
   * Stores the index of element face with the shortest distance to the interface so far.
   * (with respect to the order of calculation)
   */
  void setIndexFaceWithShortestDist(int indexC);

  /**
   * Sets the permutation of the vertices.
   * vertNum is the vertex number(with respect to the numeration on the element),
   * for that the update is calculated.
   * mTrav is an index for determining the way of mesh-traverse.
   * mTrav=0:level
   * mTrav=1:traverse
   */
  void setPermutation(int vertNum, int mTrav);

  /// Sets the permutation.
  void setPermutation(int vertNum);

  /// Sets the permutation of the vertices in 2D.
  void setPermutation_2D(int i_0, int i_1, int i_2);

  /// Sets the permutation of the vertices in 3D.
  void setPermutation_3D(int i_0, int i_1, int i_2, int i_3);

  /// Swaps two vertices in the permutation.
  void swapVertices(int i1, int i2);

 protected:
  /// Original velocity vector.
  std::vector<DOFVector<double> *> origVelDOF;

  /// Dimension of mesh.
  int dim;

  /// DOF vector with extended velocity.
  std::vector<DOFVector<double> *> velDOF;

  /// Number of velocity vectors to be extended.
  int nVelDOFs;

  /**
   * Vector with barycentric coordinates of points from which the velocity is extended.
   * Whole vector is needed for the calculation on the boundary elements.
   * For non-boundary elements the first place is used.
   * If the update is calculated with the element faces, the places two-dim are used.
   */
  VectorOfFixVecs<DimVec<double> > lamVec;

  /**
   * If the update is calculated with the element faces, indexFace stores the number of the element
   * face with the shortest distance to the interface.
   * Otherwise indexFace is -1;
   */
  int indexFace;

  /// List with the permutation of the vertices.
  DimVec<int> permutation;
};

} 

using reinit::VelocityExt;

#endif  // VELOCITYEXT_H
