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




#ifndef AMDIS_ELEMENTLEVELSET_H
#define AMDIS_ELEMENTLEVELSET_H

#include "AMDiS_fwd.h"
#include "ElementFunction.h"
#include "FixVec.h"
#include "Initfile.h"

namespace reinit {
  
  using namespace AMDiS;

// ===========================================================================
// ===== class ElementLevelSet ===============================================
// ===========================================================================
//
// Class Description:
// The class ElementLevelSet contains the functionality for
//    - calculating the intersection points resulting from an intersection
//      of the boundary (level set zero) with an element,
//    - calculating the status of an element, i.e. is the element cut
//      by the zero level set or not.
// The calculation of the intersection points is done as follows:
// The level set function is linearly approximated on the element, i.e. its 
// graph is approximated by the plane through the level set values of the 
// element vertices. We approximate the intersection points by the 
// intersection of the plane with the element edges.
// If in 3-dimensional finite element space the intersection produced 4
// intersection points S0, S1, S2, S3, the intersection points in  
// elIntersecPoints are rearranged so that S1 and S2 divides the intersection
// plane in two (dim - 1)-dimensional simplices.
//
// Constants indicating the level set status of element:
//    LEVEL_SET_INTERIOR   -  element is in domain where level set function
//                            is negative
//    LEVEL_SET_BOUNDARY   -  element is in domain where level set function
//                            is positive
//    LEVEL_SET_EXTERIOR   -  element is cut by the zero level set
//
// Main routines:
// setElementLevelSet()          - Defines the level set function for the
//                                 following calculations.
// createElementLevelSet()       - Calculates level set status of element and 
//                                 intersection points if needed.
// calculateElementLevelSetVal() - Calculates values of the level set 
//                                 function in the element vertices.
// setElement()                  - Sets elInfo.
// getElementLevelSetStatus()    - Returns status of element.
// getElementIntersecPoints()    - Returns intersection points.
// getElVertStatusVec()          - Returns vector with status information 
//                                 for each vertex.
// ===========================================================================
class ElementLevelSet
{
 public:
  ElementLevelSet(const char *name_,
		  ElementFunction<double> *lSFct_,
		  Mesh *mesh_) 
    : name(name_),
      elInfo(NULL),
      lastEl(NULL),
      level_set_domain(LEVEL_SET_UNDEFINED),
      numIntersecPoints(0),
      elStatus(LEVEL_SET_UNDEFINED),
      numElVertexInterior(0),
      numElVertexBoundary(0),
      numElVertexExterior(0),
      LS_VAL_TOL(1.e-8),
      LS_VAL_MIN(1.e-8),
      SP_BARY_TOL(1.e-7)
  {
    FUNCNAME("ElementLevelSet::ElementLevelSet()");

    TEST_EXIT(lSFct_ || mesh_)
      ("illegal initialization of ElementLevelSet!\n");

    lSFct = lSFct_;
    mesh = mesh_;
    dim = mesh->getDim();

    elIntersecPoints = new VectorOfFixVecs<DimVec<double> >(dim, 
							    MAX_INTERSECTION_POINTS, 
							    NO_INIT);
    elVertexStatusVec = new int[dim + 1];
    elVertexLevelSetVec = new double[dim + 1];

    int setElementLevelSetTol = 0;
    Parameters::get(name + "->set ElementLevelSet tolerances", 
		    setElementLevelSetTol);
    if (setElementLevelSetTol) {

      Parameters::get(name + "->LS_VAL_TOL", LS_VAL_TOL);
      Parameters::get(name + "->LS_VAL_MIN", LS_VAL_MIN);
      Parameters::get(name + "->SP_BARY_TOL", SP_BARY_TOL);

      TEST_EXIT(LS_VAL_TOL > 0)("illegal LS_VAL_TOL\n");
      TEST_EXIT(LS_VAL_MIN > 0)("illegal LS_VAL_MIN\n");
      TEST_EXIT(SP_BARY_TOL > 0)("illegal SP_BARY_TOL\n");
    }
  }

  ~ElementLevelSet()
  {
    if (elVertexStatusVec)
      delete [] elVertexStatusVec;
    if(elVertexLevelSetVec)
      delete [] elVertexLevelSetVec;
    if (elIntersecPoints)
      delete elIntersecPoints;
  }

  /**
   * Calculates LevelSet-status of element and its intersection points 
   * with the zero level set if necessary.
   *
   * Result:
   *   LEVEL_SET_BOUNDARY: Element elInfo is intersected by levelSetFct.
   *   LEVEL_SET_EXTERIOR / LEVEL_SET_INTERIOR: Element lies completely on 
   *     one side of the zero level set. 
   *
   * After proceeding this function, information about the level set 
   * status is given in:
   *   elStatus: status of element (LEVEL_SET_BOUNDARY, LEVEL_SET_INTERIOR or
   *             EXTERIOR)
   *   elVertexStatusVec: stores status of each vertex of element
   *   elVertexLevelSetVec: stores level set function value of each vertex of
   *                        element
   *   numElVertexInterior: number of vertices of element with status
   *                        LEVEL_SET_INTERIOR
   *   numElVertexExterior: number of vertices of element with status
   *                        LEVEL_SET_EXTERIOR
   *   numElVertexBoundary: number of vertices of element with status
   *                        LEVEL_SET_BOUNDARY
   *   elIntersecPoints: stores the intersection points produced by the 
   *                     intersection of element with the zero level set
   *   numIntersecPoints: number of intersection points
   */
  int createElementLevelSet(const ElInfo *elInfo_, 
			    const bool doCalcIntersecPts_ = true);

  /// Gets value of level set function at point given in barycentric coordinates.
  inline double calcLevelSetFct(const DimVec<double>& bary) 
  {
    return (*lSFct)(bary);
  }

  /// Resets level set information on element.
  inline void resetElement() 
  {
    numElVertexInterior = 0;
    numElVertexBoundary = 0;
    numElVertexExterior = 0;
    numIntersecPoints = 0;
    elStatus = LEVEL_SET_UNDEFINED;
  }

  /// Defines current element (elInfo).
  inline void setElement(const ElInfo *elInfo_) 
  {
    elInfo = elInfo_;
    resetElement();
  }

  /// Set level_set_domain.
  inline void setLevelSetDomain(int status_) 
  {
    TEST_EXIT(status_ == LEVEL_SET_INTERIOR ||
	      status_ == LEVEL_SET_EXTERIOR ||
	      status_ == LEVEL_SET_BOUNDARY)("illegal level set status !\n");
    level_set_domain = status_;
  }

  /// Functions to set tolerances for intersection point calculation.
  inline void setLsValTol(double tol) 
  {
    LS_VAL_TOL = tol;
  }

  inline void setLsValMin(double min) 
  {
    LS_VAL_MIN = min;
  }

  inline void setSpBaryTol(double tol) 
  {
    SP_BARY_TOL = tol;
  }

  /// Get level_set_domain.
  inline const int& getLevelSetDomain() const 
  {
    return level_set_domain;
  }

  /// Get LevelSet-Status of element.
  inline const int& getElementLevelSetStatus() const 
  {
    return elStatus;
  }

  /// Get number of vertices which are intersection points.
  inline const int& getNumVertIntPoints() const 
  {
    FUNCNAME("ElementLevelSet::getNumVertIntPoints");
    TEST_EXIT(numElVertexBoundary == 0)("numElVertexBoundary should be zero!\n");
    return numElVertexBoundary;
  };

  /// Get vector elVertexStatusVec.
  inline const int *getElVertStatusVec() const 
  {
    return elVertexStatusVec;
  }

  /// Get i-th component of vector elVertexLevelSetVec.
  inline double getElVertLevelSetVec(const int i) const 
  {
    return elVertexLevelSetVec[i];
  }

  /// Get vector elVertexLevelSetVec.
  inline const double *getElVertLevelSetVec() const 
  {
    return elVertexLevelSetVec;
  }

  /// Get levelSetFct.
  inline ElementFunction<double> *getLevelSetFct() const 
  {
    return lSFct;
  }

  /// Get mesh.
  inline Mesh *getMesh() const 
  {
    return mesh;
  }

  /// Get dim.
  inline int getDim() const 
  {
    return dim;
  }

  /// Get the intersection points.
  inline VectorOfFixVecs<DimVec<double> > *getElIntersecPoints() const 
  {
    return elIntersecPoints;
  }

  /// Get number of intersection points.
  inline int getNumElIntersecPoints() const 
  {
    return numIntersecPoints;
  }

  /// Calculate exterior normal to intersection plane.
  void calcIntersecNormal(WorldVector<double> &normal);

  /**
   * Gets position of point in element with barycentric coordinates
   * barCoords, i.e. whether point is in the domain with positive
   * (LEVEL_SET_EXTERIOR) or negative (LEVEL_SET_INTERIOR) level set 
   * function values. Uses level set function, thus element vertices
   * may have level set function value zero.
   */
  int getElPos(const DimVec<double> barCoords);

  /**
   * Gets position of element vertex given in barycentric coordinates,
   * i.e. whether element vertex is in the domain with positive
   * (LEVEL_SET_EXTERIOR) or negative (LEVEL_SET_INTERIOR) level set 
   * function values. Uses elVertexLevelSetVec.
   */
  int getVertexPos(const DimVec<double> barCoords);

 protected:
  /// Calculates level set value of each vertex of element.
  void calculateElementLevelSetVal();

  /**
   * Calculates the status of an element.
   *
   * Note: Uses values in elVertexLevelSetVec.
   *
   * Return value:
   *   LEVEL_SET_INTERIOR   element lies completely inside
   *   LEVEL_SET_EXTERIOR   element lies completely outside
   *   LEVEL_SET_BOUNDARY   element is cut by the zero level set
   */
  int calculateElementStatus();

  /**
   * Calculates intersection points of zero level set with element.
   *
   * Note: Uses elVertexLevelSet.
   */
  void calculateIntersecPoints();

  /**
   * Checks whether level set values of element (in elVertexLevelSetVec) 
   * are below a certain bound and corrects them if this is the case.
   *
   * Return value: number of values corrected
   */
  int checkElementLevelSetVal();

  /**
   * Checks whether barycentric coordinate of intersection point is not 
   * too small and corrects it if this is the case.
   *
   * Return value: true  - barycentric coordinate has been corrected
   *               false - barycentric coordinate is ok 
   */
  bool checkIntersecBary(double &bary);
  
  /**
   * Sort intersection points S0, S1, S2 and S3 in \ref elIntersecPoints in 
   * such a way that afterwards, a line through S1 and S2 divides the 
   * intersection plane into two (\ref dim - 1)-dimensional simplices.
   */
  void sortIntersecPoints_4IP3D();
  
  /**
   * Calculate exterior normal to intersection plane for dimension 2.
   */
  void calcIntersecNormal_2d(WorldVector<double> &normal);

  /**
   * Calculate exterior normal to intersection plane for dimension 3.
   */
  void calcIntersecNormal_3d(WorldVector<double> &normal);

 public:
  /**
   * Constants characterizing element position with respect to zero level set.
   */
  static const int LEVEL_SET_INTERIOR = -1;
  static const int LEVEL_SET_BOUNDARY = 0;
  static const int LEVEL_SET_EXTERIOR = 1;
  static const int LEVEL_SET_UNDEFINED = -2;

 protected:
  /** 
   * Name of this object.
   */
  std::string name;

  /**
   * Level set function. 
   */
  ElementFunction<double> *lSFct;

  /**
   * Mesh.
   */
  Mesh *mesh;

  /**
   * elInfo of element.
   */
  const ElInfo *elInfo;

  /**
   * Pointer to last element processed calculations on whithin this class.
   */
  Element *lastEl;

  /**
   * Indicator which can be used for example for function evaluation 
   * or integration on subelements. Indicates whether point/subelement ... 
   * is inside (LEVEL_SET_INTERIOR) or outside (LEVEL_SET_EXTERIOR) the 
   * zero level set or is cut by the zero level set (LEVEL_SET_BOUNDARY).
   */
  int level_set_domain;

  /**
   * Dimension of the problem. dim + 1 is the number of vertices 
   * of element.
   */
  int dim;

  /**
   * Vector for intersection points produced by the intersection of linearly 
   * approximated level set function with the edges of element.
   */
  VectorOfFixVecs<DimVec<double> > *elIntersecPoints;

  /**
   * Number of intersection points.
   */
  int numIntersecPoints;

  /**
   * LevelSet-Status of element.
   */
  int elStatus;

  /**
   * Holds for each vertex of element the information about the position 
   * of the vertex with respect to the zero level set.
   */
  int *elVertexStatusVec;

  /**
   * Stores for each vertex of element the level set of the vertex.
   */
  double *elVertexLevelSetVec;

  /**
   * Number of vertices in element with level set status LEVEL_SET_INTERIOR.
   */
  int numElVertexInterior;

  /**
   * Number of vertices in element with level set status LEVEL_SET_BOUNDARY.
   *
   * Note: should be zero
   */
  int numElVertexBoundary;

  /**
   * Number of vertices in element with level set status LEVEL_SET_EXTERIOR.
   */
  int numElVertexExterior;

  /**
   * Tolerance used in the calculation of the local minimal level set value.
   * The local minimal level set value depends on the gradient of the 
   * level set function.
   * Used for the calculation of intersection points.
   *
   * If intersection points are too close to vertices, they are slightly
   * moved.
   * IDEA: If d is the distance of an intersection point to vertex v, 
   *       the property
   *
   *           d > LS_VAL_TOL * h
   *
   *       must be satisfied.
   * In the implementation this results in
   *
   *           phi(v) > LS_VAL_TOL * h * grad .
   */
  double LS_VAL_TOL;

  /**
   * Lower bound for level set value on elements cut by the zero level set.
   * Used for the calculation of intersection points.
   */
  double LS_VAL_MIN;

  /**
   * Lower bound for barycentric coordinates of intersection points.
   * 
   * Each component x of the barycentric coordinates of an intersection 
   * point satisfies
   *
   *      SP_BARY_TOL < x < 1 - SP_BARY_TOL .
   */
  double SP_BARY_TOL;

  /*
   * Maximum number of intersection points.
   */
  static const int MAX_INTERSECTION_POINTS = 4;
};

}

using reinit::ElementLevelSet;

#endif  // AMDIS_ELEMENTLEVELSET_H
