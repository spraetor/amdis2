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



/** \file SubPolytope.h */

#ifndef AMDIS_SUBPOLYTOPE_H
#define AMDIS_SUBPOLYTOPE_H

#include <vector>
#include "ElInfo.h"
#include "Mesh.h"
#include "FixVec.h"
#include "SubElInfo.h"


namespace compositeFEM {
  
  using namespace AMDiS;

  // ===========================================================================
  // === class SubPolytope =====================================================
  // ===========================================================================
  //
  // Class description:
  // The class SubPolytope holds the functionality for the division of a
  // subpolytope in subelements and contains a list of these subelements
  // in subElements. The number of subelements is given in numSubElements.
  //
  // The subpolytope of the element elInfo is initially given by the 
  // intersection points (barycentric coordinates with respect to element) 
  // resulting from the intersection of a "plane" with the element.
  // There is/are
  //   1 intersection point in 1-dimensional finite element space,
  //   2 intersection points in 2-dimensional finite element space,
  //   3 or 4 intersection points in 3-dimensional finite element space.
  // A detailed description of how the division into subelements works can be 
  // found in the corresponding routines.
  //
  // Note: This functionality is restricted to the finite element dimensions 
  //       1, 2 and 3.
  //
  // Note: The intersection points are expected to lie exactly on edges. 
  //       That means in the barycentric coordinate vector there are exactly 
  //       as many components equal to zero as are expected for points on edges.
  // 
  // Main routines:
  // SubPolytope() - Creates the subelements for a subpolytope given by the
  //                 intersection points intPoints and stores them in
  //                 subElements.
  // createSubElementPolytopeIsSubElement1D() 
  //               - Creates the subelement in 1-dimensional finite element 
  //                 space.
  // createSubElementPolytopeIsSubElement2D3D()
  //               - Creates the subelement in 2-dimensional finite element
  //                 space and in 3-dimensional finite element space if there
  //                 are 3 intersection points.
  // createSubElementsForPolytope3D()
  //               - Creates the 3 subelements in 3-dimensional finite element
  //                 space if there are 4 intersection points.
  // ============================================================================

  class SubPolytope 
  {
  public:
    /**
     *  Constructor
     *  
     *  Divides the polytope produced by the intersection into subelements.
     *  In dimensions 1 and 3 (case "4 intersection points") indexElVertInPol_
     *  indicates which subpolytope to divide. The element vertice with index 
     *  indexElVertInPol_ is a vertice of the choosen subpolytope.
     *  Here the standard vertice numeration is used (e.g. in dimension 3: 
     *  (1,0,0,0) - index 0, (0,1,0,0) - index 1, ...)
     *  The subelements are stored in subElements.
     */
    SubPolytope(const ElInfo *elInfo_, 
		VectorOfFixVecs<DimVec<double> > *intPoints_,
		int numIntPoints_,
		const int &indexElVertInPol_ = 0);

    /// Destructor
    ~SubPolytope() 
    {
      for (int i = 0; i < numSubElements; i++)
	delete subElements[i];
    }

    /// Returns begin of vector subElements.
    inline std::vector<SubElInfo *>::iterator getSubElementsBegin() 
    {
      return subElements.begin();
    }

    /// Returns end of vector subElements.
    inline std::vector<SubElInfo *>::iterator getSubElementsEnd() 
    {
      return subElements.end();
    }

    /// Returns num-th subelement of polytope.
    inline SubElInfo* getSubElement(int num) 
    {
      FUNCNAME("SubPolytope::getSubElement()");
  
      TEST_EXIT(num <= numSubElements)("invalid index for subelement");
      return subElements[num];
    }

    /// Returns the number of subelements of the polytope.
    inline int getNumSubElements() 
    { 
      return numSubElements; 
    }

  protected:

    /**
     * Checks whether the elements of intPoints are really intersection points, 
     * i.e. lie on edges.
     */
    bool checkIntPoints();

    /**
     *  In 1-dimensional space the intersection of an element always produces 
     *  two subelements. indexElVertInPol indicates which subelement to take.
     *  The resulting subelement is created by this routine.
     */
    void createSubElementPolytopeIsSubElement1D(int indexElVertInPol);

    /**
     *  In 2-dimensional space the intersection of an element always produces 
     *  a subelement. The same in 3-dimensional space, if the intersection has
     *  three intersection points.
     *  The resulting subelement is created by this routine.
     */
    void createSubElementPolytopeIsSubElement2D3D();

    /**
     *  Routine used in createSubElementsForPolytope().
     *
     *  The intersection point /ref intPoints[0] lies on an edge of element
     *  and isn't a vertex of element. Thus it lies in exactly two faces.
     *  In barycentric coordinates, lying in the face opposite to e.g. (0,1,0,0)
     *  means that the second barycentric coordinate is equal to zero. We 
     *  give this face the index 1 (we start indexing by 0).
     *  In this routine we know already one face containing intPoints[0],
     *  namely indexFirstFace. The task of the routine is to get the
     *  second face. 
     */
    int getIndexSecondFaceIntPoint0(int indexFirstFace, int dim);

    /**
     *  If in 3-dimensional space the intersection of an element has four 
     *  intersection points, the resulting polytope is no subelement, but it 
     *  can be divided into three subelements. This is done by this routine.
     */
    void createSubElementsForPolytope3D(int indexElVertInPol);

  protected:
    /// elInfo of the element containing the polytope
    const ElInfo *elInfo;

    /**
     * Intersection points with the element in barycentric coordinates with
     * respect to element
     */
    VectorOfFixVecs<DimVec<double> > *intPoints;

    /// Number of intersection points
    int numIntPoints;

    /// List of the subelements of subpolytope
    std::vector<SubElInfo *> subElements;

    /// Number of subelements
    int numSubElements;

    /// Dimension of the polytope
    int dim;
  };
}

#endif  // AMDIS_SUBPOLYTOPE_H
