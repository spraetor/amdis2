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


#include "ElInfo.h"
#include "FixVec.h"
#include "SubElInfo.h"
#include "SubPolytope.h"


namespace compositeFEM {
  
  bool SubPolytope::checkIntPoints()
  {
    ////////////////////////////////////////////////////////////////////////////
    //
    //  Do the points in intPoints lie on edges ? And are they inner points,
    //  i.e. they aren't vertices ?
    //  
    //  Return value: true  -  yes
    //                false -  no
    ////////////////////////////////////////////////////////////////////////////

    for (int i = 0; i < numIntPoints; i++) {
      int zeroCounter = 0;
      for (int j = 0; j < dim + 1; j++) {
	if (fabs((*intPoints)[i][j]) <= 1.e-15 ) { 
	  zeroCounter++; 
	}
      }

      /**
       *  Dimension 1
       *
       *  "Inner" points on edges aren't equal to 0.0 in any component.
       */
      if (dim == 1  &&  zeroCounter != 0 ) { 
	return false;
      }

      /**
       *  Dimension 2 
       *
       *  "Inner" points on edges are equal to 0.0 in exactly 1 component.
       */
      if (dim == 2  &&  zeroCounter != 1) {
	return false;
      }

      /**
       *  Dimension 3 
       *
       *  "Inner" points on edges are equal to 0.0 in exactly 2 components.
       */
      if (dim == 3  &&  zeroCounter != 2) { 
	return false;
      }
    }

    return true;
  }

  SubPolytope::SubPolytope(const ElInfo *elInfo_, 
			   VectorOfFixVecs<DimVec<double> > *intPoints_,
			   int numIntPoints_,
			   const int &indexElVertInPol_)
    : elInfo(elInfo_),
      intPoints(intPoints_),
      numIntPoints(numIntPoints_)
  {
    FUNCNAME("SubPolytope::SubPolytope()");

    dim = (*intPoints_)[0].getSize() - 1;

    TEST_EXIT((dim == 1 && numIntPoints == 1) ||
	      (dim == 2 && numIntPoints == 2) ||
	      (dim == 3 && numIntPoints == 3) || 
	      (dim == 3 && numIntPoints == 4))
      ("invalid number of intersection points\n");

    /**
     *  Check whether the points in intPoints are really intersection points,
     *  i.e. lie on edges.
     */
    TEST_EXIT(checkIntPoints() == true)
      ("invalid intersection points - do not lie on edges\n");
    
    /*
     *  Create the subelements the polytope consists of.
     */
    switch ( dim ) {
    case 1:
      createSubElementPolytopeIsSubElement1D(indexElVertInPol_);
      break;
    case 2:
      createSubElementPolytopeIsSubElement2D3D();
      break;
    case 3: 
      if (numIntPoints == 3) {
	createSubElementPolytopeIsSubElement2D3D();
      } else {
	createSubElementsForPolytope3D(indexElVertInPol_);
      }
      break;
    default:
      ERROR_EXIT("invalid dimension\n");
      break;
    }
  }


  void SubPolytope::createSubElementPolytopeIsSubElement1D(int indexElVertInPol)
  {
    //**************************************************************************
    // The intersection of the one-dimensional element (interval) divides 
    // element into two subelements. indexElVertInPol indicates which 
    // subelement to take.
    //**************************************************************************

    FUNCNAME("SubPolytope::createSubElementPolytopeIsSubElement1D()");

    TEST_EXIT(dim == 1 && numIntPoints == 1)("invalid call of this routine\n");

    VectorOfFixVecs<DimVec<double> > *subElVertices = 
      new VectorOfFixVecs<DimVec<double> >(dim, dim + 1, NO_INIT);
    DimVec<double> vertex(dim, 1.0);

    /**
     *  Get the vertex which - with the intersection point in intPoints - forms
     *  a subelement of element.
     */
    if (indexElVertInPol == 0) {
      /**
       *  The vertex in element with barycentric coordinates (1,0) is in 
       *  subelement.
       */
      vertex[1] = 0.0;
    } else {
      /**
       *  The vertex in element with barycentric coordinates (0,1) is in 
       *  subelement.
       */
      vertex[0] = 0.0;
    }

    /**
     *  Collect the vertices of the subelement in subElVertices.
     *
     *  Note: The routines CompositeFEMOperator::getElementMatrix and 
     *        CompositeFEMOperator::getElementVector expect the first vertice 
     *        of a subelement to be a vertice of the corresponding element and
     *        not to be an intersection point.
     */
    (*subElVertices)[0] = vertex;
    (*subElVertices)[dim] = (*intPoints)[0];


    /**
     *  Create a new ElInfo for the subelement.
     */
    subElements.push_back( new SubElInfo(subElVertices, elInfo) );

    TEST_EXIT( subElements.size() == 1 )("error in creating subelements");
    numSubElements = 1;

    delete subElVertices;
  }


  void SubPolytope::createSubElementPolytopeIsSubElement2D3D()
  {
    //**************************************************************************
    // The intersection with element produced dim intersection points.
    // Thus there is exactly one vertice in element which - with the 
    // intersection points - forms a subelement of element.
    // This routine determines this vertex and creates a new object of SubElInfo
    // for the subelement stored in subElements.
    //
    // How does it work ?
    // All intersection points lie on edges of element and aren't vertices. The
    // missing vertex is the intersection point of these edges.  
    // Lying on edges means, that each intersection point has exactly one 
    // component equal to 0.0. The missing vertex is equal to 0.0 in all these
    // components.
    //**************************************************************************

    FUNCNAME("SubPolytope::createSubElementPolytopeIsSubElement2D3D()");

    TEST_EXIT((dim == 2 && numIntPoints == 2) || 
	      (dim == 3 && numIntPoints == 3))
      ("invalid call of this routine\n");

    VectorOfFixVecs<DimVec<double> >*subElVertices = 
      new VectorOfFixVecs<DimVec<double> >(dim, dim + 1, NO_INIT);
    DimVec<double> vertex(dim, 1.0);

    /**
     *  Get the vertex which - with the intersection points intPoints - forms
     *  a subelement of element.
     */
    for (int i = 0; i < numIntPoints; i++) {
      for (int j = 0; j < dim + 1; j++) {
	if ( fabs((*intPoints)[i][j]) <= 1.e-15 ) {
	  vertex[j] = 0.0; 
	};
      }
    }

    /**
     *  Collect the vertices of the subelement in subElVertices.
     *
     *  Note: The routines CompositeFEMOperator::getElementMatrix and 
     *        CompositeFEMOperator::getElementVector expect the first vertice 
     *        of a subelement to be a vertice of the corresponding element and
     *        not to be an intersection point.
     */
    (*subElVertices)[0] = vertex;
    for (int i = 0; i < numIntPoints; i++) {
      (*subElVertices)[i+1] = (*intPoints)[i];
    }


    /**
     *  Create a new ElInfo for the subelement.
     */
    subElements.push_back( new SubElInfo(subElVertices, elInfo) );

    TEST_EXIT( subElements.size() == 1 )("error in creating subelements");
    numSubElements = 1;

    delete subElVertices;
  }


  int SubPolytope::getIndexSecondFaceIntPoint0(int indexFirstFace, int dim)
  {
    for (int i = 0; i < dim + 1; i++) {
      if ( fabs((*intPoints)[0][i]) <= 1.e-15  &&  i != indexFirstFace ) {
	return i;
      }
    }

    ERROR_EXIT("couldn't determine the second face for IntPoint0\n");
    return -1;
  }

 
  void SubPolytope::createSubElementsForPolytope3D(int indexElVertInPol1)
  {
    //**************************************************************************
    // The intersection with element produced four intersection points. Thus the
    // intersection doesn't induce a subelement. This routine divides the
    // subpolytope given by the intersection into three subelements.
    //
    // How does it work ?
    // The intersection points and two vertices of element form a subpolytope
    // of element. First of all, we determine these vertices, and call them
    // A and B. Then we sort the intersection points (S_0, S_1, S_2, S_3)
    // in the following way:
    // S_0 is the first intersection point in the creation-list /ref intPoints_.
    // A and S_0 lie in the face of element opposite to B.
    // S_0 and S_1 are neighbours of A.
    // S_2 is opposite to S_0 in the intersection plane.
    // S_1 and S_3 are neighbours of A.
    // Then the subelements are:
    // A - B - S_0 - S_1 ,  B - S_0 - S_1 - S_2 ,  B - S_0 - S_2 - S_3 .
    //
    // The index of one vertex of element that is in subpolytope is handed to
    // this routine. The subpolytope is determined uniquely by this vertex and
    // the intersection points.
    //**************************************************************************

    FUNCNAME("SubPolytope::createSubElementForPolytope3D");

    TEST_EXIT(dim == 3 && numIntPoints == 4)("invalid call of this routine\n");

    TEST_EXIT(0 <= indexElVertInPol1  &&  indexElVertInPol1 <= 3)
      ("invalid index for vertex of a tetrahedron");

    VectorOfFixVecs<DimVec<double> > *subElVertices = 
      new VectorOfFixVecs<DimVec<double> >(dim, dim + 1, NO_INIT);
    DimVec<double> vertexA(dim, 0.0);
    DimVec<double> vertexB(dim, 0.0);

    int indexElVertInPol2 = 0;  // index of second vertex of element lying in 
    // subpolytope
    // The vertices of element (3D -> tetrahedron) are
    // indexed as usual:
    // 0: (1,0,0,0), 1: (0,1,0,0), 2: (0,0,1,0), 
    // 3: (0,0,0,1)
    bool intPointOnEdge[4][4] = {{false, false, false, false},
				 {false, false, false, false},
				 {false, false, false, false},
				 {false, false, false, false}};
    // /ref intPointOnEdge[i][j] indicates whether there
    // is an intersection point on the edge from  
    // vertice i to vertice j :
    // false : no intersection point
    // true: there is an intersection point
    int indexEdge[2];       // For a vertex lying on an edge indexEdge 
    // stores the indices of the two barycentric 
    // coordinates which are not equal to zero.
    int indexA = 0;            
    int indexB = 0;            
    int indexS_0 = 0;
    int indexS_1 = 0;
    int indexS_2 = 0;
    int indexS_3 = 0;
    int indexSecondFaceIntPoint0 = 0;

    /**
     *  Get the second vertex of element lying in subpolytope.
     *
     *  There is exactly one vertex of element which - with vertex 
     *  indexElVertInPol1 and the intersection points intPoints -
     *  forms a subpolytope of element. It is the vertex adjacent with
     *  indexElVertInPol1 whose common edge with indexElVertInPol1
     *  doesn't contain an intersection point.
     */

    // Get the edges including the intersection points.
    for (int i = 0; i < numIntPoints; i++) {
      int k = 0;
      for (int j = 0; j < dim + 1; j++) {
	if (fabs((*intPoints)[i][j]) > 1.e-15 ) {
	  indexEdge[k] = j;
	  k++;
	}
      }
      intPointOnEdge[indexEdge[0]][indexEdge[1]] = true;
      intPointOnEdge[indexEdge[1]][indexEdge[0]] = true;
    }

    // Get the vertex of element adjacent with indexElVertInPol1 whose
    // common edge with indexElVertInPol1 doesn't contain an
    // intersection point, and store it in indexElVertInPol2.
    for (int i = 0; i < dim + 1; i++) {
      if (intPointOnEdge[indexElVertInPol1][i] == false  &&  
	  i != indexElVertInPol1 ) {
	indexElVertInPol2 = i;
	break;
      }
    }

    /**
     *  Determine A and B, so that intPoint0 is a neighbour of A. 
     *
     *  In the subpolytope A and two intersection points lie in the face
     *  opposite to B. And B and the other two intersection points lie in the
     *  face opposite to A.
     */

    if (fabs((*intPoints)[0][indexElVertInPol1]) <= 1.e-15) {

      // (*intPoints)[0] lies in the face opposite to vertex 
      // /ref indexElVertInPol1.

      indexA = indexElVertInPol2;
      indexB = indexElVertInPol1;
    } else if (fabs((*intPoints)[0][indexElVertInPol2]) <= 1.e-15) {

      // (*intPoints)[0] lies in the face opposite to vertex
      // /ref indexElVertInPol2.

      indexA = indexElVertInPol1;
      indexB = indexElVertInPol2;
    } else {
      ERROR_EXIT("couldn't determine A and B\n");
    }

    /**
     *  Sort the intersection points.
     */

    // (*intPoints)[0] is a neighbour of A (A has been constructed this way).
    indexS_0 = 0;

    if (fabs((*intPoints)[1][indexB]) <= 1.e-15) {

      // (*intPoints)[1] lies in the face opposite to B, thus is a neighbour
      // of A.

      indexS_1 = 1;

      indexSecondFaceIntPoint0 = getIndexSecondFaceIntPoint0(indexB, dim);

      if (fabs((*intPoints)[2][indexSecondFaceIntPoint0]) <= 1.e-15) {

	// (*intPoints)[2] is neighbour of (*intPoints)[0] 

	indexS_2 = 3;
	indexS_3 = 2;
      } else {

	// (*intPoints)[2] is opposite to (*intPoints)[0] in the intersection
	// plane

	indexS_2 = 2;
	indexS_3 = 3;
      }
    } else if (fabs((*intPoints)[1][indexA]) <= 1.e-15) {

      // (*intPoints)[1] lies in the face opposite to A

      indexSecondFaceIntPoint0 = getIndexSecondFaceIntPoint0(indexB, dim);

      if (fabs((*intPoints)[1][indexSecondFaceIntPoint0]) <= 1.e-15) {

	// (*intPoints)[1] is neighbour of (*intPoints)[0], but isn't 
	// neighbour of A

	indexS_3 = 1;

	if (fabs((*intPoints)[2][indexB]) <= 1.e-15) {

	  // (*intPoints)[2] is neighbour of (*intPoints)[0] and neighbour of A

	  indexS_1 = 2;
	  indexS_2 = 3;
	} else {

	  // (*intPoints)[2] is opposite to (*intPoints)[0] in the intersection
	  // plane

	  indexS_1 = 3;
	  indexS_2 = 2;
	}
      } else {

	// (*intPoints)[1] isn't neighbour of (*intPoints)[0], thus lies opposite 
	// to (*intPoints)[0] in the intersection plane

	indexS_2 = 1;

	if (fabs((*intPoints)[2][indexB]) <= 1.e-15) {

	  // (*intPoints)[2] is neighbour of A

	  indexS_1 = 2;
	  indexS_3 = 3;
	} else {

	  // (*intPoints)[2] isn't neighbour of A

	  indexS_1 = 3;
	  indexS_3 = 2;
	}
      }
    } else {
      ERROR_EXIT("IntPoint1 isn't either part of the face opposite to A nor part of the face opposite to B\n");
    }

    /**
     *  For each subelement: Collect the vertices of the subelement in 
     *  subElVertices, create a new SubElInfo for the subelement and
     *  store it in subElements.
     *
     *  Note: The routines CompositeFEMOperator::getElementMatrix and 
     *        CompositeFEMOperator::getElementVector expect the first vertice 
     *        of a subelement to be a vertice of the corresponding element and
     *        not to be an intersection point.
     */

    // Create vertex A and vertex B.
    vertexA[indexA] = 1.0;
    vertexB[indexB] = 1.0;

    // Subelement 1: A - B - S_0 - S_1
    (*subElVertices)[0] = vertexA;
    (*subElVertices)[1] = vertexB;
    (*subElVertices)[2] = (*intPoints)[indexS_0];
    (*subElVertices)[3] = (*intPoints)[indexS_1];
    subElements.push_back( new SubElInfo(subElVertices, elInfo) );

    // Subelement 2: B - S_0 - S_1 - S_2
    (*subElVertices)[0] = vertexB;
    (*subElVertices)[1] = (*intPoints)[indexS_0];
    (*subElVertices)[2] = (*intPoints)[indexS_1];
    (*subElVertices)[3] = (*intPoints)[indexS_2];
    subElements.push_back( new SubElInfo(subElVertices, elInfo) );

    // Subelement 3: B - S_0 - S_2 - S_3
    (*subElVertices)[0] = vertexB;
    (*subElVertices)[1] = (*intPoints)[indexS_0];
    (*subElVertices)[2] = (*intPoints)[indexS_2];
    (*subElVertices)[3] = (*intPoints)[indexS_3];
    subElements.push_back( new SubElInfo(subElVertices, elInfo) );

    TEST_EXIT( subElements.size() == 3 )("error in creating subelements");
    numSubElements = 3;

    delete subElVertices;
  }

}
