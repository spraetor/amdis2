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



/** \file Line.h */

#ifndef AMDIS_LINE_H
#define AMDIS_LINE_H

#include "Element.h"
#include "BoundaryObject.h"

namespace AMDiS {

  /** \ingroup Triangulation 
   * \brief
   * A Line is an 1-dimensional Element.
   *
   * A Line and its refinements:
   *
   * <img src = "line.png">
   */
  class Line : public Element
  {
  public:
    /// calls base class contructor.
    Line(Mesh* aMesh) 
      : Element(aMesh) 
    {}

    /// implements Element::getVertexOfEdge
    inline int getVertexOfEdge(int i, int j) const 
    {
      return vertexOfEdge[i][j];
    }

    /// implements Element::getVertexOfPosition
    virtual int getVertexOfPosition(GeoIndex position,
				    int positionIndex,
				    int vertexIndex) const;

    virtual int getPositionOfVertex(int side, int vertex) const 
    {
      static int positionOfVertex[2][2] = {{0, -1}, {-1, 0}};
      return positionOfVertex[side][vertex];
    }

    /// implements Element::getGeo
    inline int getGeo(GeoIndex i) const 
    {
      switch (i) {
      case VERTEX: case PARTS: case NEIGH:
	return 2;
	break;
      case EDGE: case FACE:
	return 0;
	break;
      case CENTER: case DIMEN:
	return 1;
	break;
      case PROJECTION: case BOUNDARY:
	return 2;
	break;
      default:
	ERROR_EXIT("invalid geo-index\n");
	return 0;
      }
    }

    inline int getEdgeOfFace(int /* face */, int /*edge*/ ) const 
    {
      ERROR_EXIT("This does not work in 1D!\n");
      return 0;
    }

    DofEdge getEdge(int localEdgeIndex) const
    {
      ERROR_EXIT("This does not work in 1D!\n");
      return DofEdge();
    }

    DofFace getFace(int localFaceIndex) const
    {
      ERROR_EXIT("This does not work in 1D!\n");
      return DofFace();
    }

    /// implements Element::sortFaceIndices
    void sortFaceIndices(int face, FixVec<int,WORLD> &vec) const;
  

    /// implements Element::clone
    inline Element *clone() 
    { 
      return new Line(mesh); 
    }

    /// implements Element::getSideOfChild()
    int getSideOfChild(int child, int side, int) const 
    {
      FUNCNAME_DBG("Line::getSideOfChild()");
      TEST_EXIT_DBG(child == 0 || child == 1)("Child must be either 0 or 1!\n");
      TEST_EXIT_DBG(side >= 0 && side <= 1)("Side must be either 0 or 1!\n");
      return sideOfChild[child][side];
    }

    int getSubObjOfChild(int childnr, GeoIndex subObj, int ithObj, int elType) const
    {
      ERROR_EXIT("Not yet implemented!\n");
      return 0;
    }

    /// implements Element::getVertexOfParent()
    int getVertexOfParent(int child, int side, int) const 
    {
      FUNCNAME_DBG("Line::getVertexOfParent()");
      TEST_EXIT_DBG(child == 0 || child == 1)("Child must be either 0 or 1!\n");
      TEST_EXIT_DBG(side >= 0 && side <= 1)("Side must be between 0 and 2!\n");
      return vertexOfParent[child][side];
    }


    /// implements Element::hasSide
    inline bool hasSide(Element* /*sideElem*/) const 
    { 
      ERROR_EXIT("a Line has no side elements!\n");
      return false; 
    }

    /// Returns true because this element is a Line.
    inline bool isLine() const 
    { 
      return true; 
    }

    /// Returns false because this element is a Line.
    inline bool isTriangle() const 
    { 
      return false; 
    }

    /// Returns false because this element is a Line
    inline bool isTetrahedron() const 
    { 
      return false; 
    }

    /// Element type number is not used in 1d, so return 0.
    inline int getChildType(int) const
    {
      return 0;
    }
  
    std::string getTypeName() const 
    { 
      return "Line"; 
    }

    void getNodeDofs(const FiniteElemSpace*, BoundaryObject, 
		     DofContainer&, bool) const
    {
      FUNCNAME("Line::getNodeDofs()");
      ERROR_EXIT("Not yet implemented!\n");
    }

    void getHigherOrderDofs(const FiniteElemSpace*, BoundaryObject, 
			    DofContainer&, bool, std::vector<GeoIndex>*) const
    {
      FUNCNAME("Line::getHigherOrderDofs()");
      ERROR_EXIT("Not yet implemented!\n");
    }


    void getSubBoundary(BoundaryObject bound, 
			std::vector<BoundaryObject> &subBound) const
    {
      FUNCNAME("Line::getSubBoundary()");
      ERROR_EXIT("Not yet implemented!\n");
    }

  protected:
    /// vertexOfEdge[i][j] is the local number of the j-th vertex of the i-th 
    /// edge of this element.
    static const int vertexOfEdge[1][2];

    static const int sideOfChild[2][2];

    static const int vertexOfParent[2][2];
  };

}

#endif // AMDIS_LINE_H
