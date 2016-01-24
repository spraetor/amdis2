/** \file Line.h */

#pragma once

#include "Element.h"
#include "BoundaryObject.h"

namespace AMDiS
{
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

    /// Implementation of \ref Element::getVertexOfEdge
    virtual int getVertexOfEdge(int i, int j) const override
    {
      return vertexOfEdge[i][j];
    }

    /// Implementation of \ref Element::getVertexOfPosition
    virtual int getVertexOfPosition(GeoIndex position,
                                    int positionIndex,
                                    int vertexIndex) const override;

    /// Implementation of \ref Element::getPositionOfVertex
    virtual int getPositionOfVertex(int side, int vertex) const override
    {
      static constexpr int positionOfVertex[2][2] = {{0, -1}, {-1, 0}};
      return positionOfVertex[side][vertex];
    }

    /// Implementation of \ref Element::getGeo
    virtual int getGeo(GeoIndex i) const override
    {
      switch (i)
      {
      case VERTEX:
      case PARTS:
      case NEIGH:
        return 2;
        break;
      case EDGE:
      case FACE:
        return 0;
        break;
      case CENTER:
      case DIMEN:
        return 1;
        break;
      case PROJECTION:
      case BOUNDARY:
        return 2;
        break;
      default:
        ERROR_EXIT("invalid geo-index\n");
        return 0;
      }
    }

    /// Implementation of \ref Element::getEdgeOfFace (should not be called)
    virtual int getEdgeOfFace(int, int) const override
    {
      FUNCNAME("Line::getEdgeOfFace()");
      ERROR_EXIT("This does not work in 1D!\n");
      return 0;
    }

    /// Implementation of \ref Element::getEdge (should not be called)
    virtual DofEdge getEdge(int) const override
    {
      FUNCNAME("Line::getEdge()");
      ERROR_EXIT("This does not work in 1D!\n");
      return {};
    }

    /// Implementation of \ref Element::getFace (should not be called)
    virtual DofFace getFace(int) const override
    {
      FUNCNAME("Line::getFace()");
      ERROR_EXIT("This does not work in 1D!\n");
      return {};
    }

    /// Implementation of \ref Element::sortFaceIndices
    virtual void sortFaceIndices(int face, FixVec<int,WORLD>& vec) const override;


    /// Implementation of \ref Element::clone
    virtual Element* clone() const override
    {
      return new Line(mesh);
    }

    /// Implementation of \ref Element::getSideOfChild
    virtual int getSideOfChild(int child, int side, int) const override
    {
      FUNCNAME_DBG("Line::getSideOfChild()");
      TEST_EXIT_DBG(child == 0 || child == 1)("Child must be either 0 or 1!\n");
      TEST_EXIT_DBG(side >= 0 && side <= 1)("Side must be either 0 or 1!\n");
      return sideOfChild[child][side];
    }

    /// Implementation of \ref Element::getSubObjOfChild (not implemented)
    virtual int getSubObjOfChild(int, GeoIndex, int, int) const override
    {
      FUNCNAME("Line::getSubObjOfChild()");
      ERROR_EXIT("Not yet implemented!\n");
      return 0;
    }

    /// Implementation of \ref Element::getVertexOfParent
    virtual int getVertexOfParent(int child, int side, int) const override
    {
      FUNCNAME_DBG("Line::getVertexOfParent()");
      TEST_EXIT_DBG(child == 0 || child == 1)("Child must be either 0 or 1!\n");
      TEST_EXIT_DBG(side >= 0 && side <= 1)("Side must be between 0 and 2!\n");
      return vertexOfParent[child][side];
    }

    /// Implementation of \ref Element::isLine
    /// Returns true because this element is a Line.
    virtual bool isLine() const override
    {
      return true;
    }

    /// Implementation of \ref Element::isTriangle
    /// Returns false because this element is a Line.
    virtual bool isTriangle() const override
    {
      return false;
    }

    /// Implementation of \ref Element::isTetrahedron
    /// Returns false because this element is a Line
    virtual bool isTetrahedron() const override
    {
      return false;
    }


    /// Implementation of \ref Element::hasSide (should not be called)
    virtual bool hasSide(Element*) const override
    {
      FUNCNAME("Line::hasSide()");
      ERROR_EXIT("A Line has no side elements!\n");
      return false;
    }

    /// Implementation of \ref Element::getChildType
    /// Element type number is not used in 1d, so return 0.
    virtual int getChildType(int) const override
    {
      return 0;
    }

    /// Implementation of \ref Element::getNodeDofs (not implemented)
    virtual void getNodeDofs(FiniteElemSpace const*, BoundaryObject,
			     DofContainer&, bool) const override
    {
      FUNCNAME("Line::getNodeDofs()");
      ERROR_EXIT("Not yet implemented!\n");
    }

    /// Implementation of \ref Element::getHigherOrderDofs (not implemented)
    virtual void getHigherOrderDofs(FiniteElemSpace const*, BoundaryObject,
		    DofContainer&, bool, std::vector<GeoIndex>*) const override
    {
      FUNCNAME("Line::getHigherOrderDofs()");
      ERROR_EXIT("Not yet implemented!\n");
    }


    /// Implementation of \ref Element::getSubBoundary (snot implemented)
    virtual void getSubBoundary(BoundaryObject,
				std::vector<BoundaryObject>&) const override
    {
      FUNCNAME("Line::getSubBoundary()");
      ERROR_EXIT("Not yet implemented!\n");
    }

  protected:
    /// vertexOfEdge[i][j] is the local number of the j-th vertex of the i-th
    /// edge of this element.
    static constexpr int vertexOfEdge[1][2] = {{0, 1}};

    static constexpr int sideOfChild[2][2] = {{ 0,-1}, {-1, 1}};

    static constexpr int vertexOfParent[2][2] = {{ 0,-1}, {-1, 1}};
  };

} // end namespace AMDiS
