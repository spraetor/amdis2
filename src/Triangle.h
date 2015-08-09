/** \file Triangle.h */

#pragma once

#include "Element.h"

namespace AMDiS 
{
  /** \ingroup Triangulation 
   * \brief
   * A Triangle is a 2-dimensional Element.
   *
   * A Triangle and its refinements:
   *
   * <img src = "triangle.png">
   */
  class Triangle : public Element
  {
  public:
    /// calls base class contructor.
    Triangle(Mesh* aMesh) 
      : Element(aMesh) 
    {}

    /// implements Element::clone
    virtual Element *clone() const override
    { 
      return new Triangle(mesh); 
    }

    /// implements Element::getVertexOfEdge
    virtual int getVertexOfEdge(int i, int j) const override
    {
      return vertexOfEdge[i][j];
    }

    /// implements Element::getVertexOfPosition
    virtual int getVertexOfPosition(GeoIndex position,
                        				    int positionIndex,
                        				    int vertexIndex) const override;

    /// implements Element::getGeo
    virtual int getGeo(GeoIndex i) const override
    {
      switch (i) {
      case VERTEX: case PARTS: case NEIGH:
      	return 3;
      	break;
      case EDGE:
      	return 3;
      case FACE:  
      	return 0;
      case CENTER:
      	return 1;
      	break;
      case DIMEN: 
      	return 2;
      	break;
      case BOUNDARY:
      	return 6;
      	break;
      case PROJECTION:
      	return 3;
      	break;
      default:
      	ERROR_EXIT("invalid geo-index\n");
      	return 0;
      }
    }

    /// implements Element::hasSide
    virtual bool hasSide(Element* sideElem) const override;

    /// implements Element::sortFaceIndices
    virtual void sortFaceIndices(int face, FixVec<int, WORLD>& vec) const override;

    /// implements Element::isLine. Returns false because this element is a Triangle
    virtual bool isLine() const override
    { 
      return false; 
    }
 
    /// implements Element::isTriangle. Returns true because this element is a Triangle
    virtual bool isTriangle() const override
    { 
      return true; 
    }

    /// implements Element::isTetrahedron. Returns false because this element is a Triangle
    virtual bool isTetrahedron() const override
    { 
      return false; 
    }

    /// Element type number is not used in 2d, so return 0.
    virtual int getChildType(int) const override
    {
      return 0;
    }

    /// implements Element::getSideOfChild()
    virtual int getSideOfChild(int child, int side, int) const override
    {
      FUNCNAME_DBG("Triangle::getSideOfChild()");
      TEST_EXIT_DBG(child == 0 || child == 1)("child must be in (0,1)\n");
      TEST_EXIT_DBG(side >= 0 && side <= 2)("side must be between 0 and 2\n");

      return sideOfChild[child][side];
    }

    virtual int getSubObjOfChild(int childnr, GeoIndex subObj, 
		                             int ithObj, int elType) const override
    {
      FUNCNAME_DBG("Triangle::getSubObjOfChild()");
      TEST_EXIT_DBG(subObj == EDGE)("Not yet implemented!\n");

      return getSideOfChild(childnr, ithObj, elType);
    }

    /// implements Element::getVertexOfParent()
    virtual int getVertexOfParent(int child, int side, int = 0) const override
    {
      FUNCNAME_DBG("Triangle::getVertexOfParent()");
      TEST_EXIT_DBG(child == 0 || child == 1)("child must be in (0,1)\n");
      TEST_EXIT_DBG(side >= 0 && side <= 2)("side must be between 0 and 2\n");

      return vertexOfParent[child][side];
    }

    virtual int getPositionOfVertex(int side, int vertex) const override
    {
      FUNCNAME_DBG("Triangle::getPositionOfVertex()");
      TEST_EXIT_DBG(side >= 0 && side <= 2)("Wrong side number %d!\n", side);
      TEST_EXIT_DBG(vertex >= 0 && vertex <= 2)("Wrong vertex number %d!\n", vertex);

      static const int positionOfVertex[3][3] = {{-1, 0, 1}, {1, -1, 0}, {0, 1, -1}};
      return positionOfVertex[side][vertex];
    }

    virtual int getEdgeOfFace(int face, int edge) const override
    {
      FUNCNAME_DBG("Triangle::getEdgeOfFace()");
      TEST_EXIT_DBG(face == 0)("face must be zero at triangle\n");
      TEST_EXIT_DBG(edge >= 0 && edge < 3)("invalid edge\n");

      return edge;
    }

    virtual DofEdge getEdge(int localEdgeIndex) const override
    {
      FUNCNAME_DBG("Triangle::getEdge()");
      TEST_EXIT_DBG(localEdgeIndex >= 0 && localEdgeIndex < 3)("invalid edge\n");

      DegreeOfFreedom dof0 = dof[vertexOfEdge[localEdgeIndex][0]][0];
      DegreeOfFreedom dof1 = dof[vertexOfEdge[localEdgeIndex][1]][0];
      return std::minmax({dof0, dof1});
//       return {std::min(dof0, dof1), std::max(dof0, dof1)};
    }

    virtual DofFace getFace(int localFaceIndex) const override
    {
      ERROR_EXIT("This does not work in 2D!\n");
      return {};
    }

    virtual void getNodeDofs(const FiniteElemSpace* feSpace, 
                  			     BoundaryObject bound,
                  			     DofContainer& dofs,
                  			     bool baseDofPtr = false) const override;
    
    virtual void getHigherOrderDofs(const FiniteElemSpace* feSpace, 
                        				    BoundaryObject bound,
                        				    DofContainer& dofs,
                        				    bool baseDofPtr = false,
                        				    std::vector<GeoIndex>* dofGeoIndex = NULL) const override;

    virtual void getSubBoundary(BoundaryObject bound, 
				                        std::vector<BoundaryObject>& subBound) const override;
				

    std::string getTypeName() const
    { 
      return "Triangle"; 
    }
    
    void prepareNextBound(BoundaryObject &bound, int ithChild) const;

  protected:
    /// vertexOfEdge[i][j] is the local number of the j-th vertex of the i-th 
    /// edge of this element.
    static constexpr int vertexOfEdge[3][2] = {{1, 2}, {2, 0}, {0, 1}};

    static constexpr int sideOfChild[2][3] = {{-1, 2, 0}, {2, -1, 1}};

    static constexpr int vertexOfParent[2][3] = {{2, 0, -1}, {1, 2, -1}};
  };

} // end namespace AMDiS
