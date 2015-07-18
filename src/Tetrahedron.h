/** \file Tetrahedron.h */

#pragma once

#include "Element.h"

namespace AMDiS 
{
  /** \ingroup Triangulation 
   * \brief
   * A Tetrahedron is a 3-dimensional Element. 
   *
   * A Tetrahedron and its refinements:
   *
   * <img src = "tetrahedron.png">
   */
  class Tetrahedron : public Element
  {
  public:
    
    /// calls base class contructor.
    Tetrahedron(Mesh* aMesh) 
      : Element(aMesh) 
    {}

    /// implements Element::clone
    virtual Element *clone() const override
    { 
      return new Tetrahedron(mesh); 
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


    virtual int getPositionOfVertex(int side, int vertex) const override
    {
      static const int positionOfVertex[4][4] = {{-1, 0, 1, 2},
					   {0, -1, 1, 2},
					   {0, 1, -1, 2},
					   {0, 1, 2, -1}};
      return positionOfVertex[side][vertex];
    }

    /// implements Element::getGeo
    virtual int getGeo(GeoIndex i) const override
    {
      switch (i) {
      case VERTEX: case PARTS: case NEIGH:
	return 4;
	break;
      case EDGE:
	return 6;
      case FACE:
	return 4;
      case CENTER:
	return 1;
	break;
      case DIMEN:
	return 3;
	break;
      case BOUNDARY:
	return 14;
	break;
      case PROJECTION:
	return 10;
	break;
      default:
	ERROR_EXIT("invalid geo-index\n");
	return 0;
      }
    }

    /// implements Element::hasSide
    virtual bool hasSide(Element* sideElem) const override;

    /// implements Element::sortFaceIndices
    virtual void sortFaceIndices(int face, FixVec<int,WORLD> &vec) const override;

    /// Returns false because this element is a Tetrahedron.
    virtual bool isLine() const override
    { 
      return false; 
    }

    /// Returns false because this element is a Tetrahedron.
    virtual bool isTriangle() const override
    { 
      return false; 
    }
 
    /// Returns true because this element is a Tetrahedron.
    virtual bool isTetrahedron() const override
    { 
      return true; 
    }

    /// Return the new element type of the children.
    virtual int getChildType(int elType) const override
    {
      return (elType + 1) % 3;
    }

    virtual void getNodeDofs(const FiniteElemSpace* feSpace, 
			    BoundaryObject bound,
			    DofContainer& dofs,
			    bool baseDofPtr = false) const override;

    virtual void getNodeDofsAtFace(const FiniteElemSpace* feSpace, 
				  BoundaryObject bound,
				  DofContainer& dofs,
				  bool baseDofPtr) const override;

    virtual void getNodeDofsAtEdge(const FiniteElemSpace* feSpace, 
				  BoundaryObject bound,
				  DofContainer& dofs,
				  bool baseDofPtr) const override;

    virtual void getHigherOrderDofs(const FiniteElemSpace* feSpace,
				    BoundaryObject bound,
				    DofContainer& dofs,
				    bool baseDofPtr = false,
				    std::vector<GeoIndex>* dofGeoIndex = NULL) const override;

    virtual void getSubBoundary(BoundaryObject bound, 
				std::vector<BoundaryObject> &subBound) const override;


    /// implements Element::getSideOfChild()
    virtual int getSideOfChild(int child, int side, int elType = 0) const override
    {
      FUNCNAME_DBG("Tetrahedron::getSideOfChild()");
      TEST_EXIT_DBG(child == 0 || child == 1)("Child must be in (0,1)!\n");
      TEST_EXIT_DBG(side >= 0 && side <= 3)("Side must be between 0 and 3!\n");
      TEST_EXIT_DBG(elType >= 0 && elType <= 2)
	("ElType must be between 0 and 2!\n");

      return sideOfChild[elType][child][side];
    }

    virtual int getSubObjOfChild(int childnr, GeoIndex subObj, 
				 int ithObj, int elType) const override
    {
      FUNCNAME_DBG("Tetrahedron::getSubObjOfChild()");
      TEST_EXIT_DBG(subObj == EDGE || subObj == FACE)("Not yet implemented!\n");

      if (subObj == FACE)
	return getSideOfChild(childnr, ithObj, elType);
      else
	return getEdgeOfChild(childnr, ithObj, elType);
    }

    /// implements Element::getVertexOfParent()
    virtual int getVertexOfParent(int child, int side, int elType = 0) const override
    {
      FUNCNAME_DBG("Tetrahedron::getVertexOfParent()");
      TEST_EXIT_DBG(child == 0 || child == 1)("Child must be in (0,1)!\n");
      TEST_EXIT_DBG(side >= 0 && side <= 3)("Side must be between 0 and 3!\n");
      TEST_EXIT_DBG(elType >= 0 && elType <= 2)
	("ElType must be between 0 and 2!\n");

      return vertexOfParent[elType][child][side];
    }

    virtual int getEdgeOfFace(int face, int edge) const override
    {
      FUNCNAME_DBG("Tetrahedron::getEdgeOfFace()");
      TEST_EXIT_DBG(face >= 0 && face < 4)("Invalid face number!\n");
      TEST_EXIT_DBG(edge >= 0 && edge < 3)("Invalid edge number!\n");

      return edgeOfFace[face][edge];
    }
   
    virtual DofEdge getEdge(int localEdgeIndex) const override
    {
      FUNCNAME_DBG("Tetrahedron::getEdge()");
      TEST_EXIT_DBG(localEdgeIndex >= 0 && localEdgeIndex < 6)("Invalid edge!\n");

      DegreeOfFreedom dof0 = dof[vertexOfEdge[localEdgeIndex][0]][0];
      DegreeOfFreedom dof1 = dof[vertexOfEdge[localEdgeIndex][1]][0];
#if 0
      DofEdge edge = std::make_pair(std::min(dof0, dof1), std::max(dof0, dof1));
#endif
      return std::minmax({dof0, dof1});
    }

    virtual DofFace getFace(int localFaceIndex) const override
    {
      FUNCNAME_DBG("Tetrahedron::getFace()");
      TEST_EXIT_DBG(localFaceIndex >= 0 && localFaceIndex < 4)("Invalid face!\n");

      // Get the three DOFs of the face.
      DegreeOfFreedom dof0 = dof[vertexOfFace[localFaceIndex][0]][0];
      DegreeOfFreedom dof1 = dof[vertexOfFace[localFaceIndex][1]][0];
      DegreeOfFreedom dof2 = dof[vertexOfFace[localFaceIndex][2]][0];
      
      // Sort the three DOFs of the face with respect to their values.
      DegreeOfFreedom dofMin0, dofMin1, dofMin2;
      if (dof0 < dof1 && dof0 < dof2) {
	dofMin0 = dof0;
	dofMin1 = std::min(dof1, dof2);
	dofMin2 = std::max(dof1, dof2);
      } else if (dof1 < dof0 && dof1 < dof2) {
	dofMin0 = dof1;
	dofMin1 = std::min(dof0, dof2);
	dofMin2 = std::max(dof0, dof2);
      } else {
	TEST_EXIT_DBG(dof2 < dof0 && dof2 < dof1)("Should not happen!\n");
	dofMin0 = dof2;
	dofMin1 = std::min(dof0, dof1);
	dofMin2 = std::max(dof0, dof1);
      }
      
      return {dofMin0, dofMin1, dofMin2};
    }

    /// Returns for an edge number its local edge number on a child element. See
    /// \ref edgeOfChild for mor information.
    inline int getEdgeOfChild(int child, int edge, int elType) const
    {
      FUNCNAME_DBG("Tetrahedron::getEdgeOfChild()");
      TEST_EXIT_DBG(child == 0 || child == 1)("Child must be in (0,1)!\n");
      TEST_EXIT_DBG(edge >= 0 && edge <= 5)("Side must be between 0 and 3!\n");
      TEST_EXIT_DBG(elType >= 0 && elType <= 2)
	("ElType must be between 0 and 2!\n");

      return edgeOfChild[elType][child][edge];
    }
    
    void prepareNextBound(BoundaryObject &bound, int ithChild) const;
 
    std::string getTypeName() const 
    { 
      return "Tetrahedron"; 
    }
    
  protected:
    /// nChildEdge[el_type][ichild][dir] gives local index of new edge on 
    /// child[ichild] part of face [2 + dir] on the parent.
    static const unsigned char nChildEdge[3][2][2];

    /** \brief
     * nChildFace[el_type][ichild][dir]                                     
     * gives local index of sub-face on child[ichild] part of face [2+dir] on
     * the parent   
     */
    static const unsigned char nChildFace[3][2][2];

    /** \brief
     * childVertex[el_type][child][i] =
     * parent's local vertex index of new vertex i.
     * 4 stands for the newly generated vertex                            
     */
    static const int childVertex[3][2][4];
 
    /** \brief
     * childEdge[el_type][child][i] =
     * parent's local edge index of new edge i
     * new edge 2 is half of old edge 0,
     * new edges 4,5 are really new edges, and value is different:
     * childEdge[][][4,5] = index of same edge in other child
     *
     * Note: el_type, i.e., the first index of the array, defines the element type
     * of the parent elements.
     */
    static const int childEdge[3][2][6];

    /** \brief
     * adjacentChild[position][ichild]                                      
     * gives number of the adjacent child on a neighbour element:             
     *   position = 0  same position of element and neigh at refinement edge, 
     *   position = 1  different ...                                         
     */
    static const unsigned char adjacentChild[2][2];

    /** \brief
     * childOrientation[el_type][child] =
     *   +1 if orientation is not changed during refinement, 
     *   -1 if orientation is changed during refinement
     */
    static const signed char childOrientation[3][2];

    /// edgeOfDOFs[i][j]: gives the local index of edge with vertices i and j 
    static const unsigned char edgeOfDofs[4][4];

    ///
    static const int edgeOfFace[4][3];

    /// vertexOfEdge[i][j] is the local number of the j-vertex of edge i
    static const int vertexOfEdge[6][2]; 

    /// vertexOfFace[i][j] is the local number of the j-vertex of face i
    static const int vertexOfFace[4][3];

    /// See \ref Element::getSideOfChild for more information.
    static const int sideOfChild[3][2][4];

    /// edgeOfChild[elType][i][j] is the local edge number of the j-th edge within
    /// the i-th children of an element of elType. If the value is -1, the edge is
    /// not included in the element's child. Note that the 0 edge is included in
    /// both children only by its half.
    static const int edgeOfChild[3][2][6];

    /// See \ref Element::getVertexOfParent for more information.
    static const int vertexOfParent[3][2][4];
  };

} // end namespace AMDiS
