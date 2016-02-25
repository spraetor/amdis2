/** \file ElInfo.h */

#pragma once

#include <boost/numeric/mtl/matrix/dense2D.hpp>

#include "Flag.hpp"
#include "Boundary.hpp"
#include "Global.hpp"
#include "FixVec.hpp"
#include "Element.hpp"
#include "AMDiS_fwd.hpp"

namespace AMDiS
{

  /** \ingroup Traverse
   * \brief
   * An ElInfo object holds informations wich are not stored in the corresponding
   * element. It is filled during mesh traversal by the traversal routines.
   * A fill flag determines which informations should be filled and which elements
   * should be visited. Since it is a
   * pure virtual base class for the dimension speciefic ElInfo classes, it must
   * not be instantiated directly.
   * \see ElInfo1d \see ElInfo2d \see ElInfo3d
   */

  class ElInfo
  {
  protected:
    /// Protected constructor. Avoids instatiation of the basis class
    ElInfo();

    /// Protected constructor. Avoids instatiation of the basis class.
    /// \param mesh pointer to the corresponding mesh.
    ElInfo(Mesh* mesh);

  public:
    /// Virtual destructor because ElInfo is pure virtual.
    virtual ~ElInfo();

    /// Assignement operator.
    /// \param rhs right hand side.
    ElInfo& operator=(const ElInfo& rhs);

  public:
    /** \name getting methods
     * \{
     */

    /// Get ElInfo's \ref mesh
    Mesh* getMesh() const
    {
      return mesh;
    }

    /// Get ElInfo's \ref macroElement
    MacroElement* getMacroElement() const
    {
      return macroElement;
    }

    /// Get ElInfo's \ref element
    Element* getElement() const
    {
      return element;
    }

    /// Get ElInfo's \ref parent
    Element* getParent() const
    {
      return parent;
    }

    /// Get ElInfo's \ref fillFlag
    Flag getFillFlag() const
    {
      return fillFlag;
    }

    /// Get ElInfo's \ref level
    int getLevel() const
    {
      return level;
    }

    /// Get ElInfo's \ref iChild
    int getIChild() const
    {
      return iChild;
    }

    /// Get ElInfo's \ref coord[i]. This is a WorldVector<double> filled with
    /// the coordinates of the i-th vertex of element \ref el.
    WorldVector<double>& getCoord(int i)
    {
      return coord[i];
    }

    /// Get ElInfo's \ref coord[i]. This is a WorldVector<double> filled with the
    /// coordinates of the i-th vertex of element \ref el.
    WorldVector<double> const& getCoord(int i) const
    {
      return coord[i];
    }

    /// Get ElInfo's \ref coord. This is a FixVec<WorldVector<double> > filled
    /// with the coordinates of the all vertice of element \ref el.
    FixVec<WorldVector<double>, VERTEX>& getCoords()
    {
      return coord;
    }

    /// Get ElInfo's \ref coord. This is a FixVec<WorldVector<double> > filled
    /// with the coordinates of the all vertice of element \ref el.
    const FixVec<WorldVector<double>, VERTEX>& getCoords() const
    {
      return coord;
    }

    /// Get ElInfo's \ref oppCoord[i]
    WorldVector<double>& getOppCoord(int i)
    {
      return oppCoord[i];
    }

    WorldVector<double> const& getOppCoord(int i) const
    {
      return oppCoord[i];
    }

    /// Get ElInfo's \ref boundary[i]
    BoundaryType getBoundary(int i) const
    {
      return boundary[i];
    }

    /// Get boundary type of i-th vertex/edge/face (pos).
    BoundaryType getBoundary(GeoIndex pos, int i);

    FixVec<Element*, NEIGH> const& getNeighbours() const
    {
      return neighbour;
    }

    /// Get ElInfo's \ref neighbour[i]
    Element* getNeighbour(int i) const
    {
      return neighbour[i];
    }

    /// Get ElInfo's \ref neighbourCoord[i]
    FixVec<WorldVector<double>, VERTEX> getNeighbourCoord(int i) const
    {
      return neighbourCoord[i];
    }

    /// Get ElInfo's \ref oppVertex[i]
    int getOppVertex(int i) const
    {
      return oppVertex[i];
    }

    int getSideOfNeighbour(int i) const
    {
      return oppVertex[i];
    }

    /// Get ElInfo's \ref det
    double getDet() const
    {
      return det;
    }

    /// Returns \ref grdLambda
    const DimVec<WorldVector<double>>& getGrdLambda() const
    {
      return grdLambda;
    }

    void getGrdLambda(mtl::dense2D<double>& grd_lam);

    /// Returns \ref projection[i]
    Projection* getProjection(int i) const
    {
      return projection[i];
    }

    /// Returns \ref parametric
    bool getParametric() const
    {
      return parametric;
    }

    /// Returns \ref refinementPath
    unsigned long getRefinementPath() const
    {
      return refinementPath;
    }

    /// Get \ref refinementPathLength
    int getRefinementPathLength() const
    {
      return refinementPathLength;
    }

    /// Get ElInfo's \ref elType.
    int getType() const
    {
      return elType;
    }

    /// get transfomation matrix for multi-mesh ethod for given polynomial degree
    virtual mtl::dense2D<double>& getSubElemCoordsMat(int degree) const
    {
      return subElemMatrices[degree][std::make_pair(refinementPathLength, refinementPath)];
    }

    /** \} */

    /** \name setting methods
     * \{
     */

    /// Set ElInfo's \ref mesh
    void setMesh(Mesh* aMesh)
    {
      mesh = aMesh;
    }

    /// Set ElInfo's \ref macroElement
    void setMacroElement(MacroElement* mel)
    {
      macroElement = mel;
    }

    /// Set ElInfo's \ref element
    void setElement(Element* elem)
    {
      element = elem;
    }

    /// Set ElInfo's \ref parent
    void setParent(Element* elem)
    {
      parent = elem;
    }

    /// Set ElInfo's \ref fillFlag
    void setFillFlag(Flag flag)
    {
      fillFlag = flag;
    }

    /// Sets ElInfo's \ref coord[i].
    void setCoord(int i, WorldVector<double>& c)
    {
      coord[i] = c;
    }

    /// Sets ElInfo's \ref coord.
    void setCoords(FixVec<WorldVector<double>, VERTEX>& c)
    {
      coord = c;
    }

    /// Set ElInfo's \ref level
    void setLevel(int l)
    {
      level = l;
    }

    /// Set ElInfo's \ref boundary[i]
    void setBoundary(int i, BoundaryType t)
    {
      boundary[i] = newBound(boundary[i], t);
    }

    /// Set \ref projection[i] = p
    void setProjection(int i, Projection* p)
    {
      projection[i] = p;
    }

    /// Set ElInfo's \ref boundary[i]
    void setNeighbor(int i, Element* neigh)
    {
      neighbour[i] = neigh;
    }

    /// Set \ref det = d
    void setDet(double d)
    {
      det = d;
    }

    /// Set \ref parametric = param
    void setParametric(bool param)
    {
      parametric = param;
    }

    /// Set \ref refinementPath
    void setRefinementPath(unsigned long rPath)
    {
      refinementPath = rPath;
    }

    /// Set \ref refinementPathLength
    void setRefinementPathLength(int length)
    {
      refinementPathLength = length;
    }

    /// Set ElInfo's \ref elType.
    void setType(int t)
    {
      elType = t;
    }

    /** \} */


    /// Returns the absolute value of the determinant of the affine linear
    /// parametrization's Jacobian
    double calcDet() const;

    /// Used by non static method \ref calcDet(). Calculates the determinant
    /// for a given vector of vertex coordinates.
    double calcDet(const FixVec<WorldVector<double>, VERTEX>& coords) const;

    /// from CFE_Integration
    double calcSurfaceDet(VectorOfFixVecs<DimVec<double>>& surfVert) const;

    /// Checks whether flag is set in ElInfo's \ref fillFlag. If not, the
    /// program exits.
    void testFlag(const Flag& DBG_VAR(flag)) const
    {
      TEST_EXIT_DBG(fillFlag.isSet(flag))("flag not set\n");
    }

    /// Transforms local barycentric coordinates of a point defined on this
    /// element to global world coordinates.
    void coordToWorld(const DimVec<double>& lambda,
                      WorldVector<double>& world) const;

    /// Fills ElInfo's \ref det and \ref grdLambda entries.
    void fillDetGrdLambda();

    /// Returns a pointer to a vector, which contains the barycentric coordinates
    /// with respect to \ref element of a point with world coordinates world.
    /// The barycentric coordinates are stored in lambda.
    /// pure virtual => must be overriden in sub-class.
    virtual int worldToCoord(const WorldVector<double>& world,
                             DimVec<double>& lambda) const = 0;

    /// Fills this ElInfo with macro element information of mel.
    /// pure virtual => must be overriden in sub-class.
    virtual void fillMacroInfo(const MacroElement* mel) = 0;

    /// Fills this ElInfo for the child ichild using hierarchy information and
    /// parent data parentInfo.
    /// pure virtual => must be overriden in sub-class.
    virtual void fillElInfo(int ichild, const ElInfo* parentInfo) = 0;


    void fillElInfo(const MacroElement* mel,
                    int refinementPathLength,
                    unsigned long refinementPath);

    /// Calculates the Jacobian of the barycentric coordinates on \element and
    /// stores the matrix in grd_lam. The return value of the function is the
    /// absolute value of the determinant of the affine linear paraetrization's
    /// Jacobian.
    /// pure virtual => must be overriden in sub-class.
    virtual double calcGrdLambda(DimVec<WorldVector<double>>& grd_lam) = 0;

    /// calculates a normal of the given side (1d, 2d: edge, 3d: face) of \ref element.
    /// Returns the absolute value of the determinant of the
    /// transformation to the reference element.
    /// pure virtual => must be overriden in sub-class.
    virtual double getNormal(int side, WorldVector<double>& normal) const = 0;

    /// calculates a normal of the element in dim of world = dim + 1.
    /// Returns the absolute value of the determinant of the
    /// transformation to the reference element.
    /// pure virtual => must be overriden in sub-class.
    virtual double getElementNormal(WorldVector<double>& elementNormal) const = 0;

  protected:
    /// Pointer to the current mesh
    Mesh* mesh;

    /// Pointer to the current element
    Element* element;

    /// \ref element is child of element parent
    Element* parent;

    /// \ref element is an element of the binary tree located at MacroElement
    /// macroElement
    MacroElement* macroElement;

    /// Indicates wich elements will be called and wich information should be
    /// present while mesh traversal.
    Flag fillFlag;

    /// Level of the element. The level is zero for macro elements and the level
    /// of the children is (level of the parent + 1). level_ is filled always by
    /// the traversal routines.
    unsigned char level;

    /// Elements type index. This is used only for 3d, where type can be either 0, 1 or
    /// 2. In all other cases type is not used and the variable is set to 0.
    /// In 3d, it is filled automatically by the traversal routines.
    int elType;

    /// This ElInfo is the iChild-th child of the parent element.
    int iChild;

    /// \ref coord[i] is a WorldVector<double> storing the world coordinates of the
    /// i-th vertex of element \ref element.
    FixVec<WorldVector<double>, VERTEX> coord;

    /** \brief
     * boundary[i] is the BoundaryType of the i-th edge/face
     * for i=0,...,N_NEIGH - 1. In 3d
     * (*boundary)[N_FACES + i] is a pointer to the Boundary
     * object of the i-th edge, for i=0,..,N_EDGES - 1. It is
     * a NULL for an interior edge/face.
     */
    FixVec<BoundaryType, BOUNDARY> boundary;

    /// Vector storing pointers to projections for each face, edge, vertex.
    FixVec<Projection*, PROJECTION> projection;

    /// oppCoord[i] coordinates of the i-th neighbour vertex opposite the
    /// common edge/face.
    FixVec<WorldVector<double>, NEIGH> oppCoord;

    /// neighbour[i] pointer to the element at the edge/face with local index i.
    /// It is a NULL for boundary edges/faces.
    FixVec<Element*, NEIGH> neighbour;

    /// neighbourCoord[i][j] are the coordinate of the j-th vertex of the i-th
    /// neighbour element with the common edge/face.
    FixVec<FixVec<WorldVector<double>, VERTEX>, NEIGH> neighbourCoord;

    /// oppVertex[i] is undefined if neighbour[i] is a pointer to NULL.
    /// Otherwise it is the local index of the neighbour's vertex opposite the
    /// common edge/face.
    FixVec<int, NEIGH> oppVertex;

    /// Elements determinant.
    double det;

    /// Gradient of lambda.
    DimVec<WorldVector<double>> grdLambda;

    /// True, if this elInfo stores parametrized information. False, otherwise.
    bool parametric;

    /// Stores the world dimension.
    int dimOfWorld;

    unsigned long refinementPath;

    int refinementPathLength;

  public:
    static std::vector<std::map<std::pair<int, unsigned long>, mtl::dense2D<double>>> subElemMatrices;

    static std::vector<std::map<std::pair<int, unsigned long>, mtl::dense2D<double>>> subElemGradMatrices;

    /// child_vertex[el_type][child][i] = father's local vertex index of new
    /// vertex i. 4 stands for the newly generated vertex .
    static const int childVertex[3][2][4];

    /// child_edge[el_type][child][i] = father's local edge index of new edge i.
    /// new edge 2 is half of old edge 0, new edges 4,5 are really new edges, and
    /// value is different: child_edge[][][4,5] = index of same edge in other
    /// child.
    static const int childEdge[3][2][6];
  };

} // end namespace AMDiS

#include "ElInfo1d.hpp"
#include "ElInfo2d.hpp"
#include "ElInfo3d.hpp"
