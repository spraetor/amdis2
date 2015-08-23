/** \file Element.h */

#pragma once

#include "AMDiS_fwd.h"
#include "MatrixVector_fwd.h"
#include "Global.h"
#include "RefinementManager.h"
#include "ElementData.h"
#include "LeafData.h"

namespace AMDiS 
{
#define AMDIS_UNDEFINED  5 // used in MacroReader (???)

  /// Accessor class used as passkey pattern in the setIndex method
  /// to provide access to this specific method only.
  class MeshAccessor
  {
    friend class Mesh;
    MeshAccessor() = default;
    MeshAccessor(MeshAccessor const&) = default;
  };
  
  
  /** \ingroup Triangulation 
   * \brief
   * Base class for Line, Triangle, Tetrahedron
   *
   * Elements in AMDiS are always simplices (a simplex is a Line in 1d, a 
   * Triangle in 2d and a Tetrahedron in 3d). 
   * We restrict ourselves here to simplicial meshes, for several reasons:
   * -# A simplex is one of the most simple geometric types and complex domains 
   *    may be approximated by a set of simplices quite easily.
   * -# Simplicial meshes allow local refinement without the need of 
   *    nonconforming meshes (hanging nodes), parametric elements, or mixture of
   *    element types (which is the case for quadrilateral meshes).
   * -# Polynomials of any degree are easily represented on a simplex using 
   *    local (barycentric) coordinates.
   *
   * A Line element and its refinement:
   *
   * <img src="line.png">
   *
   * A Triangle element and its refinement:
   *
   * <img src="triangle.png">
   *
   * A Tetrahedron element and its refinements:
   *
   * <img src="tetrahedron.png">
   */
  class Element
  {
  private:
    /// private standard constructor because an Element must know his Mesh
    Element() {}

  public:
    /// constructs an Element which belongs to Mesh
    Element(Mesh *);

    /// copy constructor
    Element(const Element& old);

    /// destructor
    virtual ~Element();

    ///
    void deleteElementDOFs();

    /// Clone this Element and return a reference to it. Because also the DOFs
    /// are cloned, \ref Mesh::serializedDOfs must be used.
    Element* cloneWithDOFs(std::map<std::pair<DegreeOfFreedom, int>, DegreeOfFreedom*>& serializedDOFs);

    /** \name getting methods
     * \{
     */

    /// Returns \ref child[0]
    Element* getFirstChild() const 
    {
      return child[0];
    }

    /// Returns \ref child[1]
    Element* getSecondChild() const 
    {
      return child[1];
    }

    /// Returns \ref child[i], i=0,1
    Element* getChild(int i) const 
    {
      TEST_EXIT_DBG(i == 0 || i == 1)("There is only child 0 or 1! (i = %d)\n", i);
      return child[i];
    }

    /// Returns true if Element is a leaf element (\ref child[0] == NULL), returns
    /// false otherwise.
    bool isLeaf() const 
    { 
      return (child[0] == NULL); 
    }

    /// Returns \ref dof[i][j] which is the j-th DOF of the i-th node of Element.
    DegreeOfFreedom getDof(int i, int j) const 
    { 
      TEST_EXIT_DBG(dof != NULL)("DOFs are not valid in element %d!\n", index);

      return dof[i][j];
    }

    /// Returns \ref dof[i] which is a pointer to the DOFs of the i-th node.
    const DegreeOfFreedom* getDof(int i) const 
    {
      TEST_EXIT_DBG(dof != NULL)("DOFs are not valid in element %d!\n", index);

      return dof[i];
    }

    /// Returns a pointer to the DOFs of this Element
    const DegreeOfFreedom** getDof() const 
    {
      TEST_EXIT_DBG(dof != NULL)("DOFs are not valid in element %d!\n", index);

      return const_cast<const DegreeOfFreedom**>(dof);
    }

    /// Returns \ref mesh of Element
    Mesh* getMesh() const 
    { 
      return mesh; 
    }

    /// Returns \ref elementData's error estimation, if Element is a leaf element
    /// and has leaf data. 
    double getEstimation(int row) const
    {
      if (isLeaf()) {
        TEST_EXIT_DBG(elementData)("leaf element without leaf data\n");
        ElementData *ld = elementData->getElementData(ESTIMATABLE);
        TEST_EXIT_DBG(ld)("leaf data not estimatable!\n");

        return dynamic_cast<LeafDataEstimatableInterface*>(ld)->getErrorEstimate(row);
      }	
      
      return 0.0;
    }

    /// Returns Element's coarsening error estimation, if Element is a leaf 
    /// element and if it has leaf data and if this leaf data are coarsenable.
    double getCoarseningEstimation(int row) 
    {
      if (isLeaf()) {
        TEST_EXIT_DBG(elementData)("leaf element without leaf data\n");
        ElementData *ld = elementData->getElementData(COARSENABLE);
        TEST_EXIT_DBG(ld)("element data not coarsenable!\n");

        return dynamic_cast<LeafDataCoarsenableInterface*>(ld)->getCoarseningErrorEstimate(row);
      }
      
      return 0.0;
    }

    /// Returns region of element if defined, -1 else.
    int getRegion() const;
    
    /// Returns Element's \ref mark
    int getMark() const 
    { 
      return mark;
    }

    /// Returns \ref newCoord[i]
    double getNewCoord(int i) const 
    {
      TEST_EXIT_DBG(newCoord)("newCoord = NULL\n");
      return (*newCoord)[i];
    }

    /// Returns Element's \ref index
    int getIndex() const 
    { 
      return index; 
    }

    /// Returns \ref newCoord
    WorldVector<double>* getNewCoord() const 
    { 
      return newCoord; 
    }

    /** \} */

    /** \name setting methods
     * \{
     */

    /// Sets \ref child[0]
    virtual void setFirstChild(Element *aChild) 
    {
      child[0] = aChild;
    }

    /// Sets \ref child[1]
    virtual void setSecondChild(Element *aChild) 
    {
      child[1] = aChild;
    }

    /// Sets \ref elementData of Element
    void setElementData(ElementData* ed) 
    {
      elementData = ed;
    }

    /// Sets \ref newCoord of Element. Needed by refinement, if Element has a
    /// boundary edge on a curved boundary.
    void setNewCoord(WorldVector<double>* coord) 
    {
      newCoord = coord;
    }

    /// Sets \ref mesh.
    void setMesh(Mesh *m) 
    {
      mesh = m;
    }

    /// Sets the pointer to the DOFs of the i-th node of Element
    void setDof(int pos, DegreeOfFreedom* p) 
    {
      dof[pos] = p;
    }

    /// Delets the main DOF pointer. The caller must ensure that the DOFs are
    /// either freed correctly elsewhere or that they are still used by other
    /// elements.
    void delDofPtr()
    {
      delete [] dof;
      dof = NULL;
    }

    /// Checks whether Element is a leaf element and whether it has leaf data.
    /// If the checks don't fail, leaf data's error estimation is set to est.
    void setEstimation(double est, int row)
    {
      FUNCNAME("Element::setEstimation()");

      if (isLeaf()) {
        TEST_EXIT_DBG(elementData)("Leaf element %d without leaf data!\n", index);
        ElementData *ld = elementData->getElementData(ESTIMATABLE);
        TEST_EXIT_DBG(ld)("Leaf data %d not estimatable!\n", index);

        dynamic_cast<LeafDataEstimatableInterface*>(ld)->
          setErrorEstimate(row, est);
      } else {
        ERROR_EXIT("setEstimation only for leaf elements!\n");
      }
    }

    /// Sets Element's coarsening error estimation, if Element is a leaf element
    /// and if it has leaf data and if this leaf data are coarsenable.
    void setCoarseningEstimation(double est, int row)
    {
      if (isLeaf()) {
        TEST_EXIT_DBG(elementData)("leaf element without leaf data\n");
        ElementData *ld = elementData->getElementData(COARSENABLE);
        TEST_EXIT_DBG(ld)("leaf data not coarsenable\n");

        dynamic_cast<LeafDataCoarsenableInterface*>(ld)->
          setCoarseningErrorEstimate(row, est);
      } else {
        ERROR_EXIT("setEstimation only for leaf elements!\n");
      }
    }

    /// Sets Elements \ref mark = mark + 1;
    void incrementMark() 
    {
      mark++;
    }

    /// Sets Elements \ref mark = mark - 1;
    void decrementMark() 
    {
      if (0 < mark) 
        mark--;
    }

    /// Sets Element's \ref mark
    void setMark(int m) 
    {
      mark = m;
    }

    /** \} */


    /** \name pure virtual methods 
     * \{ 
     */

    /// Returns local vertex number of the j-th vertex of the i-th edge
    virtual int getVertexOfEdge(int i, int j) const = 0; 

    /// Returns local vertex number of the vertexIndex-th vertex of the
    /// positionIndex-th part of type position (vertex, edge, face)
    virtual int getVertexOfPosition(GeoIndex position,
                                    int positionIndex,
                                    int vertexIndex) const = 0;

    ///
    virtual int getPositionOfVertex(int side, int vertex) const = 0;

    ///
    virtual int getEdgeOfFace(int face, int edge) const = 0;

    ///
    virtual DofEdge getEdge(int localEdgeIndex) const = 0;

    ///
    virtual DofFace getFace(int localFaceIndex) const = 0;

    /// Returns the number of parts of type i in this element
    virtual int getGeo(GeoIndex i) const = 0;

    /// Orient the vertices of edges/faces. Used
    /// by Estimator for the jumps => same quadrature nodes from both sides!
    virtual void sortFaceIndices(int face, FixVec<int, WORLD> &vec) const = 0;

    /// Returns a copy of itself. Needed by Mesh to create Elements by 
    /// a prototype. 
    virtual Element *clone() const = 0;

    /// Returns which side of child[childnr] corresponds to side sidenr of 
    /// this Element. If the child has no corresponding side, the return value 
    /// is negative.
    virtual int getSideOfChild(int childnr, int sidenr, int elType = 0) const = 0;

    /** \brief
     * Generalization of \ref getSideOfChild to arbitrary subObject. Thus, 
     * e.g., in 3d we can ask for the local id of a verte, edge or face 
     * on the elements children.
     *
     * \param[in]  childnr    Either 0 or 1 for the left or right children.
     * \param[in]  subObj     Defines whether we ask for VERTEX, EDGE or FACE.
     * \param[in]  ithObj     Number of the object on the parent.
     * \param[in]  elType     Type of the element. Important only in 3D.
     */
    virtual int getSubObjOfChild(int childnr, GeoIndex subObj, int ithObj, 
                                 int elType = 0) const = 0;

    /// Returns which vertex of elements parent corresponds to the vertexnr of
    /// the element, if the element is the childnr-th child of the parent.
    /// If the vertex is the ner vertex at the refinement edge, -1 is returned.
    virtual int getVertexOfParent(int childnr, 
                                  int vertexnr, 
                                  int elType = 0) const = 0;

    /// Returns whether Element is a Line
    virtual bool isLine() const = 0;

    /// Returns whether Element is a Triangle
    virtual bool isTriangle() const = 0;

    /// Returns whether Element is a Tetrahedron
    virtual bool isTetrahedron() const = 0;

    /// Returns whether Element has sideElem as one of its sides.
    virtual bool hasSide(Element *sideElem) const = 0;

    /** \brief
     * Returns for a given element type number the element type number of the children.
     * For 1d and 2d this is always 0, because element type number are used in the 
     * 3d case only.
     */
    virtual int getChildType(int elType) const = 0;

    /** \brief
     * Traverses a vertex/edge/face of a given element (this includes also all
     * children of the element having the same edge/face). All DOFs on mesh
     * nodes alonge this vertex/edge/face are assembled and put together to 
     * a list.
     *
     * \param[in]  feSpace     FE space which is used to get the dofs.
     * \param[in]  bound       Defines the vertex/edge/face of the element on
     *                         which all vertex dofs are assembled.
     * \param[out] dofs        List of dofs, where the result is stored.
     * \param[in]  baseDofPtr  If true, the base DOF pointes are stored. Thus,
     *                         dof* [\ref dof] of the element is inserted. If 
     *                         false, &(dof[.][n0]) is put to the result vector, 
     *                         with n0 beging the number of predofs.
     */
    virtual void getNodeDofs(const FiniteElemSpace* feSpace, 
                             BoundaryObject bound,
                             DofContainer& dofs,
                             bool baseDofPtr = false) const = 0;

    /** \brief
     * Traverses a vertex/edge/face of a given element (this includes also all
     * children of the element having the same edge/face). All DOFs belonging
     * to higher order basis functions alonge this vertex/edge/face are 
     * assembled and put together to a list.
     *
     * \param[in]  feSpace     FE space which is used to get the dofs.
     * \param[in]  bound       Defines the edge/face of the element on which
     *                         all non vertex dofs are assembled.
     * \param[out] dofs        All dofs are put to this dof list.
     * \param[in]  baseDofPtr  If true, the base DOF pointes are stored. Thus,
     *                         dof* [\ref dof] of the element is inserted. If 
     *                         false, &(dof[.][n0]) is put to the result vector, 
     *                         with n0 beging the number of predofs.
     * \param[out] dofGeoIndex Optional, the function can store to each DOF in
     *                         the DofContainer dofs the geometric index, thus
     *                         identifing the DOF to be a vertex, edge, face or
     *                         center DOF.
     */
    virtual void getHigherOrderDofs(const FiniteElemSpace* feSpace, 
                                    BoundaryObject bound,
                                    DofContainer& dofs,
                                    bool baseDofPtr = false,
                                    std::vector<GeoIndex>* dofGeoIndex = NULL) const = 0;

    virtual void getSubBoundary(BoundaryObject bound, 
                                std::vector<BoundaryObject> &subBound) const = 0;

    /** \} */

    // ===== other public methods =================================================

    /// Combines \ref getNodeDofs and \ref getHigherOrderDofs to one function. 
    /// See parameter description there.
    void getAllDofs(const FiniteElemSpace* feSpace, 
                    BoundaryObject bound, 
                    DofContainer& dofs,
                    bool baseDofPtr = false,
                    std::vector<GeoIndex>* dofGeoIndex = NULL);
    
    /// assignment operator
    Element& operator=(const Element& el);

    /// Checks whether the face with vertices dof[0],..,dof[DIM-1] is
    /// part of mel's boundary. returns the opposite vertex if true, -1 else
    int oppVertex(FixVec<DegreeOfFreedom*, DIMEN> pdof) const;

    /// Refines Element's leaf data
    void refineElementData(Element* child1, Element* child2, int elType = 0) 
    {
      if (elementData) {
        bool remove = elementData->refineElementData(this, child1, child2, elType);
        if (remove) {
          ElementData *tmp = elementData->getDecorated();
          delete elementData;
          elementData = tmp;
        }
      }
    }

    /// Coarsens Element's leaf data
    void coarsenElementData(Element* child1, Element* child2, 
				   int elType = 0) 
    {
      ElementData *childData;
      childData = child1->getElementData();
      if (childData) {
        childData->coarsenElementData(this, child1, child2, elType);
        delete childData;
        child1->setElementData(NULL);
      }
      childData = child2->getElementData();
      if (childData) {
        childData->coarsenElementData(this, child2, child1, elType);
        delete childData;
        child2->setElementData(NULL);
      }
    }

    /// Returns pointer to \ref elementData
    ElementData* getElementData() const 
    {
      return elementData;
    }

    ///
    ElementData* getElementData(int typeID) const 
    {
      if (elementData)
        return elementData->getElementData(typeID);

      return NULL;
    }

    /// Deletes the \ref elementData with a specific typeID.
    bool deleteElementData(int typeID);

    /** \brief
     * Returns whether element is refined at side side
     * el1, el2 are the corresponding children. 
     * (not neccessarly the direct children!)
     * elementTyp is the type of this element (comes from ElInfo)
     */
    bool isRefinedAtSide(int side, Element *el1, Element *el2, 
                         unsigned char elementTyp = 255);

    /// Returns whether Element's \ref newCoord is set
    bool isNewCoordSet() const 
    { 
      return (newCoord != NULL);
    }

    /// Frees memory for \ref newCoord
    void eraseNewCoord();

    /// Sets Element's \ref dof pointer.
    void createNewDofPtrs(bool setDofs = false);
    
  /* protected: */
    
    // this method can only be accessed by a friend of MeshAccessor
    void setIndex(MeshAccessor, int i) 
    {
      index = i;
    }
    
    static void clearDeletedDofs(MeshAccessor)
    {
      Element::deletedDOFs.clear();
    }

    /// Used by friend class Mesh while dofCompress
    void newDofFct1(MeshAccessor, const DOFAdmin*, std::vector<DegreeOfFreedom>&);

    /// Used by friend class Mesh while dofCompress
    void newDofFct2(MeshAccessor, const DOFAdmin*);
 
  protected:
    /// Changes old dofs to negative new dofs
    void changeDofs1(const DOFAdmin* admin, std::vector<DegreeOfFreedom>& newDofIndex,
                     int n0, int nd0, int nd, int pos);

    /// Changes negative new dofs to positive
    void changeDofs2(int n0, int nd0, int nd, int pos);

  protected:
    /// Pointers to the two children of interior elements of the tree. Pointers
    /// to NULL for leaf elements.
    Element *child[2];

    /// Vector of pointers to DOFs. These pointers must be available for elements
    /// vertices (for the geometric description of the mesh). There my be pointers
    /// for the edges, for faces and for the center of an element. They are 
    /// ordered the following way: The first N_VERTICES entries correspond to the
    /// DOFs at the vertices of the element. The next ones are those at the edges,
    /// if present, then those at the faces, if present, and then those at the 
    /// barycenter, if present.
    DegreeOfFreedom **dof;

    /// Unique global index of the element. these indices are not strictly ordered
    /// and may be larger than the number of elements in the binary tree (the list
    /// of indices may have holes after coarsening).
    int index;

    /// Marker for refinement and coarsening. if mark is positive for a leaf
    /// element, this element is refined mark times. if mark is negative for
    /// a leaf element, this element is coarsened -mark times.
    int mark;
 
    /// If the element has a boundary edge on a curved boundary, this is a pointer
    /// to the coordinates of the new vertex that is created due to the refinement
    /// of the element, otherwise it is a NULL pointer. Thus coordinate 
    /// information can be also produced by the traversal routines in the case of 
    /// curved boundary.
    WorldVector<double> *newCoord;

    /// Pointer to the Mesh this element belongs to
    Mesh* mesh;

    /// Pointer to Element's leaf data
    ElementData* elementData;

    /// This map is used for deletion of all DOFs of all elements of a mesh. Once
    /// a DOF-vector (all DOFS at a node, edge, etc.) is deleted, its address is
    /// added to this map to note not to delete it a second time.
    static std::map<DegreeOfFreedom*, bool> deletedDOFs;
  };


  /// Writes the element hierarchie to a Graphviz dot-file. Using the dot-tool from
  /// Graphvis, this dot-file can be converted to a ps-file. Useful for debugging!
  void writeDotFile(Element *el, std::string filename, int maxLevels = -1);
  
} // end namespace AMDiS
