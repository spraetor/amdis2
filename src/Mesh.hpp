/** \defgroup Triangulation Triangulation module
 * @{ <img src="triangulation.png"> @}
 *
 * Example:
 *
 * @{ <img src="hierarchicalMesh.png"> @}
 *
 * \brief
 * Contains all triangulation classes.
 */

#pragma once

#include <deque>
#include <set>
#include <stdio.h>

#include "AMDiS_fwd.hpp"
#include "BoundaryCondition.hpp"
#include "DOFAdmin.hpp"
#include "ElInfo.hpp"
#include "Element.hpp"
#include "FixVec.hpp"
#include "Line.hpp"
#include "Tetrahedron.hpp"
#include "Triangle.hpp"

namespace AMDiS
{
  /** \ingroup Triangulation
   * \brief
   * A Mesh holds all information about a triangulation.
   */
  class Mesh
  {
  public:
    /// Creates a mesh with the given name of dimension dim
    Mesh(std::string name, int dim);

    /// Destructor
    ~Mesh();

    /// Reads macro triangulation.
    void initialize();

    /// Assignment operator
    Mesh& operator=(Mesh const&);

    /** \name getting methods
     * \{
     */

    /// Returns geometric information about this mesh. With GeoIndex p it is
    /// specified which information is requested.
    int getGeo(GeoIndex p) const
    {
      return Global::getGeo(p, dim);
    }

    /// Returns \ref name of the mesh
    std::string getName() const
    {
      return name;
    }

    /// Returns \ref dim of the mesh
    int getDim() const
    {
      return dim;
    }

    /// Returns \ref nDofEl of the mesh
    int getNumberOfAllDofs() const
    {
      return nDofEl;
    }

    /// Returns \ref nNodeEl of the mesh
    int getNumberOfNodes() const
    {
      return nNodeEl;
    }

    /// Returns \ref nVertices of the mesh
    int getNumberOfVertices() const
    {
      return nVertices;
    }

    /// Returns \ref nEdges of the mesh
    int getNumberOfEdges() const
    {
      return nEdges;
    }

    /// Returns \ref nFaces of the mesh
    int getNumberOfFaces() const
    {
      return nFaces;
    }

    /// Returns \ref nLeaves of the mesh
    int getNumberOfLeaves() const
    {
      return nLeaves;
    }

    /// Returns \ref nElements of the mesh
    int getNumberOfElements() const
    {
      return nElements;
    }

    /// Returns \ref maxEdgeNeigh of the mesh
    int getMaxEdgeNeigh() const
    {
      return maxEdgeNeigh;
    }

    /// Returns \ref parametric of the mesh
    Parametric* getParametric() const
    {
      return parametric;
    }

    /// Returns \ref diam of the mesh
    WorldVector<double> const& getDiameter() const
    {
      return diam;
    }

    /// Returns \ref diam of the mesh
    std::pair<WorldVector<double>, WorldVector<double>> const& getBoundingBox() const
    {
      return boundingBox;
    }

    /// Returns nDof[i] of the mesh
    int getNumberOfDofs(int i) const
    {
      TEST_EXIT_DBG(i <= dim)("Wrong index: %d %d\n", i, dim);
      return nDof[i];
    }

    /// Returns \ref elementPrototype of the mesh
    Element* getElementPrototype()
    {
      return elementPrototype;
    }

    /// Returns \ref leafDataPrototype of the mesh
    ElementData* getElementDataPrototype()
    {
      return elementDataPrototype;
    }

    /// Returns node[i] of the mesh
    int getNode(int i) const
    {
      return node[i];
    }

    /// Allocates the number of DOFs needed at position and registers the DOFs
    /// at the DOFAdmins. The number of needed DOFs is the sum over the needed
    /// DOFs of all DOFAdmin objects belonging to this mesh.
    /// The return value is a pointer to the first allocated DOF.
    DegreeOfFreedom* getDof(GeoIndex position);

    /// Returns *(\ref admin[i]) of the mesh
    DOFAdmin const& getDofAdmin(int i) const
    {
      return *(admin[i]);
    }

    /// Creates a DOFAdmin with name lname. nDof specifies how many DOFs
    /// are needed at the different positions (see \ref DOFAdmin::nrDOF).
    /// A pointer to the created DOFAdmin is returned.
    DOFAdmin const* createDOFAdmin(std::string lname, DimVec<int> nDof);

    /// Returns the size of \ref admin which is the number of the DOFAdmins
    /// belonging to this mesh
    int getNumberOfDOFAdmin() const
    {
      return int(admin.size());
    }

    /// Returns the size of \ref macroElements which is the number of
    /// of macro elements of this mesh
    int getNumberOfMacros() const
    {
      return int(macroElements.size());
    }

    /// Returns a DOFAdmin which at least manages vertex DOFs
    DOFAdmin const* getVertexAdmin() const;

    /// Allocates an array of DOF pointers. The array holds one pointer for
    /// each node.
    DegreeOfFreedom** createDofPtrs();

    /// Returns \ref preserveCoarseDOFs of the mesh
    bool queryCoarseDOFs() const
    {
      return preserveCoarseDOFs;
    }

    /// Returns an iterator to the begin of \ref macroElements
    std::deque<MacroElement*>::iterator firstMacroElement()
    {
      return macroElements.begin();
    }

    /// Returns macroElements[i].
    MacroElement* getMacroElement(int i)
    {
      return macroElements[i];
    }

    /// Returns an iterator to the end of \ref macroElements
    std::deque<MacroElement*>::iterator endOfMacroElements()
    {
      return macroElements.end();
    }

    /// Returns \ref macroElements, the list of all macro elements in the mesh.
    std::deque<MacroElement*>& getMacroElements()
    {
      return macroElements;
    }

    /** \} */

    /** \name setting methods
     * \{
     */

    /// Sets \ref name of the mesh
    void setName(std::string aName)
    {
      name = aName;
    }

    /// Sets \ref nVertices of the mesh
    void setNumberOfVertices(int n)
    {
      nVertices = n;
    }

    /// Sets \ref nFaces of the mesh
    void setNumberOfFaces(int n)
    {
      nFaces = n;
    }

    /// Increments \ref nVertices by inc
    void incrementNumberOfVertices(int inc)
    {
      nVertices += inc;
    }

    /// Sets \ref nEdges of the mesh
    void setNumberOfEdges(int n)
    {
      nEdges = n;
    }

    /// Increments \ref nEdges by inc
    void incrementNumberOfEdges(int inc)
    {
      nEdges += inc;
    }

    /// Increments \ref nFaces by inc
    void incrementNumberOfFaces(int inc)
    {
      nFaces += inc;
    }

    /// Sets \ref nLeaves of the mesh
    void setNumberOfLeaves(int n)
    {
      nLeaves = n;
    }

    /// Increments \ref nLeaves by inc
    void incrementNumberOfLeaves(int inc)
    {
      nLeaves += inc;
    }

    /// Sets \ref nElements of the mesh
    void setNumberOfElements(int n)
    {
      nElements = n;
    }

    /// Increments \ref nElements by inc
    void incrementNumberOfElements(int inc)
    {
      nElements += inc;
    }

    /// Sets *\ref diam to w
    void setDiameter(WorldVector<double> const& w);

    /// Sets (*\ref diam)[i] to d
    void setDiameter(int i, double d);
    
    void setBoundingBox(WorldVector<double> const& min_corner, 
                        WorldVector<double> const& max_corner);

    /// Sets \ref preserveCoarseDOFs = true
    void retainCoarseDOFs()
    {
      preserveCoarseDOFs = true;
    }

    /// Sets \ref preserveCoarseDOFs = b
    void setPreserveCoarseDOFs(bool b)
    {
      preserveCoarseDOFs = b;
    }

    /// Sets \ref preserveCoarseDOFs = false
    void noCoarseDOFs()
    {
      preserveCoarseDOFs = false;
    }

    /// Sets \ref elementPrototype of the mesh
    void setElementPrototype(Element* prototype)
    {
      elementPrototype = prototype;
    }

    /// Sets \ref elementDataPrototype of the mesh
    void setElementDataPrototype(ElementData* prototype)
    {
      elementDataPrototype = prototype;
    }

    ///
    void setParametric(Parametric* param)
    {
      parametric = param;
    }

    ///
    void setMaxEdgeNeigh(int m)
    {
      maxEdgeNeigh = m;
    }

    /** \} */

    /// Creates a new Element by cloning \ref elementPrototype
    Element* createNewElement(Element* parent = NULL);

    /// Creates a new ElInfo dependent of \ref dim of the mesh
    ElInfo* createNewElInfo();

    /// Frees DOFs at the given position pointed by dof
    void freeDof(DegreeOfFreedom* dof, GeoIndex position);

    /// Frees memory for the given element el
    void freeElement(Element* el);

    /// Performs DOF compression for all DOFAdmins (see \ref DOFAdmin::compress)
    void dofCompress();

    /// Adds a DOFAdmin to the mesh
    void addDOFAdmin(DOFAdmin* admin);

    /// Recalculates the number of leave elements.
    void updateNumberOfLeaves();

    /// Clears \ref macroElements
    void clearMacroElements()
    {
      macroElements.clear();
    }

    /// Adds a macro element to the mesh
    void addMacroElement(MacroElement* me);

    /// Removes a set of macro elements from the mesh. This works only for the
    /// case, that there are no global or local refinements, i.e., all macro
    /// elements have no children.
    void removeMacroElements(std::set<MacroElement*>& macros,
                             std::vector<FiniteElemSpace const*>& feSpaces);

    void removeAllMacroElements();

    /// Frees the array of DOF pointers (see \ref createDofPtrs)
    void freeDofPtrs(DegreeOfFreedom** ptrs);

    /// Used by \ref findElementAtPoint.
    bool findElInfoAtPoint(WorldVector<double> const& xy,
                           ElInfo* el_info,
                           DimVec<double>& bary,
                           MacroElement const* start_mel,
                           WorldVector<double> const* xy0,
                           double* sp);

    /** \brief
     * Access to an element at world coordinates xy. Some applications need the
     * access to elements at a special location in world coordinates. Examples
     * are characteristic methods for convection problems, or the implementation
     * of a special right hand side like point evaluations or curve integrals.
     * For such purposes, a routine is available which returns an element pointer
     * and corresponding barycentric coordinates.
     *
     * \param xy world coordinates of point
     * \param elp return address for a pointer to the element at xy
     * \param pary returns barycentric coordinates of xy
     * \param start_mel initial guess for the macro element containing xy or NULL
     * \param xy0 start point from a characteristic method, see below, or NULL
     * \param sp return address for relative distance to domain boundary in a
     *        characteristic method, see below, or NULL
     * \return true is xy is inside the domain , false otherwise
     *
     * For a characteristic method, where \f$ xy = xy_0 - V\tau \f$, it may be
     * convenient to know the point on the domain's boundary which lies on the
     * line segment between the old point xy0 and the new point xy, in case that
     * xy is outside the domain. Such information is returned when xy0 and a
     * pointer sp!=NULL are supplied: *sp is set to the value s such that
     * \f$ xy_0 +s (xy -xy_0) \in \partial Domain \f$, and the element and local
     * coordinates corresponding to that boundary point will be returned via elp
     * and bary.
     *
     * The implementation of findElementAtPoint() is based on the transformation
     * from world to local coordinates, available via the routine worldToCoord(),
     * At the moment, findElementAtPoint() works correctly only for domains with
     * non-curved boundary. This is due to the fact that the implementation first
     * looks for the macro-element containing xy and then finds its path through
     * the corresponding element tree based on the macro barycentric coordinates.
     * For non-convex domains, it is possible that in some cases a point inside
     * the domain is considered as external.
     */
    bool findElementAtPoint(WorldVector<double> const& xy,
                            Element** elp,
                            DimVec<double>& bary,
                            MacroElement const* start_mel,
                            WorldVector<double> const* xy0,
                            double* sp);

    /** \brief
     * Returns for a given dof its world coordinates in this mesh. Because we do
     * not have any direct connection between dofs and coordinates, this function
     * has to search for the element in this mesh, that contains the dof. Than the
     * coordinates can be computed. Therefore, this function is very costly and
     * should be used for debugging purpose only.
     *
     * @param[in]    dof       A pointer to the dof we have to search for.
     * @param[in]    feSpace   The fe space to be used for the search.
     * @param[out]   coords    World vector that stores the coordinates of the dof.
     *
     * The function returns true, if the dof was found, otherwise false.
     */
    bool getDofIndexCoords(DegreeOfFreedom const* dof,
                           FiniteElemSpace const* feSpace,
                           WorldVector<double>& coords)
    {
      return getDofIndexCoords(*dof, feSpace, coords);
    }


    /// This function is equal to \ref getDofIndexCoords as defined above, but
    /// takes a DOF index instead of a DOF pointer.
    bool getDofIndexCoords(DegreeOfFreedom dof,
                           FiniteElemSpace const* feSpace,
                           WorldVector<double>& coords);

    /** \brief
     * Traverse the whole mesh and stores to each DOF the coordinates in a given
     * DOFVector. Works in the same way as the function \ref getDofIndexCoords
     * defined above.
     *
     * @param[out]  coords    DOF vector that stores the coordinates to each DOF.
     */
    void getDofIndexCoords(DOFVector<WorldVector<double>>& coords);


    /** \brief
     * Traverse the mesh and get all DOFs in this mesh for a given FE space.
     *
     * @param[in]   feSpace   The FE space to be used for collecting DOFs.
     * @param[out]  allDofs   The set which is filled with all DOFs.
     */
    void getAllDofs(FiniteElemSpace const* feSpace,
                    std::set<DegreeOfFreedom const*>& allDofs);

    /// Returns FILL_ANY_?D
    static Flag const& getFillAnyFlag(int dim)
    {
      switch (dim)
      {
      case 1:
        return FILL_ANY_1D;
        break;
      case 2:
        return FILL_ANY_2D;
        break;
      case 3:
        return FILL_ANY_3D;
        break;
      default:
        ERROR_EXIT("invalid dim\n");
        return FILL_ANY_1D;
      }
    }

    /// Returns \ref elementIndex and increments it by 1.
    int getNextElementIndex()
    {
      return elementIndex++;
    }

    /// Returns \ref initialized.
    bool isInitialized() const
    {
      return initialized;
    }

    ///
    std::map<BoundaryType, VertexVector*>& getPeriodicAssociations()
    {
      return periodicAssociations;
    }

    /// Returns the periodic association for a specific boundary type.
    VertexVector& getPeriodicAssociations(BoundaryType b)
    {
      FUNCNAME_DBG("Mesh::getPeriodicAssociations()");

      TEST_EXIT_DBG(periodicAssociations.count(b) == 1)
      ("There are no periodic assoications for boundary type %d!\n", b);

      return (*(periodicAssociations[b]));
    }

    void setPeriodicAssociations(BoundaryType b, VertexVector* vec)
    {
      periodicAssociations[b] = vec;
    }


    /// Returns whether the given boundary type is periodic, i.e., if there is
    /// a periodic association for this boundary type.
    bool isPeriodicAssociation(BoundaryType b) const
    {
      return (periodicAssociations.count(b) == 1 ? true : false);
    }

    ///
    bool associated(DegreeOfFreedom dof1, DegreeOfFreedom dof2);

    ///
    bool indirectlyAssociated(DegreeOfFreedom dof1, DegreeOfFreedom dof2);

    /// Returns \macroFileInfo
    MacroInfo* getMacroFileInfo()
    {
      return macroFileInfo;
    }

    /// Increment the value of mesh change index, see \ref changeIndex.
    void incChangeIndex()
    {
      changeIndex++;
    }

    /// Returns the mesh change index, see \ref changeIndex.
    long getChangeIndex() const
    {
      return changeIndex;
    }

    ///
    void clearMacroFileInfo();

    ///
    int calcMemoryUsage();

    ///
    void deleteMeshStructure();

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    /// In parallel computations the level of all macro elements is equal to the
    /// number of global pre refinements, \ref nParallelPreRefinements.
    int getMacroElementLevel() const
    {
      return nParallelPreRefinements;
    }
#else
    /// In sequentiel computations the level of all macro elements is always 0.
    int getMacroElementLevel() const
    {
      return 0;
    }
#endif

    /// Creates a map for all elements in mesh that maps from element indices
    /// to the corresponding pointers.
    void getElementIndexMap(std::map<int, Element*>& elIndexMap);

  public:
    ///
    static const Flag FILL_NOTHING;

    ///
    static const Flag FILL_COORDS;

    ///
    static const Flag FILL_BOUND;

    ///
    static const Flag FILL_NEIGH;

    ///
    static const Flag FILL_OPP_COORDS;

    ///
    static const Flag FILL_ORIENTATION;

    ///
    static const Flag FILL_ADD_ALL;

    ///
    static const Flag FILL_ANY_1D;

    ///
    static const Flag FILL_ANY_2D;

    ///
    static const Flag FILL_ANY_3D;

    ///
    static const Flag FILL_DET;

    ///
    static const Flag FILL_GRD_LAMBDA;

    //**************************************************************************
    //  flags for Mesh traversal
    //**************************************************************************

    ///
    static const Flag CALL_EVERY_EL_PREORDER;

    ///
    static const Flag CALL_EVERY_EL_INORDER;

    ///
    static const Flag CALL_EVERY_EL_POSTORDER;

    ///
    static const Flag CALL_LEAF_EL;

    ///
    static const Flag CALL_LEAF_EL_LEVEL;

    ///
    static const Flag CALL_EL_LEVEL;

    ///
    static const Flag CALL_MG_LEVEL;

    /// If set, left and right children are swapped in traverse.
    static const Flag CALL_REVERSE_MODE;

  protected:
    ///
    bool findElementAtPointRecursive(ElInfo* elinfo,
                                     DimVec<double> const& lambda,
                                     int outside,
                                     ElInfo* final_el_info);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    /** \brief
     * This functions is called in parallel computations by the function \ref
     * Mesh::initialize(). It checks that the macro file has enough macro elements
     * for the number of used processors and that all macro elements are of type 0.
     * If this is not the case, that macro mesh is globally refined in an
     * apropriate way and is written to a new macro file.
     *
     * The function overwrittes the macro and periodic filenames, if a new macro
     * fule was created for the current parallel usage.
     *
     * \param[in/out]  macroFilename      Name of the macro mesh file.
     * \param[in/out]  periodicFilename   If periodic boundaries are used, name of the
     *                                    periodicity file. Otherwise, the string must
     *                                    be empty.
     * \param[in]      check              If the mesh should be checked to be a correct
     *                                    AMDiS macro mesh, the value must be 1 and 0
     *                                    otherwise.
     */
    void checkParallelMacroFile(std::string& macroFilename,
                                std::string& periodicFilename,
                                int check);
#endif

  protected:
    /// maximal number of DOFs at one position
    static const int MAX_DOF;

    /// Name of this Mesh
    std::string name;

    /// Dimension of this Mesh. Doesn't have to be equal to dimension of world.
    int dim;

    /// Number of vertices in this Mesh
    int nVertices;

    /// Number of Edges in this Mesh
    int nEdges;

    /// Number of leaf elements in this Mesh
    int nLeaves;

    /// Total number of elements in this Mesh
    int nElements;

    /// Number of faces in this Mesh
    int nFaces;

    /// Maximal number of elements that share one edge; used to allocate memory
    /// to store pointers to the neighbour at the refinement/coarsening edge
    /// (only 3d);
    int maxEdgeNeigh;

    /// Diameter of the mesh in the DIM_OF_WORLD directions
    WorldVector<double> diam;
    
    /// Corner coordinates of the bounding-box of the mesh
    std::pair<WorldVector<double>, WorldVector<double>> boundingBox;

    /// Is pointer to NULL if mesh contains no parametric elements else pointer
    /// to a Parametric object containing coefficients of the parameterization
    /// and related information
    Parametric* parametric;

    /// When an element is refined, not all dofs of the coarse element must be
    /// part of the new elements. An example are centered dofs when using higher
    /// lagrange basis functions. The midpoint dof of the parents element is not
    /// a dof of the both children elements. Therefore, the dof can be deleted.
    /// In some situation, e.g., when using multigrid techniques, it can be
    /// necessary to store this coarse dofs. Then this variable must be set to
    /// true. If false, the not required coarse dofs will be deleted.
    bool preserveCoarseDOFs;

    /// Number of all DOFs on a single element
    int nDofEl;

    /** \brief
     * Number of DOFs at the different positions VERTEX, EDGE, (FACE,) CENTER on
     * an element:
     *
     * - nDof[VERTEX]: number of DOFs at a vertex (>= 1)
     *
     * - nDof[EDGE]: number of DOFs at an edge; if no DOFs are associated to
     *   edges, then this value is 0
     *
     * - nDof[FACE]: number of DOFs at a face; if no DOFs are associated to
     *   faces, then this value is 0 (only 3d)
     *
     * - nDof[CENTER]: number of DOFs at the barycenter; if no DOFs are
     *   associated to the barycenter, then this value is 0
     */
    DimVec<int> nDof;

    /// Number of nodes on a single element where DOFs are located. Needed for
    /// the (de-) allocation of the DOF-vector on the element (\ref Element::dof).
    /// Here "node" is equivalent to the number of basis functions on the element.
    int nNodeEl;

    /** \brief
     * Gives the index of the first node at vertex, edge, face (only 3d), and
     * barycenter:
     *
     * - node[VERTEX]: has always value 0; dof[0],...,dof[N_VERTICES-1] are
     *   always DOFs at the vertices;
     *
     * - node[EDGE]: dof[node[EDGE]],..., dof[node[EDGE]+N_EDGES-1] are the DOFs
     *   at the N_EDGES edges, if DOFs are located at edges;
     *
     * - node[FACE]: dof[node[FACE]],..., dof[node[FACE]+N_FACES-1] are the DOFs
     *   at the N_FACES faces, if DOFs are located at faces (only 3d);
     *
     * - node[CENTER]: dof[node[CENTER]] are the DOFs at the barycenter, if DOFs
     *   are located at the barycenter;
     */
    DimVec<int> node;

    /// List of all DOFAdmins
    std::vector<DOFAdmin*> admin;

    /// List of all MacroElements of this Mesh
    std::deque<MacroElement*> macroElements;

    /// Used by check functions
    static std::vector<DegreeOfFreedom> dof_used;

    static std::set<std::string> refinedMeshNames;

    /// This map is used for serialization and deserialization of mesh elements.
    /// During the serialization process, all elements are visited and their
    /// DOF indices are written to the file. If a dof index at a position, i.e.
    /// vertex, line or face, was written to file, the combination of dof index
    /// and position is inserted to this map. That ensures that the same dof at
    /// the same position, but being part of another element, is not written
    /// twice to the file. When a state should be deserialized, the information
    /// can be used to construct exactly the same dof structure.
    static std::map<std::pair<DegreeOfFreedom, int>, DegreeOfFreedom*> serializedDOFs;

    /// Used while mesh refinement. To create new elements
    /// elementPrototype->clone() is called, which returns a Element of the
    /// same type as elementPrototype. So e.g. Elements of the different
    /// dimensions can be created in a uniform way.
    Element* elementPrototype;

    /// Prototype for leaf data. Used for creation of new leaf data while
    /// refinement.
    ElementData* elementDataPrototype;

    /// Used for enumeration of all mesh elements
    int elementIndex;

    /// True if the mesh is already initialized, false otherwise.
    bool initialized;

    /// Map of managed periodic vertex associations.
    std::map<BoundaryType, VertexVector*> periodicAssociations;

    /// If the mesh has been created by reading a macro file, here the information
    /// are stored about the content of the file.
    MacroInfo* macroFileInfo;

    /// This index is incremented every time the mesh is changed, e.g. by the
    /// refinement or the coarsening manager. It can be used by other object if
    /// the mesh has been changed by first copying this variable elsewhere and
    /// comparing its values.
    long changeIndex;

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    /// In parallel computations the mesh may be globally prerefined to achieve a
    /// fine enought starting mesh for the given number of ranks. The value of the
    /// variable will be defined in function \ref checkParallelMacroFile.
    int nParallelPreRefinements;
#endif

  protected:
    /// for findElement-Fcts
    DimVec<double> final_lambda;

    /// Temporary variables that are used in functions \ref findElInfoAtPoint
    /// and \ref findElementAtPointRecursive.
    WorldVector<double> const *g_xy0, *g_xy;

    /// Temporary variable that is used in functions \ref findElInfoAtPoint and
    /// \ref findElementAtPointRecursive.
    double* g_sp;

    friend class MacroInfo;
    friend class io::MacroReader;
    friend struct io::MacroWriter;
  };


  /// Accessor class used as passkey pattern in the setIndex method
  /// to provide access to this specific method only.
  class MeshAccessor
  {
    MeshAccessor() = default;
    MeshAccessor(MeshAccessor const&) = default;

    friend class Mesh;
  };

  class BoundaryWrapper
  {
  public:
    BoundaryWrapper(BoundaryType nr, Mesh* mesh = NULL)
      : nr(nr),
        mesh(mesh)
    {}
    
    Mesh* getMesh() 
    { 
      return mesh; 
    }
    
    Mesh const* getMesh() const
    { 
      return mesh; 
    }
    
    BoundaryType getBoundary() const
    { 
      return nr; 
    }
    
  private:
    BoundaryType nr;
    Mesh* mesh;
  };
  
  
  // generator function for \ref BoundaryWrapper
  inline auto boundary(BoundaryType nr, Mesh* mesh = NULL)
    -> BoundaryWrapper
  {
    return {nr, mesh};
  }
  
  
} // end namespace AMDiS
