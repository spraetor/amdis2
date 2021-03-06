#include "Mesh.hpp"

#include <algorithm>
#include <ctime>
#include <map>
#include <set>
#include <stack>

#include "AdaptInstationary.hpp"
#include "AdaptStationary.hpp"
#include "CoarseningManager.hpp"
#include "DOFIterator.hpp"
#include "DOFVector.hpp"
#include "ElInfoStack.hpp"
#include "ElementData.hpp"
#include "ElementDofIterator.hpp"
#include "FiniteElemSpace.hpp"
#include "FixVec.hpp"
#include "Initfile.hpp"
#include "Lagrange.hpp"
#include "MacroElement.hpp"
#include "Projection.hpp"
#include "Traverse.hpp"
#include "VertexVector.hpp"
#include "io/MacroInfo.hpp"
#include "io/MacroReader.hpp"
#include "io/MacroWriter.hpp"
#include "io/Reader.hpp"

using namespace std;

namespace AMDiS
{

  //**************************************************************************
  //  flags, which information should be present in the elInfo structure
  //**************************************************************************

  const Flag Mesh::FILL_NOTHING    = 0X00L;
  const Flag Mesh::FILL_COORDS     = 0X01L;
  const Flag Mesh::FILL_BOUND      = 0X02L;
  const Flag Mesh::FILL_NEIGH      = 0X04L;
  const Flag Mesh::FILL_OPP_COORDS = 0X08L;
  const Flag Mesh::FILL_ORIENTATION= 0X10L;
  const Flag Mesh::FILL_DET        = 0X20L;
  const Flag Mesh::FILL_GRD_LAMBDA = 0X40L;
  const Flag Mesh::FILL_ADD_ALL    = 0X80L;


  const Flag Mesh::FILL_ANY_1D = (0X01L|0X02L|0X04L|0X08L|0x20L|0X40L|0X80L);
  const Flag Mesh::FILL_ANY_2D = (0X01L|0X02L|0X04L|0X08L|0x20L|0X40L|0X80L);
  const Flag Mesh::FILL_ANY_3D = (0X01L|0X02L|0X04L|0X08L|0X10L|0x20L|0X40L|0X80L);

  //**************************************************************************
  //  flags for Mesh traversal
  //**************************************************************************

  const Flag Mesh::CALL_EVERY_EL_PREORDER  = 0X0100L;
  const Flag Mesh::CALL_EVERY_EL_INORDER   = 0X0200L;
  const Flag Mesh::CALL_EVERY_EL_POSTORDER = 0X0400L;
  const Flag Mesh::CALL_LEAF_EL            = 0X0800L;
  const Flag Mesh::CALL_LEAF_EL_LEVEL      = 0X1000L;
  const Flag Mesh::CALL_EL_LEVEL           = 0X2000L;
  const Flag Mesh::CALL_MG_LEVEL           = 0X4000L;
  const Flag Mesh::CALL_REVERSE_MODE       = 0X8000L;


  vector<DegreeOfFreedom> Mesh::dof_used;
  const int Mesh::MAX_DOF = 100;
  map<pair<DegreeOfFreedom, int>, DegreeOfFreedom*> Mesh::serializedDOFs;
  std::set<string> Mesh::refinedMeshNames;

  Mesh::Mesh(string aName, int dimension)
    : name(aName),
      dim(dimension),
      nVertices(0),
      nEdges(0),
      nLeaves(0),
      nElements(0),
      parametric(NULL),
      preserveCoarseDOFs(false),
      nDofEl(0),
      nDof(dimension, 0),
      nNodeEl(0),
      node(dimension, 0),
      elementPrototype(NULL),
      elementDataPrototype(NULL),
      elementIndex(-1),
      initialized(false),
      macroFileInfo(NULL),
      changeIndex(0),
      final_lambda(dimension, 0.0)
  {
    FUNCNAME("Mesh::Mesh()");

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    nParallelPreRefinements = 0;
#endif

    // set default element prototype
    switch(dim)
    {
    case 1:
      elementPrototype = new Line(this);
      break;
    case 2:
      elementPrototype = new Triangle(this);
      break;
    case 3:
      elementPrototype = new Tetrahedron(this);
      break;
    default:
      ERROR_EXIT("invalid dimension\n");
    }

    elementPrototype->setIndex({}, -1);

    elementIndex = 0;
  }


  Mesh::~Mesh()
  {
    deleteMeshStructure();

    if (macroFileInfo != NULL)
      delete macroFileInfo;
    if (elementPrototype)
      delete elementPrototype;
    if (elementDataPrototype)
      delete elementDataPrototype;
  }


  Mesh& Mesh::operator=(Mesh const& m)
  {
    FUNCNAME("Mesh::operator=()");

    if (this == &m)
      return *this;

    TEST_EXIT(dim == m.dim)("operator= works only on meshes with equal dim!\n");

    name = m.name;
    nVertices = m.nVertices;
    nEdges = m.nEdges;
    nLeaves = m.nLeaves;
    nElements = m.nElements;
    nFaces = m.nFaces;
    maxEdgeNeigh = m.maxEdgeNeigh;
    diam = m.diam;
    parametric = NULL;

    preserveCoarseDOFs = m.preserveCoarseDOFs;
    nDofEl = m.nDofEl;
    nDof = m.nDof;
    nNodeEl = m.nNodeEl;
    node = m.node;
    elementIndex = m.elementIndex;
    initialized = m.initialized;
    final_lambda = m.final_lambda;

    /* ====================== Create new DOFAdmins ================== */
    admin.resize(m.admin.size());
    for (size_t i = 0; i < admin.size(); ++i)
    {
      admin[i] = new DOFAdmin(this);
      *(admin[i]) = *(m.admin[i]);
      admin[i]->setMesh(this);
    }


    /* ====================== Copy macro elements =================== */

    // mapIndex[i] is the index of the MacroElement element in the vector
    // macroElements, for which holds: element->getIndex() = i
    map<int, int> mapIndex;

    // We use this map for coping the DOFs of the Elements within the
    // MacroElements objects.
    Mesh::serializedDOFs.clear();

    int insertCounter = 0;

    macroElements.clear();

    // Go through all MacroElements of mesh m, and create for every a new
    // MacroElement in this mesh.
    for (MacroElement* macroElem : m.macroElements)
    {
      // Create new MacroElement.
      MacroElement* el = new MacroElement(dim);

      // Use copy operator to copy all the data to the new MacroElement.
      *el = *macroElem;

      // Make a copy of the Element data, together with all DOFs
      el->setElement(macroElem->getElement()->cloneWithDOFs(Mesh::serializedDOFs));

      // Insert the new MacroElement in the vector of all MacroElements.
      macroElements.push_back(el);

      // Update the index map.
      mapIndex.insert(pair<int, int>(el->getIndex(), insertCounter++));
    }

    // Now we have to go through all the new MacroElements, and update the neighbour
    // connections.
    insertCounter = 0;
    for (MacroElement* macroElem : m.macroElements)
    {
      // Go through all neighbours.
      for (int i = 0; i < dim; i++)
      {
        // 1. Get index of the old MacroElement for its i-th neighbour.
        // 2. Because the index in the new MacroElement is the same, search
        //    for the vector index the corresponding element is stored in.
        // 3. Get this element from macroElements, and set it as the i-th
        //    neighbour for the current element.
        if (macroElem->getNeighbour(i)!=NULL)
        {
          macroElements[insertCounter]->
          setNeighbour(i, macroElements[mapIndex[macroElem->getNeighbour(i)->getIndex()]]);
        }
      }
      insertCounter++;
    }

    // Cleanup
    Mesh::serializedDOFs.clear();

    /* ================== Things will be done when required ============ */

    TEST_EXIT(elementDataPrototype == NULL)("TODO\n");
    TEST_EXIT(m.parametric == NULL)("TODO\n");
    TEST_EXIT(periodicAssociations.size() == 0)("TODO\n");

    return *this;
  }


  void Mesh::updateNumberOfLeaves()
  {
    nLeaves = 0;

    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(this, -1, Mesh::CALL_LEAF_EL);
    while (elInfo)
    {
      nLeaves++;
      elInfo = stack.traverseNext(elInfo);
    }
  }


  void Mesh::addMacroElement(MacroElement* me)
  {
    macroElements.push_back(me);
    me->setIndex(int( macroElements.size() ));
  }
  

  void Mesh::removeAllMacroElements()
  {
    // Delete all the dofs
    Element::clearDeletedDofs({});
    for (MacroElement const* macroEl : macroElements)
      macroEl->getElement()->deleteElementDOFs();

    Element::clearDeletedDofs({});

    // Set all neighbors null
    for (MacroElement* macroEl : macroElements)
    {
      for (int i = 0; i < getGeo(NEIGH); i++)
        macroEl->setNeighbour(i, NULL);

      Element* mel = macroEl->getElement();
      // Delete element hierarchie
      if (!(mel->isLeaf()))
      {
        delete mel->getChild(0);
        delete mel->getChild(1);

        mel->setFirstChild(NULL);
        mel->setSecondChild(NULL);

        mel->setElementData(elementDataPrototype->clone());
      }
    }

    macroElements.clear();
    nLeaves = 0;
    nElements = 0;
    nVertices = 0;

    for (size_t i = 0; i < admin.size(); i++)
    {
      TEST_EXIT_DBG(admin[i]->getUsedSize() == admin[i]->getHoleCount())
      ("All macro elements has been removed. But not all dofs are cleaned. (UsedSize = %d, HoleCount = %d)\n",
       admin[i]->getUsedSize(), admin[i]->getHoleCount());

      admin[i]->reset();
    }
  }
  

  void Mesh::removeMacroElements(std::set<MacroElement*>& delMacros,
                                 vector<const FiniteElemSpace*>& feSpaces)
  {
    FUNCNAME("Mesh::removeMacroElement()");

    using DofElMap  = map<const DegreeOfFreedom*, std::set<MacroElement*>>;
    using DofPosMap = map<const DegreeOfFreedom*, GeoIndex>;

    TEST_EXIT(!feSpaces.empty())("FeSpace-Vector is empty. This should not happen!\n");

    // Search for the FE space with the highest degree of polynomials. Using this
    // FE space ensures that deleting DOFs defined on it, also DOFs of lower
    // order FE spaces will be deleted correctly.
    const FiniteElemSpace* feSpace = FiniteElemSpace::getHighest(feSpaces);


    // === Determine to all DOFs in mesh the macro elements where the DOF  ===
    // === is part of.                                                     ===

    // Map that stores for each DOF pointer (which may have a list of DOFs)
    // all macro element indices that own this DOF.
    DofElMap dofsOwner;
    DofPosMap dofsPosIndex;

    ElementDofIterator elDofIter(feSpace);
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(this, -1, Mesh::CALL_LEAF_EL);
    while (elInfo)
    {
      elDofIter.reset(elInfo->getElement());
      do
      {
        dofsOwner[elDofIter.getBaseDof()].insert(elInfo->getMacroElement());
        dofsPosIndex[elDofIter.getBaseDof()] = elDofIter.getPosIndex();
      }
      while (elDofIter.nextStrict());

      elInfo = stack.traverseNext(elInfo);
    }


    // === Remove macro elements from mesh macro element list. ===

    // Removing arbitrary elements from an deque is very slow. Therefore, we
    // create a new deque with all macro elements that should not be deleted. The
    // macro element deque is than replaced by the new created one.

    deque<MacroElement*> newMacroElements;
    for (MacroElement* macroEl :  macroElements)
    {
      // If the current mesh macro element should not be deleted, i.e., it is not
      // a member of the list of macro elements to be deleted, is is inserted to
      // the new macro element list.
      if (delMacros.find(macroEl) == delMacros.end())
        newMacroElements.push_back(macroEl);
    }

    // And replace the macro element list with the new one.
    macroElements.clear();
    macroElements = newMacroElements;


    // === For all macro elements to be deleted, delete them also to be       ===
    // === neighbours of some other macro elements. Furtheremore, delete the  ===
    // === whole element hierarchie structure of the macro element.           ===

    for (MacroElement* macroEl : delMacros)
    {
      // Go through all neighbours of the macro element and remove this macro
      // element to be neighbour of some other macro element.
      for (int i = 0; i < getGeo(NEIGH); i++)
        if (macroEl->getNeighbour(i))
          for (int j = 0; j < getGeo(NEIGH); j++)
            if (macroEl->getNeighbour(i)->getNeighbour(j) == macroEl)
              macroEl->getNeighbour(i)->setNeighbour(j, NULL);


      Element* mel = macroEl->getElement();
      // Delete element hierarchie
      if (!(mel->isLeaf()))
      {
        delete mel->getChild(0);
        delete mel->getChild(1);

        mel->setFirstChild(NULL);
        mel->setSecondChild(NULL);

        mel->setElementData(elementDataPrototype->clone());
      }

      mel->delDofPtr();
    }


    // === Check now all the DOFs that have no owner anymore and therefore  ===
    // === have to be removed.                                              ===

    for (auto dofsIt : dofsOwner)
    {
      bool deleteDof = true;
      for (MacroElement* macroEl : dofsIt.second)
      {
        if (delMacros.find(macroEl) == delMacros.end())
        {
          deleteDof = false;
          break;
        }
      }

      if (deleteDof)
        freeDof(const_cast<DegreeOfFreedom*>(dofsIt.first),
                dofsPosIndex[dofsIt.first]);
    }


    // === Update number of elements, vertices, etc. ===

    nLeaves = 0;
    nElements = 0;
    nVertices = 0;

    if (!macroElements.empty())
    {
      std::set<const DegreeOfFreedom*> allVertices;

      elInfo = stack.traverseFirst(this, -1, Mesh::CALL_EVERY_EL_PREORDER);
      while (elInfo)
      {
        nElements++;

        if (elInfo->getElement()->isLeaf())
        {
          nLeaves++;

          for (int i = 0; i < getGeo(VERTEX); i++)
            allVertices.insert(elInfo->getElement()->getDof(i));
        }

        elInfo = stack.traverseNext(elInfo);
      }

      nVertices = int( allVertices.size() );
    }
    else
    {

      for (size_t i = 0; i < admin.size(); i++)
      {
        TEST_EXIT_DBG(admin[i]->getUsedSize() == admin[i]->getHoleCount())
        ("All macro elements has been removed. But not all dofs are cleaned. (UsedSize = %d, HoleCount = %d)\n",
         admin[i]->getUsedSize(), admin[i]->getHoleCount());

        admin[i]->reset();
      }
    }
    // === Note: Although the macro elements are removed from the mesh,   ===
    // === they are not deleted from memory. The macro elements are still ===
    // === stored in macroInfo structure. They are needed, if the mesh is ===
    // === redistributed between the ranks.                               ===
  }


  void Mesh::addDOFAdmin(DOFAdmin* localAdmin)
  {
    FUNCNAME("Mesh::addDOFAdmin()");

    localAdmin->setMesh(this);

    TEST_EXIT(find(admin.begin(), admin.end(), localAdmin) == admin.end())
    ("admin %s is already associated to mesh %s\n",
     localAdmin->getName().c_str(), this->getName().c_str());

    admin.push_back(localAdmin);

    nDofEl = 0;

    localAdmin->setNumberOfPreDofs(VERTEX, nDof[VERTEX]);
    nDof[VERTEX] += localAdmin->getNumberOfDofs(VERTEX);
    nDofEl += getGeo(VERTEX) * nDof[VERTEX];

    if (dim > 1)
    {
      localAdmin->setNumberOfPreDofs(EDGE, nDof[EDGE]);
      nDof[EDGE] += localAdmin->getNumberOfDofs(EDGE);
      nDofEl += getGeo(EDGE) * nDof[EDGE];
    }

    localAdmin->setNumberOfPreDofs(CENTER, nDof[CENTER]);
    nDof[CENTER] += localAdmin->getNumberOfDofs(CENTER);
    nDofEl += nDof[CENTER];

    TEST_EXIT_DBG(nDof[VERTEX] > 0)("no vertex dofs\n");

    node[VERTEX] = 0;
    nNodeEl = getGeo(VERTEX);

    if (dim > 1)
    {
      node[EDGE] = nNodeEl;
      if (nDof[EDGE] > 0)
        nNodeEl += getGeo(EDGE);
    }

    if (dim == 3)
    {
      localAdmin->setNumberOfPreDofs(FACE, nDof[FACE]);
      nDof[FACE] += localAdmin->getNumberOfDofs(FACE);
      nDofEl += getGeo(FACE) * nDof[FACE];
      node[FACE] = nNodeEl;
      if (nDof[FACE] > 0)
        nNodeEl += getGeo(FACE);
    }

    node[CENTER] = nNodeEl;
    if (nDof[CENTER] > 0)
      nNodeEl += 1;
  }


  void Mesh::dofCompress()
  {
    FUNCNAME_DBG("Mesh::dofCompress()");

    for (size_t iadmin = 0; iadmin < admin.size(); ++iadmin)
    {
      DOFAdmin* compressAdmin = admin[iadmin];

      TEST_EXIT_DBG(compressAdmin)("no admin[%d] in mesh\n", iadmin);

      int size = compressAdmin->getSize();
      if (size < 1 ||
          compressAdmin->getUsedDofs() < 1 ||
          compressAdmin->getHoleCount() < 1)
        continue;

      vector<DegreeOfFreedom> newDofIndex(size);
      compressAdmin->compress(newDofIndex);

      Flag fill_flag = (preserveCoarseDOFs ?
                        Mesh::CALL_EVERY_EL_PREORDER | Mesh::FILL_NOTHING :
                        Mesh::CALL_LEAF_EL | Mesh::FILL_NOTHING);
      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(this, -1, fill_flag);
      while (elInfo)
      {
        elInfo->getElement()->newDofFct1({}, compressAdmin, newDofIndex);
        elInfo = stack.traverseNext(elInfo);
      }

      elInfo = stack.traverseFirst(this, -1, fill_flag);
      while (elInfo)
      {
        elInfo->getElement()->newDofFct2({}, compressAdmin);
        elInfo = stack.traverseNext(elInfo);
      }
    }
  }


  DegreeOfFreedom* Mesh::getDof(GeoIndex position)
  {
    FUNCNAME_DBG("Mesh::getDof()");

    TEST_EXIT_DBG(position >= CENTER && position <= FACE)
    ("unknown position %d\n", position);

    int ndof = getNumberOfDofs(position);
    if (ndof <= 0)
      return NULL;

    DegreeOfFreedom* dof = new DegreeOfFreedom[ndof];

    for (int i = 0; i < getNumberOfDOFAdmin(); i++)
    {
      const DOFAdmin* localAdmin = &getDofAdmin(i);
      TEST_EXIT_DBG(localAdmin)("no admin[%d]\n", i);

      int n  = localAdmin->getNumberOfDofs(position);
      int n0 = localAdmin->getNumberOfPreDofs(position);

      TEST_EXIT_DBG(n + n0 <= ndof)
      ("n = %d, n0 = %d too large: ndof = %d\n", n, n0, ndof);

      for (int j = 0; j < n; j++)
        dof[n0 + j] = const_cast<DOFAdmin*>(localAdmin)->getDOFIndex(MeshAccessor());
    }

    return dof;
  }


  DegreeOfFreedom** Mesh::createDofPtrs()
  {
    if (nNodeEl <= 0)
      return NULL;

    DegreeOfFreedom** ptrs = new DegreeOfFreedom*[nNodeEl];
    for (int i = 0; i < nNodeEl; i++)
      ptrs[i] = NULL;

    return ptrs;
  }


  void Mesh::freeDofPtrs(DegreeOfFreedom** ptrs)
  {
    FUNCNAME_DBG("Mesh::freeDofPtrs()");

    TEST_EXIT_DBG(ptrs)("ptrs is NULL!\n");

    if (nNodeEl <= 0)
      return;

    delete [] ptrs;
    ptrs = NULL;
  }


  const DOFAdmin* Mesh::createDOFAdmin(string lname, DimVec<int> lnDof)
  {
    DOFAdmin* localAdmin = new DOFAdmin(this, lname);

    for (int i = 0; i < dim + 1; i++)
      localAdmin->setNumberOfDofs(i, lnDof[i]);

    addDOFAdmin(localAdmin);

    return localAdmin;
  }


  const DOFAdmin* Mesh::getVertexAdmin() const
  {
    const DOFAdmin* localAdmin = NULL;

    for (unsigned int i = 0; i < admin.size(); i++)
    {
      if (admin[i]->getNumberOfDofs(VERTEX))
      {
        if (!localAdmin)
          localAdmin = admin[i];
        else if (admin[i]->getSize() < localAdmin->getSize())
          localAdmin = admin[i];
      }
    }

    return localAdmin;
  }


  void Mesh::freeDof(DegreeOfFreedom* dof, GeoIndex position)
  {
    FUNCNAME_DBG("Mesh::freeDof()");

    TEST_EXIT_DBG(position >= CENTER && position <= FACE)
    ("unknown position %d\n", position);

    if (nDof[position])
    {
      TEST_EXIT_DBG(dof != NULL)("dof = NULL, but ndof = %d\n", nDof[position]);
    }
    else
    {
      TEST_EXIT_DBG(dof == NULL)("dof != NULL, but ndof = 0\n");
    }

    TEST_EXIT_DBG(nDof[position] <= MAX_DOF)
    ("ndof too big: ndof = %d, MAX_DOF = %d\n", nDof[position], MAX_DOF);

    for (unsigned int i = 0; i < admin.size(); i++)
    {
      int n = admin[i]->getNumberOfDofs(position);
      int n0 = admin[i]->getNumberOfPreDofs(position);

      TEST_EXIT_DBG(n + n0 <= nDof[position])
      ("n = %d, n0 = %d too large: ndof = %d\n", n, n0, nDof[position]);

      for (int j = 0; j < n; j++)
        admin[i]->freeDofIndex(MeshAccessor(), dof[n0 + j]);
    }
    delete [] dof;
  }


  void Mesh::freeElement(Element* el)
  {
    freeDofPtrs(const_cast<DegreeOfFreedom**>(el->getDof()));
    delete el;
  }


  Element* Mesh::createNewElement(Element* parent)
  {
    FUNCNAME_DBG("Mesh::createNewElement()");

    TEST_EXIT_DBG(elementPrototype)("no element prototype\n");

    Element* el = parent ? parent->clone() : elementPrototype->clone();

    if (!parent && elementDataPrototype)
      el->setElementData(elementDataPrototype->clone());
    else
      el->setElementData(NULL); // must be done in ElementData::refineElementData()

    return el;
  }


  ElInfo* Mesh::createNewElInfo()
  {
    FUNCNAME("Mesh::createNewElInfo()");

    switch (dim)
    {
    case 1:
      return new ElInfo1d(this);
      break;
    case 2:
      return new ElInfo2d(this);
      break;
    case 3:
      return new ElInfo3d(this);
      break;
    default:
      ERROR_EXIT("invalid dim [%d]\n",dim);
      return NULL;
    }
  }


  bool Mesh::findElInfoAtPoint(const WorldVector<double>& xy,
                               ElInfo* el_info,
                               DimVec<double>& bary,
                               const MacroElement* start_mel,
                               const WorldVector<double>* xy0,
                               double* sp)
  {
    static const MacroElement* mel = NULL;
    DimVec<double> lambda(dim);
    ElInfo* mel_info = createNewElInfo();

    if (start_mel != NULL)
      mel = start_mel;
    else if (mel == NULL || mel->getElement()->getMesh() != this)
      mel = *(macroElements.begin());

    mel_info->setFillFlag(Mesh::FILL_COORDS);
    g_xy = &xy;
    g_xy0 = xy0;
    g_sp = sp;

    mel_info->fillMacroInfo(mel);

    // We have the care about not to visite a macro element twice. In this case the
    // function would end up in an infinite loop. If a macro element is visited a
    // second time, what can happen with periodic boundary conditions, the point is
    // not within the mesh!
    std::set<int> macrosVisited;
    std::stack<MacroElement*> active;

    //     macrosVisited.insert(mel->getIndex());

    int k;
    while ((k = mel_info->worldToCoord(xy, lambda)) >= 0)
    {
      macrosVisited.insert(mel->getIndex());

      if (mel->getNeighbour(k) && !macrosVisited.count(mel->getNeighbour(k)->getIndex()))
      {
        // look for next macro-element in the direction of the coordinates xy
        mel = mel->getNeighbour(k);
        mel_info->fillMacroInfo(mel);
        continue;
      }
      else
      {
        // consider all neighbors of the current macro-element to visit next
        for (int i = 0; i < dim + 1; ++i)
        {
          if (i !=  k && mel->getNeighbour(i) && !macrosVisited.count(mel->getNeighbour(i)->getIndex()))
            active.push(mel->getNeighbour(i));
        }
        // if all neighbors are visited already
        if (active.empty())
        {
          if (macrosVisited.size() == static_cast<size_t>(getNumberOfMacros()))
          {
            // if all macro-elements are visited -> no element found!
            delete mel_info;
            return false;
          }
          else
          {
            // go to an arbitrary macro-element to continue the search
            deque<MacroElement*>::iterator it;
            bool found = false;
            for (it = firstMacroElement(); it != endOfMacroElements(); it++)
            {
              if (!macrosVisited.count((*it)->getIndex()))
              {
                active.push(*it);
                found = true;
              }
            }
            if (!found)
            {
              delete mel_info;
              return false;
            }
          }
        }
        mel = active.top();
        active.pop();

        mel_info->fillMacroInfo(mel);
      }
    }

    /* now, descend in tree to find leaf element at point */
    bool inside = findElementAtPointRecursive(mel_info, lambda, k, el_info);
    for (int i = 0; i <= dim; i++)
      bary[i] = final_lambda[i];

    delete mel_info;

    return inside;
  }


  bool Mesh::findElementAtPoint(const WorldVector<double>&  xy,
                                Element** elp,
                                DimVec<double>& bary,
                                const MacroElement* start_mel,
                                const WorldVector<double>* xy0,
                                double* sp)
  {
    ElInfo* el_info = createNewElInfo();
    bool inside = findElInfoAtPoint(xy, el_info, bary, start_mel, xy0, sp);

    *elp = el_info->getElement();

    delete el_info;

    return inside;
  }


  bool Mesh::findElementAtPointRecursive(ElInfo* el_info,
                                         const DimVec<double>& lambda,
                                         int outside,
                                         ElInfo* final_el_info)
  {
    FUNCNAME("Mesh::findElementAtPointRecursive()");

    Element* el = el_info->getElement();
    DimVec<double> c_lambda(dim);
    bool inside;
    int ichild, c_outside;

    if (el->isLeaf())
    {
      *final_el_info = *el_info;
      if (outside < 0)
      {
        for (int i = 0; i <= dim; i++)
          final_lambda[i] = lambda[i];

        return true;
      }
      else       /* outside */
      {
        if (g_xy0)   /* find boundary point of [xy0, xy] */
        {
          el_info->worldToCoord(*(g_xy0), c_lambda);
          double s = lambda[outside] / (lambda[outside] - c_lambda[outside]);
          for (int i = 0; i <= dim; i++)
            final_lambda[i] = s * c_lambda[i] + (1.0 - s) * lambda[i];

          if (g_sp)
            *(g_sp) = s;

          if (dim == 3)
            MSG("Outside finest level on el %d: s = %.3e\n", el->getIndex(), s);

          return false;  /* ??? */
        }
        else
        {
          return false;
        }
      }
    }

    ElInfo* c_el_info = createNewElInfo();

    if (dim == 1)
    {
      if (lambda[0] >= lambda[1])
      {
        c_el_info->fillElInfo(0, el_info);
        if (outside >= 0)
        {
          outside = el_info->worldToCoord(*(g_xy), c_lambda);
          TEST_EXIT(outside == 0)("point outside domain\n");
        }
        else
        {
          c_lambda[0] = lambda[0] - lambda[1];
          c_lambda[1] = 2.0 * lambda[1];
        }
      }
      else
      {
        c_el_info->fillElInfo(1, el_info);
        if (outside >= 0)
        {
          outside = el_info->worldToCoord(*(g_xy), c_lambda);
          TEST_EXIT(outside == 0)("point outside domain\n");
        }
        else
        {
          c_lambda[1] = lambda[1] - lambda[0];
          c_lambda[0] = 2.0 * lambda[0];
        }
      }
    } /* DIM == 1 */

    if (dim == 2)
    {
      if (lambda[0] >= lambda[1])
      {
        c_el_info->fillElInfo(0, el_info);
        if (el->isNewCoordSet())
        {
          outside = c_el_info->worldToCoord(*(g_xy), c_lambda);
          TEST_EXIT(outside == 0)("outside curved boundary child 0\n");
        }
        else
        {
          c_lambda[0] = lambda[2];
          c_lambda[1] = lambda[0] - lambda[1];
          c_lambda[2] = 2.0 * lambda[1];
        }
      }
      else
      {
        c_el_info->fillElInfo(1, el_info);
        if (el->isNewCoordSet())
        {
          outside = c_el_info->worldToCoord(*(g_xy), c_lambda);
          TEST_EXIT(outside == 0)("outside curved boundary child 1\n");
        }
        else
        {
          c_lambda[0] = lambda[1] - lambda[0];
          c_lambda[1] = lambda[2];
          c_lambda[2] = 2.0 * lambda[0];
        }
      }
    } /* DIM == 2 */

    if (dim == 3)
    {
      if (el->isNewCoordSet())
      {
        if (lambda[0] >= lambda[1])
          ichild = 0;
        else
          ichild = 1;
        c_el_info->fillElInfo(ichild, el_info);
        c_outside = c_el_info->worldToCoord(*(g_xy), c_lambda);

        if (c_outside>=0)    /* test is other child is better... */
        {
          DimVec<double> c_lambda2(dim);
          ElInfo* c_el_info2 = createNewElInfo();

          c_el_info2->fillElInfo(1-ichild, el_info);
          int c_outside2 = c_el_info2->worldToCoord(*(g_xy), c_lambda2);

          MSG("new_coord CHILD %d: outside=%d, lambda=(%.2f %.2f %.2f %.2f)\n",
              ichild, c_outside, c_lambda[0], c_lambda[1], c_lambda[2], c_lambda[3]);
          MSG("new_coord CHILD %d: outside=%d, lambda=(%.2f %.2f %.2f %.2f)\n",
              1 - ichild, c_outside2, c_lambda2[0], c_lambda2[1], c_lambda2[2],
              c_lambda2[3]);

          if ((c_outside2 < 0) || (c_lambda2[c_outside2] > c_lambda[c_outside]))
          {
            for (int i = 0; i <= dim; i++)
              c_lambda[i] = c_lambda2[i];
            c_outside = c_outside2;
            *c_el_info = *c_el_info2;
            ichild = 1 - ichild;
          }
          delete c_el_info2;
        }
        outside = c_outside;
      }
      else      /* no new_coord */
      {
        if (lambda[0] >= lambda[1])
        {
          c_el_info->fillElInfo(0, el_info);
          c_lambda[0] = lambda[0] - lambda[1];
          c_lambda[1] =
            lambda[Tetrahedron::childVertex[(dynamic_cast<ElInfo3d*>(el_info))->
                                            getType()][0][1]];
          c_lambda[2] =
            lambda[Tetrahedron::childVertex[(dynamic_cast<ElInfo3d*>(el_info))->
                                            getType()][0][2]];
          c_lambda[3] = 2.0 * lambda[1];
        }
        else
        {
          c_el_info->fillElInfo(1, el_info);
          c_lambda[0] = lambda[1] - lambda[0];
          c_lambda[1] =
            lambda[Tetrahedron::childVertex[(dynamic_cast<ElInfo3d*>(el_info))->
                                            getType()][1][1]];
          c_lambda[2] =
            lambda[Tetrahedron::childVertex[(dynamic_cast<ElInfo3d*>(el_info))->
                                            getType()][1][2]];
          c_lambda[3] = 2.0 * lambda[0];
        }
      }
    }  /* DIM == 3 */

    inside = findElementAtPointRecursive(c_el_info, c_lambda, outside, final_el_info);
    delete c_el_info;

    return inside;
  }


  bool Mesh::getDofIndexCoords(DegreeOfFreedom dof,
                               const FiniteElemSpace* feSpace,
                               WorldVector<double>& coords)
  {
    DimVec<double>* baryCoords;
    bool found = false;
    TraverseStack stack;
    vector<DegreeOfFreedom> dofVec(feSpace->getBasisFcts()->getNumber());

    ElInfo* elInfo = stack.traverseFirst(this, -1,
                                         Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
    while (elInfo)
    {
      feSpace->getBasisFcts()->getLocalIndices(elInfo->getElement(),
          feSpace->getAdmin(),
          dofVec);
      for (int i = 0; i < feSpace->getBasisFcts()->getNumber(); i++)
      {
        if (dofVec[i] == dof)
        {
          baryCoords = feSpace->getBasisFcts()->getCoords(i);
          elInfo->coordToWorld(*baryCoords, coords);
          found = true;
          break;
        }
      }

      if (found)
        break;

      elInfo = stack.traverseNext(elInfo);
    }

    return found;
  }


  void Mesh::getDofIndexCoords(DOFVector<WorldVector<double>>& coords)
  {
    const FiniteElemSpace* feSpace = coords.getFeSpace();
    const BasisFunction* basFcts = feSpace->getBasisFcts();
    int nBasFcts = basFcts->getNumber();
    vector<DegreeOfFreedom> dofVec(nBasFcts);

    TraverseStack stack;
    ElInfo* elInfo =
      stack.traverseFirst(this, -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
    while (elInfo)
    {
      basFcts->getLocalIndices(elInfo->getElement(), feSpace->getAdmin(), dofVec);
      for (int i = 0; i < nBasFcts; i++)
      {
        DimVec<double>* baryCoords = basFcts->getCoords(i);
        elInfo->coordToWorld(*baryCoords, coords[dofVec[i]]);
      }

      elInfo = stack.traverseNext(elInfo);
    }

  }


  void Mesh::getAllDofs(const FiniteElemSpace* feSpace,
                        std::set<const DegreeOfFreedom*>& allDofs)
  {
    ElementDofIterator elDofIt(feSpace);
    allDofs.clear();

    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(this, -1, Mesh::CALL_LEAF_EL);
    while (elInfo)
    {
      elDofIt.reset(elInfo->getElement());
      do
      {
        allDofs.insert(elDofIt.getDofPtr());
      }
      while(elDofIt.next());

      elInfo = stack.traverseNext(elInfo);
    }
  }


  void Mesh::setDiameter(const WorldVector<double>& w)
  {
    diam = w;
  }


  void Mesh::setDiameter(int i, double w)
  {
    diam[i] = w;
  }


  void Mesh::setBoundingBox(WorldVector<double> const& min_corner, 
                            WorldVector<double> const& max_corner)
  {
    boundingBox = std::make_pair(min_corner, max_corner);
  }



  void Mesh::initialize()
  {
    FUNCNAME("Mesh::initialize()");

    TEST_EXIT(admin.size() > 0)("No DOF admin defined!\n");

    string macroFilename("");
    string valueFilename("");
    string periodicFilename("");
    int check = 1;
    int strategy = 0;
    bool preserveMacroFileInfo = false;

    Parameters::get(name + "->macro file name", macroFilename);
    Parameters::get(name + "->value file name", valueFilename);
    Parameters::get(name + "->periodic file", periodicFilename);
    Parameters::get(name + "->check", check);
    Parameters::get(name + "->preserve coarse dofs", preserveCoarseDOFs);
    Parameters::get(name + "->preserve macroFileInfo", preserveMacroFileInfo);


    Parameters::get("parallel->repartitioning->strategy", strategy);
    if (strategy && !preserveMacroFileInfo)
    {
      MSG("Preserve macroFileInfo.\n");
      preserveMacroFileInfo = true;
    }

    TEST_EXIT(macroFilename.length())
    ("No mesh defined for parameter %s->macro file name !\n", name.c_str());

    // In parallel computations, check if a finer macro mesh is required.
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    checkParallelMacroFile(macroFilename, periodicFilename, check);
#endif

    macroFileInfo =
      io::MacroReader::readMacro(macroFilename, this, periodicFilename, check);

    if (!valueFilename.length() && !preserveMacroFileInfo)
      clearMacroFileInfo();

    initialized = true;

    string arhFilename("");
    Parameters::get(name + "->arh file name", arhFilename);
    if (arhFilename != "")
      io::readFile(arhFilename, this);
  }


#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
  void Mesh::checkParallelMacroFile(string& macroFilename,
                                    string& periodicFilename,
                                    int check)
  {
    FUNCNAME("Mesh::checkParallelMacroFile()");

    // === Create a temporary mesh and load the macro file to it. ===

    Mesh testMesh(name, dim);
    testMesh.setElementDataPrototype(new LeafDataEstimatableVec(new LeafDataCoarsenableVec));
    DOFAdmin* localAdmin = new DOFAdmin(&testMesh, admin[0]->getName());
    localAdmin->setNumberOfDofs(admin[0]->getNumberOfDofs());
    testMesh.addDOFAdmin(localAdmin);

    MacroInfo* testMacroInfo =
      io::MacroReader::readMacro(macroFilename, &testMesh, periodicFilename, check);
    testMacroInfo->clear();
    delete testMacroInfo;


    // === Check the mesh structure. ===

    int nMacroElements = 0;
    int elType = -1;
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(&testMesh, -1, Mesh::CALL_LEAF_EL);
    while (elInfo)
    {
      if (elType == -1)
      {
        elType = elInfo->getType();
      }
      else
      {
        TEST_EXIT(elType == elInfo->getType())
        ("All elements in mesh must have the same element type!\n");
      }

      nMacroElements++;

      elInfo = stack.traverseNext(elInfo);
    }


    // === Calculate the number of global refinements such that the new mesh ===
    // === would fulfill all requirements.                                   ===

    // There should be at least 10 macro Elements per processor, therefore:
    //        nMacroElements * 2^gr >= nProcs * 10
    //    =>  gr = log_2(nProcs * 10 / nMacroElements)
    int minElemNum = 10;
    Parameters::get("parallel->minion number of macro elements per processor", minElemNum);
    double scale = static_cast<double>(minElemNum) * MPI::COMM_WORLD.Get_size() / nMacroElements;
    nParallelPreRefinements =
      static_cast<int>(std::max(0.0, ceil(log(scale) / log(2))));

    if (dim == 3)
    {
      int newElType = (elType + nParallelPreRefinements) % 3;
      switch (newElType)
      {
      case 1:
        if (nParallelPreRefinements > 0)
          nParallelPreRefinements--;
        else
          nParallelPreRefinements = 2;
        break;
      case 2:
        nParallelPreRefinements++;
        break;
      }

      TEST_EXIT((elType + nParallelPreRefinements) % 3 == 0)
      ("This should not happen!\n");
    }


    // === Check if number of pre refinements is set in init file. ===

    int tmp = -1;
    Parameters::get("parallel->pre refine", tmp);
    if (tmp > -1)
    {
      MSG("Calculated %d pre refines to be useful, but %d are set in init file!\n",
          nParallelPreRefinements, tmp);
      nParallelPreRefinements = tmp;
    }

    // === If we do not need to refine the mesh, return back. ===

    if (nParallelPreRefinements == 0)
      return;


    // === If macro weights are explicitly given, we cannot change the mesh. ===

    string macroWeightsFilename = "";
    Parameters::get(name + "->macro weights", macroWeightsFilename);
    if (macroWeightsFilename != "")
    {
      ERROR_EXIT("Should not happen!\n");
    }




    // === Create unique file names. ===

    int filenameRandomNumber = 0;
    if (MPI::COMM_WORLD.Get_rank() == 0)
    {
      srand(time(0));
      filenameRandomNumber = rand() % 1000000;
    }
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Bcast(&filenameRandomNumber, 1, MPI_INT, 0);
    MPI::COMM_WORLD.Barrier();

    stringstream newMacroFilename;
    newMacroFilename << macroFilename << "." << filenameRandomNumber << ".tmp";

    stringstream newPeriodicFilename;
    newPeriodicFilename << periodicFilename << "." << filenameRandomNumber << ".tmp";


    // === Rank 0 creates a new mesh file. ===

    if (MPI::COMM_WORLD.Get_rank() == 0)
    {
      RefinementManager* refManager;
      if (dim == 2)
        refManager = new RefinementManager2d();
      else if (dim == 3)
        refManager = new RefinementManager3d();
      else
      {
        ERROR_EXIT("Parallel model only support dim = 2 or 3.\n");
        return;
      }

      refManager->globalRefine(&testMesh, nParallelPreRefinements);
      delete refManager;

      Lagrange* basFcts = Lagrange::getLagrange(dim, 1);
      FiniteElemSpace* feSpace =
        FiniteElemSpace::provideFeSpace(localAdmin, basFcts, &testMesh, "tmp");

      DataCollector<> dc(feSpace);
      io::MacroWriter::writeMacro(&dc, newMacroFilename.str().c_str());

      if (periodicFilename != "")
        io::MacroWriter::writePeriodicFile(&dc, newPeriodicFilename.str().c_str());
    }


    // === All ranks must wait until rank 0 has created the new macro mesh file. ===

    MPI::COMM_WORLD.Barrier();


    // === We have refined the mesh, so reduce the number of global refinements. ===

    // This is avoid to minus parallel prerefinement multiple times
    if (refinedMeshNames.count(name + "->global refinements") == 0)
    {
      refinedMeshNames.insert(name + "->global refinements");

      int globalRefinements = 0;
      Parameters::get(name + "->global refinements", globalRefinements);

      if (globalRefinements < nParallelPreRefinements)
        globalRefinements = 0;
      else
        globalRefinements -= nParallelPreRefinements;

      Parameters::set(name + "->global refinements", globalRefinements);
    }


    // === Print a note to the screen that another mesh file will be used. ===

    MSG("The macro mesh file \"%s\" was refined %d times and stored to file \"%s\".\n",
        macroFilename.c_str(), nParallelPreRefinements, newMacroFilename.str().c_str());

    macroFilename = newMacroFilename.str();
    if (periodicFilename != "")
      periodicFilename = newPeriodicFilename.str();
  }
#endif


  bool Mesh::associated(DegreeOfFreedom dof1, DegreeOfFreedom dof2)
  {
    map<BoundaryType, VertexVector*>::iterator it;
    map<BoundaryType, VertexVector*>::iterator end = periodicAssociations.end();

    for (it = periodicAssociations.begin(); it != end; ++it)
      if ((*(it->second))[dof1] == dof2)
        return true;

    return false;
  }


  bool Mesh::indirectlyAssociated(DegreeOfFreedom dof1, DegreeOfFreedom dof2)
  {
    vector<DegreeOfFreedom> associatedToDOF1;
    map<BoundaryType, VertexVector*>::iterator it;
    map<BoundaryType, VertexVector*>::iterator end = periodicAssociations.end();
    DegreeOfFreedom dof, assDOF;

    associatedToDOF1.push_back(dof1);
    for (it = periodicAssociations.begin(); it != end; ++it)
    {
      int size = static_cast<int>(associatedToDOF1.size());
      for (int i = 0; i < size; i++)
      {
        dof = associatedToDOF1[i];
        assDOF = (*(it->second))[dof];
        if (assDOF == dof2)
        {
          return true;
        }
        else
        {
          if (assDOF != dof)
            associatedToDOF1.push_back(assDOF);
        }
      }
    }

    return false;
  }


  void Mesh::clearMacroFileInfo()
  {
    macroFileInfo->clear();
    delete macroFileInfo;
    macroFileInfo = NULL;
  }


  int Mesh::calcMemoryUsage()
  {
    int result = sizeof(Mesh);

    result += nDofEl;
    for (int i = 0; i < static_cast<int>(admin.size()); i++)
    {
      result += admin[i]->calcMemoryUsage();
      result += admin[i]->getUsedSize() * sizeof(DegreeOfFreedom);
    }

    return result;
  }


  void Mesh::deleteMeshStructure()
  {
    Element::clearDeletedDofs({});

    for (deque<MacroElement*>::const_iterator it = macroElements.begin();
         it != macroElements.end(); ++it)
    {
      (*it)->getElement()->deleteElementDOFs();
      delete *it;
    }

    Element::clearDeletedDofs({});
  }


  void Mesh::getElementIndexMap(map<int, Element*>& elIndexMap)
  {
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(this, -1, Mesh::CALL_EVERY_EL_PREORDER);
    while (elInfo)
    {
      Element* el = elInfo->getElement();
      elIndexMap[el->getIndex()] = el;
      elInfo = stack.traverseNext(elInfo);
    }
  }

} // end namespace AMDiS
