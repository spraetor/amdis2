#include "Element.hpp"
#include "DOFAdmin.hpp"
#include "Mesh.hpp"
#include "CoarseningManager.hpp"
#include "FixVec.hpp"
#include "ElementRegion_ED.hpp"
#include "MeshStructure.hpp"
#include "BasisFunction.hpp"
#include "FiniteElemSpace.hpp"

namespace AMDiS
{
  std::map<DegreeOfFreedom*, bool> Element::deletedDOFs;

  Element::Element(Mesh* aMesh)
    : dof(NULL)
  {
    mesh = aMesh;
    index = mesh ? mesh->getNextElementIndex() : -1;
    child[0] = NULL;
    child[1] = NULL;
    newCoord = NULL;
    elementData = NULL;
    mark = 0;

    if (mesh)
      createNewDofPtrs();
  }


  // call destructor through Mesh::freeElement !!!
  Element::~Element()
  {
    if (child[0])
      delete child[0];

    if (child[1])
      delete child[1];

    if (newCoord)
      delete newCoord;

    if (elementData)
    {
      elementData->deleteDecorated();
      delete elementData;
    }

    // ;TODO: eventuell nicht löschen!
    if (dof)
    {
      delete [] dof;
      dof = NULL;
    }
  }


  int Element::getRegion() const
  {
    if (!elementData)
      return -1;

    ElementRegion_ED* red =
      dynamic_cast<ElementRegion_ED*>(elementData->getElementData(ELEMENT_REGION));

    if (red)
      return red->getRegion();

    return -1;
  }


  void Element::createNewDofPtrs(bool setDofs)
  {
    FUNCNAME_DBG("Element::setDofPtrs()");

    TEST_EXIT_DBG(mesh)("no mesh!\n");

    dof = mesh->createDofPtrs();

    if (setDofs)
    {
      for (int i = 0; i < mesh->getGeo(VERTEX); i++)
        dof[i] = mesh->getDof(VERTEX);

      if (mesh->getNumberOfDofs(EDGE))
        for (int i = 0; i < mesh->getGeo(EDGE); i++)
          dof[mesh->getNode(EDGE) + i] = mesh->getDof(EDGE);

      if (mesh->getDim() == 3 && mesh->getNumberOfDofs(FACE))
        for (int i = 0; i < mesh->getGeo(FACE); i++)
          dof[mesh->getNode(FACE) + i] = mesh->getDof(FACE);

      if (mesh->getNumberOfDofs(CENTER))
        for (int i = 0; i < mesh->getGeo(CENTER); i++)
          dof[mesh->getNode(CENTER) + i] = mesh->getDof(CENTER);
    }
  }


  bool Element::deleteElementData(int typeID)
  {
    if (elementData)
    {
      if (elementData->isOfType(typeID))
      {
        ElementData* tmp = elementData;
        elementData = elementData->getDecorated();
        delete tmp;
        tmp = NULL;
        return true;
      }
      else
      {
        return elementData->deleteDecorated(typeID);
      }
    }
    return false;
  }


  void Element::deleteElementDOFs()
  {
    int dim = mesh->getDim();
    int j = 0;

    for (int pos = 0; pos <= dim; pos++)
    {
      GeoIndex position = INDEX_OF_DIM(pos, dim);
      int ndof = 0;

      for (int i = 0; i < mesh->getNumberOfDOFAdmin(); i++)
        ndof += mesh->getDofAdmin(i).getNumberOfDofs(position);

      if (ndof > 0)
      {
        for (int i = 0; i < mesh->getGeo(position); i++)
        {
          if (dof[j])
          {
            if (deletedDOFs.count(dof[j]) == 0)
            {
              deletedDOFs[dof[j]] = true;
              // 	      delete [] dof[j];
              mesh->freeDof(dof[j], position);
            }
          }
          j++;
        }
      }
    }

    delete [] dof;
    dof = NULL;

    if (child[0])
      child[0]->deleteElementDOFs();
    if (child[1])
      child[1]->deleteElementDOFs();
  }


  Element* Element::cloneWithDOFs(std::map<std::pair<DegreeOfFreedom, int>, DegreeOfFreedom*>& serializedDOFs)
  {
    Element* el;

    if (isLine())
    {
      el = new Line(NULL);
    }
    else if (isTriangle())
    {
      el = new Triangle(NULL);
    }
    else
    {
      el = new Tetrahedron(NULL);
    }

    el->mesh = mesh;
    el->index = index;
    el->mark = mark;
    if (newCoord)
    {
      WorldVector<double>* nc = new WorldVector<double>();
      *nc = *newCoord;
      el->newCoord = nc;
    }

    /* =========== And here we clone the DOFs =========== */

    el->dof = new DegreeOfFreedom*[mesh->getNumberOfNodes()];

    int dim = mesh->getDim();
    int j = 0;

    for (int pos = 0; pos <= dim; pos++)
    {
      GeoIndex position = INDEX_OF_DIM(pos, dim);
      int ndof = 0;

      for (int i = 0; i < mesh->getNumberOfDOFAdmin(); i++)
        ndof += mesh->getDofAdmin(i).getNumberOfDofs(position);

      if (ndof > 0)
      {
        for (int i = 0; i < mesh->getGeo(position); i++)
        {
          if (dof[j] != NULL)
          {
            std::pair<DegreeOfFreedom, int> idx = std::make_pair(dof[j][0], pos);

            if (serializedDOFs[idx] == NULL)
            {
              el->dof[j] = new DegreeOfFreedom[ndof];
              for (int k = 0; k < ndof; k++)
                el->dof[j][k] = dof[j][k];

              serializedDOFs[idx] = el->dof[j];
            }
            else
            {
              el->dof[j] = serializedDOFs[idx];
            }
          }
          else
          {
            el->dof[j] = NULL;
          }
          j++;
        }
      }
    }

    /* =========== And clone the children ============= */

    if (child[0])
      el->child[0] = child[0]->cloneWithDOFs(serializedDOFs);
    if (child[1])
      el->child[1] = child[1]->cloneWithDOFs(serializedDOFs);

    return el;
  }


  /****************************************************************************/
  /*  ATTENTION:                                                              */
  /*  new_dof_fct() destroys new_dof !!!!!!!!!!                               */
  /*  should be used only at the end of dof_compress()!!!!!                   */
  /****************************************************************************/

  void Element::newDofFct1(MeshAccessor const&, const DOFAdmin* admin, std::vector<DegreeOfFreedom>& newDofIndex)
  {
    int n0, nd, nd0;

    if ((nd = admin->getNumberOfDofs(VERTEX)))
    {
      int vertices = mesh->getGeo(VERTEX);
      nd0 = admin->getNumberOfPreDofs(VERTEX);
      n0 = admin->getMesh()->getNode(VERTEX);
      for (int i = 0; i < vertices; i++)
        changeDofs1(admin, newDofIndex, n0, nd0, nd, i);
    }

    if (mesh->getDim() > 1)
    {
      if ((nd = admin->getNumberOfDofs(EDGE)))
      {
        int edges = mesh->getGeo(EDGE);
        nd0 = admin->getNumberOfPreDofs(EDGE);
        n0 = admin->getMesh()->getNode(EDGE);
        for (int i = 0; i < edges; i++)
          changeDofs1(admin, newDofIndex, n0, nd0, nd, i);
      }
    }

    if (mesh->getDim() == 3)
    {
      if ((nd = admin->getNumberOfDofs(FACE)))
      {
        int faces = mesh->getGeo(FACE);
        nd0 = admin->getNumberOfPreDofs(FACE);
        n0 = admin->getMesh()->getNode(FACE);
        for (int i = 0; i < faces; i++)
          changeDofs1(admin, newDofIndex, n0, nd0, nd, i);
      }
    }

    if ((nd = admin->getNumberOfDofs(CENTER)))
    {
      nd0 = admin->getNumberOfPreDofs(CENTER);
      n0 = admin->getMesh()->getNode(CENTER);
      changeDofs1(admin, newDofIndex, n0, nd0, nd, 0);
    }
  }


  void Element::newDofFct2(MeshAccessor const&, const DOFAdmin* admin)
  {
    int n0, nd0;

    int nd = admin->getNumberOfDofs(VERTEX);
    if (nd)
    {
      int vertices = mesh->getGeo(VERTEX);
      nd0 = admin->getNumberOfPreDofs(VERTEX);
      n0 = admin->getMesh()->getNode(VERTEX);
      for (int i = 0; i < vertices; i++)
        changeDofs2(n0, nd0, nd, i);
    }

    if (mesh->getDim() > 1)
    {
      nd = admin->getNumberOfDofs(EDGE);
      if (nd)
      {
        int edges = mesh->getGeo(EDGE);
        nd0 = admin->getNumberOfPreDofs(EDGE);
        n0 = admin->getMesh()->getNode(EDGE);
        for (int i = 0; i < edges; i++)
          changeDofs2(n0, nd0, nd, i);
      }
    }

    if (mesh->getDim() == 3)
    {
      nd = admin->getNumberOfDofs(FACE);
      if (nd)
      {
        int faces = mesh->getGeo(FACE);
        nd0 = admin->getNumberOfPreDofs(FACE);
        n0 = admin->getMesh()->getNode(FACE);
        for (int i = 0; i < faces; i++)
          changeDofs2(n0, nd0, nd, i);
      }
    }

    nd = admin->getNumberOfDofs(CENTER);
    if (nd)
    {
      nd0 = admin->getNumberOfPreDofs(CENTER);
      n0 = admin->getMesh()->getNode(CENTER);
      // only one center
      changeDofs2(n0, nd0, nd, 0);
    }
  }


  void Element::changeDofs1(const DOFAdmin* /*admin*/, std::vector<DegreeOfFreedom>& newDofIndex,
                            int n0, int nd0, int nd, int pos)
  {
    DegreeOfFreedom* ldof = dof[n0 + pos] + nd0;
    for (int j = 0; j < nd; j++)
    {
      DegreeOfFreedom k = ldof[j];
      if (k >= 0)
        ldof[j] = -newDofIndex[k] - 1;
    }
  }


  void Element::changeDofs2(int n0, int nd0, int nd, int pos)
  {
    DegreeOfFreedom* ldof = dof[n0 + pos] + nd0;
    for (int j = 0; j < nd; j++)
    {
      DegreeOfFreedom k = ldof[j];
      if (k < 0)
        ldof[j] = -k - 1;
    }
  }


  /****************************************************************************/
  /* opp_vertex checks whether the face with vertices dof[0],..,dof[DIM-1] is */
  /* part of mel's boundary. returns the opposite vertex if true, -1 else     */
  /****************************************************************************/

  int Element::oppVertex(FixVec<DegreeOfFreedom*, DIMEN> pdof) const
  {
    FUNCNAME("Element::oppVertex()");

    int nv = 0;
    int ov = 0;
    int vertices = mesh->getGeo(VERTEX);
    int dim = mesh->getDim();

    for (int i = 0; i < vertices; i++)
    {
      if (nv < i - 1)
        return(-1);

      for (int j = 0; j < dim; j++)
      {
        if (dof[i] == pdof[j])
        {
          // i is a common vertex
          ov += i;
          nv++;
          break;
        }
      }

    }

    if (nv != mesh->getDim())
      return(-1);

    // The opposite vertex is 3(6) - (sum of indices of common vertices) in 2D(3D)

    switch(mesh->getDim())
    {
    case 1:
      return ov;
      break;
    case 2:
      return 3 - ov;
      break;
    case 3:
      return 6 - ov;
      break;
    default:
      ERROR_EXIT("invalid dim\n");
      return 0;
    }
  }


  void Element::eraseNewCoord()
  {
    if (newCoord != NULL)
    {
      delete newCoord;
      newCoord = NULL;
    }
  }


  void writeDotFile(Element* el, std::string filename, int maxLevels)
  {
    std::vector<int> graphNodes;
    std::vector<std::pair<int, int>> graphEdges;

    std::vector<Element*> traverseElements, nextTraverseElements;
    traverseElements.push_back(el);

    int nLevels = 0;

    while (!traverseElements.empty() && (maxLevels == -1 || nLevels < maxLevels))
    {
      nextTraverseElements.clear();

      for (unsigned int i = 0; i < traverseElements.size(); i++)
      {
        graphNodes.push_back(traverseElements[i]->getIndex());

        if (!traverseElements[i]->isLeaf() &&
            (maxLevels == -1 || nLevels + 1 < maxLevels))
        {
          graphEdges.push_back(std::make_pair(traverseElements[i]->getIndex(),
                                              traverseElements[i]->getChild(0)->getIndex()));
          graphEdges.push_back(std::make_pair(traverseElements[i]->getIndex(),
                                              traverseElements[i]->getChild(1)->getIndex()));
          nextTraverseElements.push_back(traverseElements[i]->getChild(0));
          nextTraverseElements.push_back(traverseElements[i]->getChild(1));
        }
      }

      traverseElements = nextTraverseElements;
      nLevels++;
    }

    std::ofstream file;
    file.open(filename.c_str());

    file << "digraph G\n";
    file << "{\n";

    for (unsigned int i = 0; i < graphNodes.size(); i++)
      file << "   node" << graphNodes[i]
           << " [ label = \"" << graphNodes[i] << "\"];\n";

    for (unsigned int i = 0; i < graphEdges.size(); i++)
      file << "   \"node" << graphEdges[i].first << "\" -> \"node"
           << graphEdges[i].second << "\";\n";

    file << "}\n";

    file.close();
  }


  void Element::getAllDofs(const FiniteElemSpace* feSpace,
                           BoundaryObject bound,
                           DofContainer& dofs,
                           bool baseDofPtr,
                           std::vector<GeoIndex>* dofGeoIndex)
  {
    FUNCNAME_DBG("Element::getAllDofs()");

    getNodeDofs(feSpace, bound, dofs, baseDofPtr);

    if (dofGeoIndex != NULL)
    {
      // In the case dofGeoIndex should be filled, set all node DOFs to be
      // vertex DOFs.
      dofGeoIndex->resize(dofs.size());
      for (unsigned int i = 0; i < dofs.size(); i++)
        (*dofGeoIndex)[i] = VERTEX;
    }

    if (feSpace->getBasisFcts()->getDegree() > 1)
      getHigherOrderDofs(feSpace, bound, dofs, baseDofPtr, dofGeoIndex);

    if (dofGeoIndex)
    {
      TEST_EXIT_DBG(dofs.size() == dofGeoIndex->size())
      ("Arrays do not fit together: %d %d\n", dofs.size(), dofGeoIndex->size());
    }
  }

} // end namespace AMDiS
