#include "RefinementManager.h"
#include "Mesh.h"
#include "Traverse.h"
#include "ElInfo.h"
#include "DOFAdmin.h"
#include "AdaptStationary.h"
#include "AdaptInstationary.h"
#include "FixVec.h"
#include "RCNeighbourList.h"
#include "ProblemStatBase.h"
#include "DOFIndexed.h"
#include "Projection.h"
#include "LeafData.h"
#include "VertexVector.h"
#include "Debug.h"

namespace AMDiS
{
  ElInfo* RefinementManager2d::refineFunction(ElInfo* elInfo)
  {
    if (elInfo->getElement()->getMark() <= 0)
      return elInfo;

    bool bound = false;
    DegreeOfFreedom* edge[2];
    RCNeighbourList refineList(2);
    refineList.setElement(0, elInfo->getElement());
    int n_neigh = 1;

    if (elInfo->getProjection(0) &&
        elInfo->getProjection(0)->getType() == VOLUME_PROJECTION)
      newCoords = true;

    // === Give the refinement edge the right orientation. ===

    if (elInfo->getElement()->getDof(0, 0) < elInfo->getElement()->getDof(1, 0))
    {
      edge[0] = const_cast<DegreeOfFreedom*>(elInfo->getElement()->getDof(0));
      edge[1] = const_cast<DegreeOfFreedom*>(elInfo->getElement()->getDof(1));
    }
    else
    {
      edge[1] = const_cast<DegreeOfFreedom*>(elInfo->getElement()->getDof(0));
      edge[0] = const_cast<DegreeOfFreedom*>(elInfo->getElement()->getDof(1));
    }

    // === Get the refinement patch. ===

    getRefinePatch(&elInfo, edge, 0, refineList, &n_neigh);


    // === Check for periodic boundary ===

    DegreeOfFreedom* next_edge[2] = {NULL, NULL};
    DegreeOfFreedom* first_edge[2] = {edge[0], edge[1]};
    DegreeOfFreedom* last_edge[2] = {NULL, NULL};
    int n_neigh_periodic;

    DegreeOfFreedom newDOF = -1;
    DegreeOfFreedom lastNewDOF = -1;
    DegreeOfFreedom firstNewDOF = -1;

    RCNeighbourList periodicList;

    while (edge[0] != NULL)
    {
      refineList.periodicSplit(edge, next_edge,
                               &n_neigh, &n_neigh_periodic,
                               periodicList);

      newDOF = refinePatch(edge, periodicList, n_neigh_periodic, bound);

      if (firstNewDOF == -1)
        firstNewDOF = newDOF;

      if (lastNewDOF != -1)
      {
        for (std::map<int, VertexVector*>::iterator it = mesh->getPeriodicAssociations().begin();
             it != mesh->getPeriodicAssociations().end(); ++it)
        {
          if (it->second)
          {
            if (((*(it->second))[edge[0][0]] == last_edge[0][0] &&
                 (*(it->second))[edge[1][0]] == last_edge[1][0]) ||
                ((*(it->second))[edge[0][0]] == last_edge[1][0] &&
                 (*(it->second))[edge[1][0]] == last_edge[0][0]))
            {
              (*(it->second))[lastNewDOF] = newDOF;
              (*(it->second))[newDOF] = lastNewDOF;
            }
          }
        }
      }
      lastNewDOF = newDOF;

      last_edge[0] = edge[0];
      last_edge[1] = edge[1];

      edge[0] = next_edge[0];
      edge[1] = next_edge[1];
    }

    if (lastNewDOF != firstNewDOF)
    {
      for (std::map<int, VertexVector*>::iterator it = mesh->getPeriodicAssociations().begin();
           it != mesh->getPeriodicAssociations().end(); ++it)
      {
        if (it->second)
        {
          if (((*(it->second))[first_edge[0][0]] == last_edge[0][0] &&
               (*(it->second))[first_edge[1][0]] == last_edge[1][0]) ||
              ((*(it->second))[first_edge[0][0]] == last_edge[1][0] &&
               (*(it->second))[first_edge[1][0]] == last_edge[0][0]))
          {
            (*(it->second))[lastNewDOF] = firstNewDOF;
            (*(it->second))[firstNewDOF] = lastNewDOF;
          }
        }
      }
    }

    return elInfo;
  }


  void RefinementManager2d::newCoordsFct(ElInfo* elInfo)
  {
    Element* el = elInfo->getElement();
    int dow = Global::getGeo(WORLD);

    Projection* projector = elInfo->getProjection(0);

    if (!projector || projector->getType() != VOLUME_PROJECTION)
      projector = elInfo->getProjection(2);

    if (el->getFirstChild() && projector && (!el->isNewCoordSet()))
    {
      WorldVector<double>* new_coord = new WorldVector<double>;

      for (int j = 0; j < dow; j++)
        (*new_coord)[j] = (elInfo->getCoord(0)[j] + elInfo->getCoord(1)[j]) * 0.5;

      projector->project(*new_coord);
      el->setNewCoord(new_coord);
    }
  }


  void RefinementManager2d::setNewCoords(int macroEl)
  {
    TraverseStack stack;
    ElInfo* elInfo;

    if (macroEl == -1)
      elInfo = stack.traverseFirst(mesh, -1,
                                   Mesh::CALL_EVERY_EL_PREORDER |
                                   Mesh::FILL_BOUND | Mesh::FILL_COORDS);
    else
      elInfo = stack.traverseFirstOneMacro(mesh, macroEl, -1,
                                           Mesh::CALL_EVERY_EL_PREORDER |
                                           Mesh::FILL_BOUND | Mesh::FILL_COORDS);

    while (elInfo)
    {
      newCoordsFct(elInfo);
      elInfo = stack.traverseNext(elInfo);
    }
  }


  DegreeOfFreedom RefinementManager2d::refinePatch(DegreeOfFreedom* edge[2],
      RCNeighbourList& refineList,
      int n_neigh, bool bound)
  {
    DegreeOfFreedom* dof[3] = {NULL, NULL, NULL};
    Triangle* el =
      dynamic_cast<Triangle*>(const_cast<Element*>(refineList.getElement(0)));
    Triangle* neigh =
      dynamic_cast<Triangle*>(const_cast<Element*>(refineList.getElement(1)));

    // === There is one new vertex in the refinement edge. ===

    dof[0] = mesh->getDof(VERTEX);

    mesh->incrementNumberOfVertices(1);
    mesh->incrementNumberOfEdges(1);

    if (mesh->getNumberOfDofs(EDGE))
    {
      // There are two additional dofs in the refinement edge.
      dof[1] = mesh->getDof(EDGE);
      dof[2] = mesh->getDof(EDGE);
    }

    // === First refine the element. ===

    bisectTriangle(el, dof);

    if (neigh)
    {
      DegreeOfFreedom* tmp = dof[1];

      // === There is a neighbour; refine it also, but first exchange dof[1] and ===
      // === dof[2]; thus, dof[1] is always added on child[0]!                   ===
      dof[1] = dof[2];
      dof[2] = tmp;

      bisectTriangle(neigh, dof);
    }
    else
    {
      newCoords = true;
    }

    // === If there are functions to interpolate data to the finer grid, do so.

    int nrAdmin = mesh->getNumberOfDOFAdmin();
    for(int iadmin = 0; iadmin < nrAdmin; iadmin++)
    {
      DOFAdmin* admin = const_cast<DOFAdmin*>(&mesh->getDofAdmin(iadmin));
      std::list<DOFIndexedBase*>::iterator end = admin->endDOFIndexed();
      for (std::list<DOFIndexedBase*>::iterator it = admin->beginDOFIndexed();
           it != end; it++)
        (*it)->refineInterpol(refineList, n_neigh);
    }


    if (!mesh->queryCoarseDOFs())
    {
      // === If there should be no dof information on interior leaf elements ===
      // === remove dofs from edges and the centers of parents.              ===
      if (mesh->getNumberOfDofs(EDGE))
      {
        int node = mesh->getNode(EDGE);

        // === The only DOF that can be freed is that in the refinement edge; all ===
        // === other DOFs are handed on the children.                             ===

        mesh->freeDof(const_cast<DegreeOfFreedom*>(el->getDof(node+2)), EDGE);
      }
      if (mesh->getNumberOfDofs(EDGE) || mesh->getNumberOfDofs(CENTER))
        refineList.removeDOFParents(n_neigh);
    }

    return dof[0][0];
  }


  void RefinementManager2d::bisectTriangle(Triangle* el, DegreeOfFreedom* newDOFs[3])
  {
    FUNCNAME_DBG("RefinementManager2d::bisectTriangle()");

    TEST_EXIT_DBG(mesh)("No mesh!\n");

    Triangle* child[2];
    child[0] = dynamic_cast<Triangle*>(mesh->createNewElement(el));
    child[1] = dynamic_cast<Triangle*>(mesh->createNewElement(el));

    int newMark = std::max(0, el->getMark() - 1);
    child[0]->setMark(newMark);
    child[1]->setMark(newMark);
    el->setMark(0);


    // === Transfer hidden data from parent to children. ===

    // call of subclass-method
    el->refineElementData(child[0], child[1]);
    el->setFirstChild(child[0]);
    el->setSecondChild(child[1]);

    if (newMark > 0)
      doMoreRecursiveRefine = true;


    // === Vertex 2 is the newest vertex. ===

    child[0]->setDof(2, newDOFs[0]);
    child[1]->setDof(2, newDOFs[0]);


    // === The other vertices are handed on from the parent. ===

    for (int i_child = 0; i_child < 2; i_child++)
    {
      child[i_child]->setDof(i_child, const_cast<DegreeOfFreedom*>(el->getDof(2)));
      child[i_child]->setDof(1 - i_child, const_cast<DegreeOfFreedom*>(el->getDof(i_child)));
    }


    // === There is one more leaf element, two hierachical elements and one ===
    // === more edge.                                                       ===

    mesh->incrementNumberOfEdges(1);
    mesh->incrementNumberOfLeaves(1);
    mesh->incrementNumberOfElements(2);

    if (mesh->getNumberOfDofs(EDGE))
    {
      DegreeOfFreedom* newEdgeDOFs = mesh->getDof(EDGE);

      // There are dofs in the midpoint of the edges.
      child[0]->setDof(4, newEdgeDOFs);
      child[1]->setDof(3, newEdgeDOFs);

      // Dofs handed on by the parent.
      child[0]->setDof(5, const_cast<DegreeOfFreedom*>(el->getDof(4)));
      child[1]->setDof(5, const_cast<DegreeOfFreedom*>(el->getDof(3)));

      // Dofs in the refinement edge.
      child[0]->setDof(3, newDOFs[1]);
      child[1]->setDof(4, newDOFs[2]);
    }

    if (mesh->getNumberOfDofs(CENTER))
    {
      int node = mesh->getNode(CENTER);

      // There are dofs at the barycenter of the triangles.
      child[0]->setDof(node, mesh->getDof(CENTER));
      child[1]->setDof(node, mesh->getDof(CENTER));
    }
  }


  void RefinementManager2d::getRefinePatch(ElInfo** elInfo,
      DegreeOfFreedom* edge[2],
      int dir, RCNeighbourList& refineList,
      int* n_neigh)
  {
    FUNCNAME_DBG("RefinementManager2d::getRefinePatch()");

    if ((*elInfo)->getNeighbour(2) && (*elInfo)->getOppVertex(2) != 2)
    {
      // Neighbour is not compatible devisible; refine neighbour first, store the
      // opp_vertex to traverse back to el.
      int opp_vertex = (*elInfo)->getOppVertex(2);

      ElInfo* neigh_info = stack->traverseNeighbour2d(*elInfo, 2);
      neigh_info->getElement()->setMark(std::max(neigh_info->getElement()->getMark(), 1));
      neigh_info = refineFunction(neigh_info);

      // === Now go back to the original element and refine the patch. ===


      // In Debug mode we test if traverNeighbour2d returns the expected element.
      int testIndex = neigh_info->getNeighbour(opp_vertex)->getIndex();

      *elInfo = neigh_info = stack->traverseNeighbour2d(neigh_info, opp_vertex);


      TEST_EXIT_DBG(testIndex == (*elInfo)->getElement()->getIndex())
      ("Got wrong neighbour! Should be %d, but is %d!\n",
       testIndex, (*elInfo)->getElement()->getIndex());

      TEST_EXIT_DBG(neigh_info->getElement() ==
                    dynamic_cast<Triangle*>(const_cast<Element*>((*elInfo)->getElement())))
      ("invalid traverse_neighbour1\n");
    }

    if (refineList.setElement(1, (*elInfo)->getNeighbour(2)))
    {
      TEST_EXIT_DBG((*elInfo)->getOppVertex(2) == 2)
      ("no compatible ref. edge after recursive refinement of neighbour\n");
      *n_neigh = 2;
    }
  }

} // end namespace AMDiS
