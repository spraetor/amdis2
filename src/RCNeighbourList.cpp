#include "RCNeighbourList.hpp"

#include "CoarseningManager.hpp"
#include "DOFVector.hpp"
#include "Element.hpp"
#include "FixVec.hpp"
#include "LeafData.hpp"
#include "Line.hpp"
#include "MacroElement.hpp"
#include "Mesh.hpp"
#include "PeriodicBC.hpp"
#include "Tetrahedron.hpp"
#include "Traverse.hpp"
#include "Triangle.hpp"

namespace AMDiS
{
  RCNeighbourList::RCNeighbourList(int maxEdgeNeigh)
  {
    rclist.resize(maxEdgeNeigh);
    for (int i = 0; i < maxEdgeNeigh; i++)
      rclist[i] = new RCListElement;
  }


  RCNeighbourList::~RCNeighbourList()
  {
    clearList();
  }


  /****************************************************************************/
  /*  do_coarse_patch:  if patch can be coarsend return true, else false and  */
  /*  reset the element marks                                                 */
  /****************************************************************************/

  bool RCNeighbourList::doCoarsePatch(int n_neigh)
  {
    FUNCNAME("RCNeighbourList::doCoarsePatch()");

    for (int i = 0; i < n_neigh; i++)
    {
      Element* lel = rclist[i]->el;

      if (lel->getMark() >= 0 || lel->isLeaf())
      {
        /****************************************************************************/
        /*  element must not be coarsend or element is a leaf element, reset the    */
        /*  the coarsening flag on all those elements that have to be coarsend with */
        /*  this element                                                            */
        /****************************************************************************/
        lel->setMark(0);
        for (int j = 0; j < n_neigh; j++)
          if (rclist[j]->flag)
            rclist[j]->el->setMark(0);

        return false;
      }
      else if (lel->getFirstChild()->getMark() >= 0 ||
               lel->getSecondChild()->getMark() >= 0)
      {

        /****************************************************************************/
        /*  one of the element's children must not be coarsend; reset the coarsening*/
        /*  flag on all those elements that have to be coarsend with this element   */
        /****************************************************************************/
        lel->setMark(0);
        for (int j = 0; j < n_neigh; j++)
          if (rclist[j]->flag)
            rclist[j]->el->setMark(0);

        return false;
      }
      else if (!lel->getFirstChild()->isLeaf() ||
               !lel->getSecondChild()->isLeaf())
      {
        /****************************************************************************/
        /*  one of the element's children is not a leaf element;                    */
        /*  element may be coarsend after coarsening one of the children; try again */
        /****************************************************************************/
        coarseningManager->doMore = true;

        return false;
      }
      else
      {
        /****************************************************************************/
        /*  either one element is a macro element or we can coarsen the patch       */
        /****************************************************************************/
        if (rclist[i]->flag == 0)
        {
          Mesh* coarse_mesh = coarseningManager->getMesh();
          std::deque<MacroElement*>::const_iterator mel;
          // set in Mesh::coarsen()
          for (mel = coarse_mesh->firstMacroElement();
               mel!=coarse_mesh->endOfMacroElements(); ++mel)
            if ((*mel)->getElement() == lel)
              break;

          TEST_EXIT(mel != coarse_mesh->endOfMacroElements())
          ("incompatible coarsening patch found\n");
        }
      }
    }

    return true;
  }


  void  RCNeighbourList::fillNeighbourRelations(int n_neigh, int /*bound*/)
  {
    // for all RC Elements in list
    for (int i = 0; i < n_neigh; i++)
    {
      Element* el = rclist[i]->el;
      rclist[i]->ith = i;

      // for all directions of refinementEdges =? Sides of Refinement Edge
      for (int dir = 0; dir < 2; dir++)
      {
        int j = 0;

        //Test if other elements in List are Neighbours of active Element
        for (; j < n_neigh; j++)
        {
          Element* neigh = rclist[j]->el;

          //test-neighbour Element is active element -> skip tests for this elements
          if (neigh == el)
            continue;

          //Test if test-neighbour Element is neighbour of active Element
          int k = 0;
          for (; k < 2; k++)
          {
            if (neigh->getDof(2 + k) == el->getDof(3 - dir))
            {
              //Test FACE neighbourhood via oppVertex
              // Test only Neighbourhood of Faces 2 & 3(local Indexing)
              //			Faces next to refinement Edge
              // check works by comparing DOFs (rank Indexing)
              rclist[i]->neigh[dir] = rclist[j];
              rclist[i]->oppVertex[dir] = 3 - k;
              break;
            }
            else
            {
              rclist[i]->neigh[dir] = NULL;
              rclist[i]->oppVertex[dir] = -1;
            }
          }

          if (k < 2)
            break;
        }

        //none of the Elements in rcList is FACE neighbour of active Element
        if (j >= n_neigh)
        {
          rclist[i]->neigh[dir] = NULL;
          rclist[i]->oppVertex[dir] = -1;
        }
      }
    }
  }


  void RCNeighbourList::addDOFParent(int elIndex, DegreeOfFreedom* dof)  // 3d
  {
    Element* el = rclist[elIndex]->el;
    RCListElement* neighbour = NULL;
    Mesh* coarse_mesh = coarseningManager->getMesh();
    RCListElement* coarse_list = rclist[elIndex];

    if (coarse_mesh->getNumberOfDofs(EDGE))
    {
      int node = coarse_mesh->getNode(EDGE);

      /****************************************************************************/
      /*  set the dof in the coarsening edge                                      */
      /****************************************************************************/

      el->setDof(node, dof);

      /****************************************************************************/
      /*  and now those handed on by the children                                 */
      /****************************************************************************/

      el->setDof(node + 1, const_cast<DegreeOfFreedom*>(el->getFirstChild()->getDof(node)));
      el->setDof(node + 2, const_cast<DegreeOfFreedom*>(el->getFirstChild()->getDof(node + 1)));
      el->setDof(node + 5, const_cast<DegreeOfFreedom*>(el->getFirstChild()->getDof(node + 3)));

      if (coarse_list->elType)
      {
        el->setDof(node + 3, const_cast<DegreeOfFreedom*>(el->getSecondChild()->getDof(node)));
        el->setDof(node + 4, const_cast<DegreeOfFreedom*>(el->getSecondChild()->getDof(node + 1)));
      }
      else
      {
        el->setDof(node + 3, const_cast<DegreeOfFreedom*>(el->getSecondChild()->getDof(node + 1)));
        el->setDof(node + 4, const_cast<DegreeOfFreedom*>(el->getSecondChild()->getDof(node)));
      }
    }

    if (coarse_mesh->getNumberOfDofs(FACE))
    {
      int node = coarse_mesh->getNode(FACE);
      /****************************************************************************/
      /*  dof's at the faces within the patch: add new dof's if it is a boundary  */
      /*  face or neighbour is an element behind this element in the corsen list  */
      /****************************************************************************/

      neighbour = coarse_list->neigh[0];
      if (!neighbour || neighbour > coarse_list)
      {
        if (!el->getDof(node + 2))
        {
          // face 2
          el->setDof(node + 2, const_cast<DegreeOfFreedom*>(coarse_mesh->getDof(FACE)));
          if (neighbour)
            neighbour->el->setDof(node + coarse_list->oppVertex[0],
                                  const_cast<DegreeOfFreedom*>(el->getDof(node + 2)));
        }
      }

      neighbour = coarse_list->neigh[1];
      if (!neighbour || neighbour > coarse_list)
      {
        if (!el->getDof(node + 3))
        {
          // face 3
          el->setDof(node + 3, const_cast<DegreeOfFreedom*>(coarse_mesh->getDof(FACE)));
          if (neighbour)
            neighbour->el->setDof(node + coarse_list->oppVertex[1],
                                  const_cast<DegreeOfFreedom*>(el->getDof(node + 3)));
        }
      }
      /****************************************************************************/
      /*  and now those handed on by the children                                 */
      /****************************************************************************/

      el->setDof(node, const_cast<DegreeOfFreedom*>(el->getSecondChild()->getDof(node + 3)));
      el->setDof(node + 1, const_cast<DegreeOfFreedom*>(el->getFirstChild()->getDof(node + 3)));
    }

    if (coarse_mesh->getNumberOfDofs(CENTER))
    {
      int node = coarse_mesh->getNode(CENTER);
      if (!el->getDof(node))
        el->setDof(node, const_cast<DegreeOfFreedom*>(coarse_mesh->getDof(CENTER)));
    }
  }


  void RCNeighbourList::addDOFParents(int n_neigh) // 2d
  {
    Mesh* coarse_mesh = coarseningManager->getMesh();

    if (coarse_mesh->getNumberOfDofs(EDGE))
    {
      int node = coarse_mesh->getNode(EDGE);

      /****************************************************************************/
      /*  get dofs on the boundary of the coarsening patch from the children      */
      /****************************************************************************/
      for (int i = 0; i < n_neigh; i++)
      {
        rclist[i]->el->setDof(node, const_cast<DegreeOfFreedom*>(rclist[i]->el->getSecondChild()->getDof(node + 2)));
        rclist[i]->el->setDof(node + 1, const_cast<DegreeOfFreedom*>(rclist[i]->el->getFirstChild()->getDof(node + 2)));
      }
    }

    if (coarse_mesh->getNumberOfDofs(CENTER))
    {
      int node = coarse_mesh->getNode(CENTER);

      /****************************************************************************/
      /*  get new dof on parents at the barycenter                                */
      /****************************************************************************/
      for (int i = 0; i < n_neigh; i++)
        if (!rclist[i]->el->getDof(node))
          rclist[i]->el->setDof(node, const_cast<DegreeOfFreedom*>(coarse_mesh->getDof(CENTER)));
    }
  }


  void RCNeighbourList::removeDOFParents(int n_neigh)
  {
    Mesh* mesh = rclist[0]->el->getMesh();
    int edges = mesh->getGeo(EDGE);

    if (mesh->getNumberOfDofs(EDGE))
    {
      int node = mesh->getNode(EDGE);

      for (int i = 0; i < n_neigh; i++)
        for (int j = 0; j < edges; j++)
          rclist[i]->el->setDof(node + j, NULL);
    }

    if (mesh->getNumberOfDofs(CENTER))
    {
      int node = mesh->getNode(CENTER);
      for (int i = 0; i < n_neigh; i++)
      {
        mesh->freeDof(const_cast<DegreeOfFreedom*>(rclist[i]->el->getDof(node)), CENTER);
        rclist[i]->el->setDof(node, NULL);
      }
    }
  }


  void RCNeighbourList::removeDOFParent(int index)
  {
    Tetrahedron* el = dynamic_cast<Tetrahedron*>(rclist[index]->el);
    Mesh* mesh = el->getMesh();
    int edges = mesh->getGeo(EDGE);
    int faces = mesh->getGeo(FACE);


    if (mesh->getNumberOfDofs(EDGE))
    {
      int node = mesh->getNode(EDGE);
      for (int j = 0; j < edges; j++)
        el->setDof(node + j, NULL);
    }

    if (mesh->getNumberOfDofs(FACE))
    {
      int node = mesh->getNode(FACE);
      RCListElement* neigh = rclist[index]->neigh[0];

      // face 2
      if (!neigh || neigh > rclist[index])
        mesh->freeDof(const_cast<DegreeOfFreedom*>(el->getDof(node + 2)), FACE);

      neigh = rclist[index]->neigh[1];
      // face 3
      if (!neigh || neigh > rclist[index])
        mesh->freeDof(const_cast<DegreeOfFreedom*>(el->getDof(node + 3)), FACE);

      for (int j = 0; j < faces; j++)
        el->setDof(node + j, NULL);
    }

    if (mesh->getNumberOfDofs(CENTER))
    {
      int node = mesh->getNode(CENTER);
      mesh->freeDof(const_cast<DegreeOfFreedom*>(el->getDof(node)), CENTER);
      el->setDof(node, NULL);
    }
  }


  void RCNeighbourList::periodicSplit(DegreeOfFreedom* edge[2],
                                      DegreeOfFreedom* nextEdge[2],
                                      int* n_neigh,
                                      int* n_neigh_periodic,
                                      RCNeighbourList& periodicList)
  {
    FUNCNAME_DBG("RCNeighbourList::periodicSplit()");

    *n_neigh_periodic = 0;
    int count = 0;
    int n_neigh_old = *n_neigh;
    bool secondPart = false;
    bool firstSplit = true;

    nextEdge[0] = NULL;
    nextEdge[1] = NULL;

    std::vector<RCListElement*>::iterator it = rclist.begin();
    std::vector<RCListElement*>::iterator insertIt;

    while (count < n_neigh_old)
    {
      DegreeOfFreedom* dof0 = const_cast<DegreeOfFreedom*>((*it)->el->getDof(0));
      DegreeOfFreedom* dof1 = const_cast<DegreeOfFreedom*>((*it)->el->getDof(1));

      if (dof0 != edge[0] && dof0 != edge[1])
      {
        secondPart = true;
        if (!nextEdge[0])
        {
          nextEdge[0] = dof0;
          nextEdge[1] = dof1;
        }
        ++it;
      }
      else
      {
        (*n_neigh)--;
        (*n_neigh_periodic)++;

        if (firstSplit)
        {
          periodicList.clearList();
          insertIt = periodicList.rclist.end();
          periodicList.coarseningManager = coarseningManager;
          firstSplit = false;
        }
        else
        {
          if (edge[0])
          {
            TEST_EXIT_DBG((dof0 == edge[0] && dof1 == edge[1]) ||
                          (dof1 == edge[0] && dof0 == edge[1]))
            ("Invalid macro file?\n");
          }
        }

        if (secondPart)
        {
          insertIt = periodicList.rclist.begin();
          secondPart = false;
        }

        insertIt = periodicList.rclist.insert(insertIt, *it);
        ++insertIt;

        it = rclist.erase(it);
      }
      count++;
    }

    for (int i = 0; i < *n_neigh_periodic; i++)
      periodicList.rclist[i]->ith = i;
  }


  void RCNeighbourList::clearList()
  {
    for (unsigned int i = 0; i < rclist.size(); i++)
      delete rclist[i];

    rclist.clear();
  }

} // end namespace AMDiS
