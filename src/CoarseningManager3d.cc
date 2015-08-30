#include "CoarseningManager3d.h"
#include "Mesh.h"
#include "AdaptStationary.h"
#include "AdaptInstationary.h"
#include "Traverse.h"
#include "MacroElement.h"
#include "RCNeighbourList.h"
#include "FixVec.h"
#include "DOFIndexed.h"
#include "Debug.h"

using namespace std;

namespace AMDiS
{

  void CoarseningManager3d::coarsenFunction(ElInfo* elInfo)
  {
#if HAVE_PARALLEL_DOMAIN_AMDIS
    FUNCNAME_DBG("CoarseningManager3d::coarsenFunction()");
#endif

    Tetrahedron* el =
      dynamic_cast<Tetrahedron*>(const_cast<Element*>(elInfo->getElement()));

    // If element must not be coarsend, return.
    if (el->getMark() >= 0)
      return;

    // Single leaves don't get coarsened.
    if (el->isLeaf())
      return;

    if (el->getChild(0)->getMark() >= 0 || el->getChild(1)->getMark() >= 0)
    {
      // One of the children must not be coarsend, so return.

      el->setMark(0);
      return;
    }

    if (!(el->getChild(0)->isLeaf()) || !(el->getChild(1)->isLeaf()))
    {
      // One of the children is not a leaf element. Try again later on.
      doMore = true;
      return;
    }

    DegreeOfFreedom* edge[2];
    int n_neigh, bound = 0;

    /****************************************************************************/
    /*  get a list for storing all elements at the coarsening edge and fill it  */
    /****************************************************************************/
    RCNeighbourList coarsenList(mesh->getMaxEdgeNeigh());
    coarsenList.setCoarseningManager(this);

    /****************************************************************************/
    /*  give the refinement edge the right orientation                          */
    /****************************************************************************/

    if (el->getDof(0, 0) < el->getDof(1, 0))
    {
      edge[0] = const_cast<DegreeOfFreedom*>(el->getDof(0));
      edge[1] = const_cast<DegreeOfFreedom*>(el->getDof(1));
    }
    else
    {
      edge[1] = const_cast<DegreeOfFreedom*>(el->getDof(0));
      edge[0] = const_cast<DegreeOfFreedom*>(el->getDof(1));
    }

    coarsenList.setElement(0, el, true);
    n_neigh = 1;

    coarsenList.setOppVertex(0, 0, 0);
    coarsenList.setElType(0, elInfo->getType());
    bound = false;
    if (getCoarsenPatch(elInfo, edge, 0, coarsenList, &n_neigh))
    {
      getCoarsenPatch(elInfo, edge, 1, coarsenList, &n_neigh);
      bound = true;
    }


#if HAVE_PARALLEL_DOMAIN_AMDIS
    vector<FixRefinementPatch::EdgeInEl> refineEdges;
    FixRefinementPatch::getOtherEl(mesh, stack, refineEdges);

    // === If the refinement edge must be fixed, add also the other part of this ===
    // === edge to the refinement patch.                                         ===

    if (refineEdges.size())
    {
      // TODO: Remove these two lines and make something more meaningful!!
      el->setMark(0);
      return;

      Element* otherEl = refineEdges[0].first;
      TEST_EXIT_DBG(otherEl->getMesh() == mesh)("Something is wrong.\n");

      TraverseStack stack2;
      ElInfo* elInfo2 =
        stack2.traverseFirstOneMacro(mesh, otherEl->getIndex(), -1,
                                     Mesh::CALL_EVERY_EL_PREORDER |
                                     Mesh::FILL_NEIGH |
                                     Mesh::FILL_BOUND);

      bool foundEdge = false;

      while (elInfo2)
      {
        Element* el2 = elInfo2->getElement();
        for (int i = 0; i < 6; i++)
        {
          DofEdge edge2 = elInfo2->getElement()->getEdge(i);
          if (edge2.first == *(edge[0]) && edge2.second == *(edge[1]) && !el2->isLeaf())
          {

            if (!el2->getChild(0)->isLeaf() || !el2->getChild(1)->isLeaf())
            {
              int edgeNo0 = el->getEdgeOfChild(0, i, elInfo2->getType());
              int edgeNo1 = el->getEdgeOfChild(1, i, elInfo2->getType());

              bool refineChildFirst =
                !((i > 0 && (edgeNo0 >= 0 && !el2->getChild(0)->isLeaf())) ||
                  (edgeNo1 >= 0 && !el2->getChild(1)->isLeaf()));

              if (refineChildFirst)
              {
                // TODO: WAS SOLL ICH NUR MACHEN???
                el->setMark(0);
                return;
              }
            }
            else
            {
              coarsenList.setElType(n_neigh, elInfo2->getType());
              coarsenList.setElement(n_neigh, elInfo2->getElement(), true);
              n_neigh++;

              foundEdge = true;

              TraverseStack* tmpStack = stack;
              stack = &stack2;

              if (getCoarsenPatch(elInfo2, edge, 0, coarsenList, &n_neigh))
              {
                getCoarsenPatch(elInfo2, edge, 1, coarsenList, &n_neigh);
                bound = true;
              }
              stack = tmpStack;
              break;
            }
          }
        }

        elInfo2 = stack2.traverseNext(elInfo2);
      }

      TEST_EXIT_DBG(foundEdge)("Should not happen!\n");
    }
#endif

    coarsenList.fillNeighbourRelations(n_neigh, bound);

    /****************************************************************************/
    /*  check wether we can coarsen the patch or not                            */
    /****************************************************************************/

    // ==========================================================================
    // === check for periodic boundary ==========================================
    // ==========================================================================

    if (coarsenList.doCoarsePatch(n_neigh))
    {
      int n_neigh_periodic;
      DegreeOfFreedom* next_edge[2];
      RCNeighbourList periodicList;

      while (edge[0] != NULL)
      {
        coarsenList.periodicSplit(edge, next_edge,
                                  &n_neigh, &n_neigh_periodic,
                                  periodicList);

        coarsenPatch(periodicList, n_neigh_periodic, bound);

        edge[0] = next_edge[0];
        edge[1] = next_edge[1];
      }
    }
  }

  /*****************************************************************************/
  /*  coarsenTetrahedron:  coarses a single element of the coarsening patch;   */
  /*  dofs in the interior of the element are removed and dofs in the faces of */
  /*  the element are removed if neighbour has been coarsend or if the face    */
  /*  is part of the domains boundary                                          */
  /*****************************************************************************/

  void CoarseningManager3d::coarsenTetrahedron(RCNeighbourList& coarsenList,
      int index)
  {
    Tetrahedron* el =
      dynamic_cast<Tetrahedron*>(const_cast<Element*>(coarsenList.getElement(index)));
    Tetrahedron* child[2];

    child[0] = dynamic_cast<Tetrahedron*>(const_cast<Element*>(el->getChild(0)));
    child[1] = dynamic_cast<Tetrahedron*>(const_cast<Element*>(el->getChild(1)));
    int el_type = coarsenList.getType(index);

    /****************************************************************************/
    /*  Information about patch neighbours is still valid! But edge and face    */
    /*  dof's in a common face of patch neighbours have to be removed           */
    /****************************************************************************/

    for (int dir = 0; dir < 2; dir++)
    {
      Tetrahedron* neigh =
        dynamic_cast<Tetrahedron*>(const_cast<Element*>(coarsenList.getNeighbourElement(index, dir)));

      if (!neigh || neigh->isLeaf())
      {
        /****************************************************************************/
        /*  boundary face or  neigh has been coarsend: free the dof's in the face   */
        /****************************************************************************/

        if (mesh->getNumberOfDofs(EDGE))
        {
          int node = mesh->getNode(EDGE) + Tetrahedron::nChildEdge[el_type][0][dir];
          mesh->freeDof(const_cast<DegreeOfFreedom*>(child[0]->getDof(node)), EDGE);
        }
        if (mesh->getNumberOfDofs(FACE))
        {
          int node = mesh->getNode(FACE) + Tetrahedron::nChildFace[el_type][0][dir];
          mesh->freeDof(const_cast<DegreeOfFreedom*>(child[0]->getDof(node)), FACE);
          node = mesh->getNode(FACE) + Tetrahedron::nChildFace[el_type][1][dir];
          mesh->freeDof(const_cast<DegreeOfFreedom*>(child[1]->getDof(node)), FACE);
        }
      }
    }

    /****************************************************************************/
    /*  finally remove the interior dof's: in the common face of child[0] and   */
    /*  child[1] and at the two barycenter                                      */
    /****************************************************************************/

    if (mesh->getNumberOfDofs(FACE))
    {
      int node = mesh->getNode(FACE);
      mesh->freeDof(const_cast<DegreeOfFreedom*>(child[0]->getDof(node)), FACE);
    }


    if (mesh->getNumberOfDofs(CENTER))
    {
      int node = mesh->getNode(CENTER);
      for (int i = 0; i < 2; i++)
        mesh->freeDof(const_cast<DegreeOfFreedom*>(child[i]->getDof(node)), CENTER);
    }

    /****************************************************************************/
    /*  get new data on parent and transfer data from children to parent        */
    /****************************************************************************/

    el->coarsenElementData(child[0], child[1], el_type);

    el->setFirstChild(NULL);
    el->setSecondChild(NULL);

    mesh->freeElement(child[0]);
    mesh->freeElement(child[1]);

    el->incrementMark();

    mesh->incrementNumberOfLeaves(-1);
    mesh->incrementNumberOfElements(-2);
  }

  /****************************************************************************/
  /*  get_coarse_patch:  gets the patch for coarsening starting on element    */
  /*  elInfo->el in direction of neighbour [3-dir]; returns 1 if a boundary  */
  /*  reached and 0 if we come back to the starting element.                  */
  /*                                                                          */
  /*  if NEIGH_IN_EL we only can find the complete coarsening patch if the    */
  /*  can be coarsend; otherwise neighbour information is not valid for       */
  /*  parents; in such situation we stop looping around the edge and return 0 */
  /*                                                                          */
  /*  if !NEIGH_IN_EL we complete the loop also in the case of a incompatible */
  /*  coarsening patch since then all marks of patch elements are reset by    */
  /*  do_coarse_patch() and this minimizes calls of traverse_neighbour();     */
  /*  if we reach a boundary while looping around the edge we loop back to    */
  /*  the starting element before we return                                   */
  /****************************************************************************/

  bool CoarseningManager3d::getCoarsenPatch(ElInfo* elInfo,
      DegreeOfFreedom* edge[2],
      int dir,
      RCNeighbourList& coarsenList,
      int* n_neigh)
  {
    FUNCNAME_DBG("CoarseningManager3d::getCoarsenPatch()");

    static const unsigned char next_el[6][2] = {{3,2},
      {1,3},
      {1,2},
      {0,3},
      {0,2},
      {0,1}
    };

    Tetrahedron* el =
      dynamic_cast<Tetrahedron*>(const_cast<Element*>(elInfo->getElement()));
    Tetrahedron* neigh =
      dynamic_cast<Tetrahedron*>(const_cast<Element*>(elInfo->getNeighbour(3 - dir)));
    if (neigh == NULL)
      return true;

    int opp_v = elInfo->getOppVertex(3 - dir);
    ElInfo* neigh_info = stack->traverseNeighbour3d(elInfo, 3 - dir);

    TEST_EXIT_DBG(neigh == neigh_info->getElement())
    ("neigh %d and neigh_info->el %d are not identical\n",
     neigh->getIndex(), neigh_info->getElement()->getIndex());
    /****************************************************************************/
    /*  we have to go back to the starting element via opp_v values             */
    /*  correct information is produce by get_neigh_on_patch()                  */
    /****************************************************************************/
    coarsenList.setOppVertex(*n_neigh, 0, opp_v);
    coarsenList.setElement(*n_neigh, neigh);
    coarsenList.setElType(*n_neigh, neigh_info->getType());

    int nVertices = mesh->getGeo(VERTEX);

    while (neigh != el)
    {
      // === Find the coarsening edge on the current neighbour element. ===

      int i, j, k;

      // === First, try to identify the edge DOFs directly. ===

      for (j = 0; j < nVertices; j++)
        if (neigh->getDof(j) == edge[0])
          break;

      for (k = 0; k < nVertices; k++)
        if (neigh->getDof(k) == edge[1])
          break;

      // === If one of the edge DOFs was not found, try to make use of periodic ===
      // === DOF associations. First, use the direct associations between DOFs. ===
      // === If this will not work, continue with testing on indirect           ===
      // === associations.                                                      ===

      if (j >= nVertices)
      {
        for (j = 0; j < nVertices; j++)
          if (mesh->associated(neigh->getDof(j, 0), edge[0][0]))
            break;

        if (j >= nVertices)
          for (j = 0; j < nVertices; j++)
            if (mesh->indirectlyAssociated(neigh->getDof(j, 0), edge[0][0]))
              break;

        TEST_EXIT_DBG(j < nVertices)
        ("Process element %d: DOF %d not found on element %d with nodes (%d %d %d %d)\n",
         el->getIndex(), edge[0][0], neigh->getIndex(), neigh->getDof(0, 0),
         neigh->getDof(1, 0), neigh->getDof(2, 0), neigh->getDof(3, 0));
      }

      if (k >= nVertices)
      {
        for (k = 0; k < nVertices; k++)
          if (mesh->associated(neigh->getDof(k, 0), edge[1][0]))
            break;

        if (k >= nVertices)
          for (k = 0; k < nVertices; k++)
            if (mesh->indirectlyAssociated(neigh->getDof(k, 0), edge[1][0]))
              break;

        TEST_EXIT_DBG(k < nVertices)
        ("Process element %d: DOF %d not found on element %d with nodes (%d %d %d %d)\n",
         el->getIndex(), edge[1][0], neigh->getIndex(), neigh->getDof(0, 0),
         neigh->getDof(1, 0), neigh->getDof(2, 0), neigh->getDof(3, 0));
      }

      int edgeNo = Tetrahedron::edgeOfDofs[j][k];
      coarsenList.setCoarsePatch(*n_neigh, edgeNo == 0);


      // === Get the direction of the next neighbour. ===

      if (next_el[edgeNo][0] != opp_v)
        i = next_el[edgeNo][0];
      else
        i = next_el[edgeNo][1];

      ++*n_neigh;

      opp_v = neigh_info->getOppVertex(i);
      neigh =
        dynamic_cast<Tetrahedron*>(const_cast<Element*>(neigh_info->getNeighbour(i)));
      if (neigh)
      {
        neigh_info = stack->traverseNeighbour3d(neigh_info, i);
        TEST_EXIT_DBG(neigh == neigh_info->getElement())
        ("neigh %d and neigh_info->el %d are not identical\n",
         neigh->getIndex(), neigh_info->getElement()->getIndex());
        /****************************************************************************/
        /*  we have to go back to the starting element via opp_v values             */
        /*  correct information is produce by get_neigh_on_patch()                  */
        /****************************************************************************/
        coarsenList.setOppVertex(*n_neigh, 0, opp_v);
        coarsenList.setElement(*n_neigh, neigh);
        coarsenList.setElType(*n_neigh, neigh_info->getType());
      }
      else
      {
        break;
      }
    }

    if (neigh == el)
      return false;

    /****************************************************************************/
    /*  the domain's boundary is reached; loop back to the starting el          */
    /****************************************************************************/

    int i = *n_neigh - 1;
    opp_v = coarsenList.getOppVertex(i, 0);
    do
    {
      TEST_EXIT_DBG(neigh_info->getNeighbour(opp_v)  &&  i > 0)
      ("while looping back domains boundary was reached or i == 0\n");
      opp_v = coarsenList.getOppVertex(i--, 0);
      neigh_info = stack->traverseNeighbour3d(neigh_info, opp_v);
    }
    while (neigh_info->getElement() != el);

    return true;
  }

  /****************************************************************************/
  /*  coarse_patch: first rebuild the dofs on the parents then do restriction */
  /*  of data (if possible) and finally coarsen the patch elements            */
  /****************************************************************************/

  void CoarseningManager3d::coarsenPatch(RCNeighbourList& coarsenList,
                                         int n_neigh,
                                         int bound)
  {
    FUNCNAME_DBG("CoarseningManager3d::coarsenPatch()");

    Tetrahedron* el =
      dynamic_cast<Tetrahedron*>(const_cast<Element*>(coarsenList.getElement(0)));
    DegreeOfFreedom* dof = NULL;

    TEST_EXIT_DBG(el)("No element!\n");
    TEST_EXIT_DBG(el->getChild(0))("No child in element!\n");

    if (mesh->getNumberOfDofs(EDGE))
    {
      // === Get dof for coarsening edge. ===

      int node = mesh->getNode(EDGE);
      if (!(dof = const_cast<DegreeOfFreedom*>(el->getDof(node))))
        dof = mesh->getDof(EDGE);
    }

    if (mesh->getNumberOfDofs(EDGE) ||
        mesh->getNumberOfDofs(FACE) ||
        mesh->getNumberOfDofs(CENTER))
    {
      for (int i = 0; i < n_neigh; i++)
        coarsenList.addDOFParent(i, dof);
    }


    // === Restrict dof vectors to the parents on the patch. ===

    int nrAdmin = mesh->getNumberOfDOFAdmin();
    for (int iadmin = 0; iadmin < nrAdmin; iadmin++)
    {
      std::list<DOFIndexedBase*>::iterator it;
      DOFAdmin* admin = const_cast<DOFAdmin*>(&mesh->getDofAdmin(iadmin));
      std::list<DOFIndexedBase*>::iterator end = admin->endDOFIndexed();
      for (it = admin->beginDOFIndexed(); it != end; ++it)
        (*it)->coarseRestrict(coarsenList, n_neigh);
    }

    // === And now start to coarsen the patch: remove dof's of the coarsening edge. ===

    mesh->freeDof(const_cast<DegreeOfFreedom*>(el->getChild(0)->getDof(3)), VERTEX);
    mesh->incrementNumberOfVertices(-1);

    if (mesh->getNumberOfDofs(EDGE))
    {
      int node = mesh->getNode(EDGE) + 2;
      mesh->freeDof(const_cast<DegreeOfFreedom*>(el->getChild(0)->getDof(node)), EDGE);
      mesh->freeDof(const_cast<DegreeOfFreedom*>(el->getChild(1)->getDof(node)), EDGE);
    }

    if (coarsenList.getElement(0)->isNewCoordSet())
      coarsenList.getElement(0)->eraseNewCoord();

    for (int i = 0; i < n_neigh; i++)
    {
      coarsenList.getElement(i)->setNewCoord(NULL);
      coarsenTetrahedron(coarsenList, i);
    }

    // === If the coarsening edge is an interior edge there are  n_neigh + 1 edges ===
    // === and 2 * n_neigh + 1 faces removed; if it is a boundary edge it is one   ===
    // === more edge and one more face.                                            ===

    if (bound)
    {
      mesh->incrementNumberOfEdges(-(n_neigh + 2));
      mesh->incrementNumberOfFaces(-(2 * n_neigh + 1));
    }
    else
    {
      mesh->incrementNumberOfEdges(-(n_neigh + 1));
      mesh->incrementNumberOfFaces(-(2 * n_neigh));
    }
  }

} // end namespace AMDiS
