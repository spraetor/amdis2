#include "CoarseningManager2d.h"
#include "Mesh.h"
#include "AdaptStationary.h"
#include "AdaptInstationary.h"
#include "Traverse.h"
#include "MacroElement.h"
#include "RCNeighbourList.h"
#include "FixVec.h"
#include "DOFIndexed.h"
#include "ProblemStatBase.h"   // => MESH_COARSENED

namespace AMDiS
{

  /****************************************************************************/
  /*  coarseTriangle:  coarses a single element of the coarsening patch; dofs */
  /*  in the interior of the element are removed; dofs for higher order       */
  /*  at the boundary or the coarsening patch still belong to                 */
  /*  the parent. Do not remove them form the mesh!!!                         */
  /****************************************************************************/

  void CoarseningManager2d::coarsenTriangle(Triangle* el)
  {
    FUNCNAME_DBG("CoarseningManager2d::coarseTriangle()");

    Triangle* child[2];
    child[0] = dynamic_cast<Triangle*>(const_cast<Element*>(el->getChild(0)));
    child[1] = dynamic_cast<Triangle*>(const_cast<Element*>(el->getChild(1)));

    TEST_EXIT_DBG(child[0]->getMark() < 0  &&  child[1]->getMark() < 0)
    ("element %d with children[%d,%d] must not be coarsend!\n",
     el->getIndex(), child[0]->getIndex(), child[1]->getIndex());

    // remove dof from common edge of child[0] and child[1]
    if (mesh->getNumberOfDofs(EDGE))
      mesh->freeDof(const_cast<DegreeOfFreedom*>(child[0]->getDof(4)), EDGE);

    // remove dof from the barycenters of child[0] and child[1]
    if (mesh->getNumberOfDofs(CENTER))
    {
      int node = mesh->getNode(CENTER);

      mesh->freeDof(const_cast<DegreeOfFreedom*>(child[0]->getDof(node)), CENTER);
      mesh->freeDof(const_cast<DegreeOfFreedom*>(child[1]->getDof(node)), CENTER);
    }

    el->coarsenElementData(child[0], child[1]);

    el->setFirstChild(NULL);
    el->setSecondChild(NULL);

    mesh->freeElement(child[0]);
    mesh->freeElement(child[1]);

    el->incrementMark();

    mesh->incrementNumberOfLeaves(-1);
    mesh->incrementNumberOfElements(-2);
    mesh->incrementNumberOfEdges(-1);
  }

  /****************************************************************************/
  /*  coarsenPatch: first rebuild the dofs on the parents then do restriction */
  /*  of data (if possible) and finally coarsen the patch elements            */
  /****************************************************************************/

  void CoarseningManager2d::coarsenPatch(RCNeighbourList& coarsenList,
                                         int n_neigh,
                                         int /*bound*/)
  {
    Triangle* el =
      dynamic_cast<Triangle*>(const_cast<Element*>(coarsenList.getElement(0)));
    Triangle* neigh =
      dynamic_cast<Triangle*>(const_cast<Element*>(coarsenList.getElement(1)));
    DegreeOfFreedom* dof[3];

    dof[0] = const_cast<DegreeOfFreedom*>(el->getChild(0)->getDof(2));
    if (mesh->getNumberOfDofs(EDGE))
    {
      dof[1] = const_cast<DegreeOfFreedom*>(el->getChild(0)->getDof(3));
      dof[2] = const_cast<DegreeOfFreedom*>(el->getChild(1)->getDof(4));
    }
    else
    {
      dof[1] = dof[2] = 0;
    }

    if (mesh->getNumberOfDofs(EDGE))
    {
      int node = mesh->getNode(EDGE);
      // get new dof on el at the midpoint of the coarsening edge

      if (!el->getDof(node + 2))
      {
        el->setDof(node + 2, mesh->getDof(EDGE));
        if (neigh)
          neigh->setDof(node + 2, const_cast<DegreeOfFreedom*>(el->getDof(node + 2)));
      }
    }

    if (mesh->getNumberOfDofs(EDGE) || mesh->getNumberOfDofs(CENTER))
      coarsenList.addDOFParents(n_neigh);

    // restrict dof vectors to the parents on the patch

    int nrAdmin = mesh->getNumberOfDOFAdmin();
    for (int iadmin = 0; iadmin < nrAdmin; iadmin++)
    {
      DOFAdmin* admin = const_cast<DOFAdmin*>(&mesh->getDofAdmin(iadmin));
      for (std::list<DOFIndexedBase*>::iterator it = admin->beginDOFIndexed();
           it != admin->endDOFIndexed(); ++it)
        (*it)->coarseRestrict(coarsenList, n_neigh);
    }

    coarsenTriangle(el);

    if (neigh)
      coarsenTriangle(neigh);

    // now, remove those dofs in the coarcening edge

    mesh->freeDof(dof[0], VERTEX);
    if (mesh->getNumberOfDofs(EDGE))
    {
      mesh->freeDof(dof[1], EDGE);
      mesh->freeDof(dof[2], EDGE);
    }

    mesh->incrementNumberOfVertices(-1);
    mesh->incrementNumberOfEdges(-1);
  }


  void CoarseningManager2d::coarsenFunction(ElInfo* el_info)
  {
    Triangle* el = dynamic_cast<Triangle*>(const_cast<Element*>(el_info->getElement()));
    DegreeOfFreedom* edge[2];
    int n_neigh, bound = 0;
    RCNeighbourList coarse_list(2);

    coarse_list.setCoarseningManager(this);

    if (el->getMark() >= 0)
      return; // el must not be coarsend, return
    if (!(el->getChild(0)))
      return;  // single leaves don't get coarsened

    if (el->getChild(0)->getMark() >= 0  || el->getChild(1)->getMark() >= 0)
    {
      // one of the children must not be coarsend; return :(
      el->setMark(0);
      return;
    }

    if (!el->getChild(0)->isLeaf() || !el->getChild(1)->isLeaf())
    {
      // one of the children is not a leaf element; try again later on
      doMore = true;
      return;
    }

    // give the refinement edge the right orientation

    if (el->getDof(0,0) < el->getDof(1,0))
    {
      edge[0] = const_cast<DegreeOfFreedom*>(el->getDof(0));
      edge[1] = const_cast<DegreeOfFreedom*>(el->getDof(1));
    }
    else
    {
      edge[1] = const_cast<DegreeOfFreedom*>(el->getDof(0));
      edge[0] = const_cast<DegreeOfFreedom*>(el->getDof(1));
    }

    coarse_list.setElement(0, el, true);

    n_neigh = 1;
    if (coarse_list.setElement(1, el_info->getNeighbour(2)))
    {
      n_neigh = 2;
      coarse_list.setCoarsePatch(1, el_info->getOppVertex(2) == 2);
    }

    // check wether we can coarsen the patch or not

    // ==========================================================================
    // === check for periodic boundary ==========================================
    // ==========================================================================

    if (coarse_list.doCoarsePatch(n_neigh))
    {
      int n_neigh_periodic;
      DegreeOfFreedom* next_edge[2];
      RCNeighbourList periodicList;

      while (edge[0] != NULL)
      {
        coarse_list.periodicSplit(edge, next_edge,
                                  &n_neigh, &n_neigh_periodic,
                                  periodicList);

        coarsenPatch(periodicList, n_neigh_periodic, bound);

        edge[0] = next_edge[0];
        edge[1] = next_edge[1];
      }
    }
  }

} // end namespace AMDiS
