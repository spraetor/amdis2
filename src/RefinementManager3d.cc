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
#include "DOFVector.h"
#include "PeriodicBC.h"
#include "VertexVector.h"
#include "Debug.h"

using namespace std;

namespace AMDiS 
{
  map<Mesh*, FixRefinementPatch::ConnectedEdges> FixRefinementPatch::connectedEdges;

  void RefinementManager3d::bisectTetrahedron(RCNeighbourList& refineList, 
                              					      int index,
                              					      DegreeOfFreedom* dof[3],
                              					      DegreeOfFreedom *edge[2])
  {
    Tetrahedron *el = 
      dynamic_cast<Tetrahedron*>(const_cast<Element*>(refineList.getElement(index)));
    Tetrahedron *child[2];
    int el_type = refineList.getType(index);

    child[0] = dynamic_cast<Tetrahedron*>(mesh->createNewElement(el));
    child[1] = dynamic_cast<Tetrahedron*>(mesh->createNewElement(el));
  
    int mark = std::max(0, el->getMark() - 1);
    child[0]->setMark(mark);
    child[1]->setMark(mark);
    el->setMark(0);

    /****************************************************************************/
    /*  transfer hidden data from parent to children                            */
    /****************************************************************************/

    el->refineElementData(child[0], child[1], el_type);

    el->setFirstChild(child[0]);
    el->setSecondChild(child[1]);

    if (child[0]->getMark() > 0)
      doMoreRecursiveRefine = true;

    int n_vertices = mesh->getGeo(VERTEX);
    child[0]->setDof(n_vertices - 1, dof[0]);
    child[1]->setDof(n_vertices - 1, dof[0]);

    for (int i = 0; i < n_vertices - 1; i++) {
      child[0]->
	setDof(i, const_cast<DegreeOfFreedom*>(el->getDof(Tetrahedron::childVertex[el_type][0][i])));
      child[1]->
	setDof(i, const_cast<DegreeOfFreedom*>(el->getDof(Tetrahedron::childVertex[el_type][1][i])));
    }
    /****************************************************************************/
    /*  there is one more leaf element and two more hierachical elements        */
    /****************************************************************************/

    mesh->incrementNumberOfLeaves(1);
    mesh->incrementNumberOfElements(2);

    /****************************************************************************/
    /* first set those dof pointers for higher order without neighbour          */
    /* information                                                              */
    /****************************************************************************/

    if (mesh->getNumberOfDofs(EDGE)) {
      int node = mesh->getNode(EDGE);

      /****************************************************************************/
      /*  set pointers to those dof's that are handed on from the parant          */
      /****************************************************************************/
      
      child[0]->
	setDof(node, 
	       const_cast<DegreeOfFreedom*>(el->getDof(node + Tetrahedron::childEdge[el_type][0][0])));
      child[1]->
	setDof(node, 
	       const_cast<DegreeOfFreedom*>(el->getDof(node + Tetrahedron::childEdge[el_type][1][0])));
      child[0]->
	setDof(node + 1, 
	       const_cast<DegreeOfFreedom*>(el->getDof(node + Tetrahedron::childEdge[el_type][0][1])));
      child[1]->
	setDof(node + 1, 
	       const_cast<DegreeOfFreedom*>(el->getDof(node + Tetrahedron::childEdge[el_type][1][1])));
      child[0]->
	setDof(node + 3, 
	       const_cast<DegreeOfFreedom*>(el->getDof(node + Tetrahedron::childEdge[el_type][0][3])));
      child[1]->
	setDof(node + 3, 
	       const_cast<DegreeOfFreedom*>(el->getDof(node + Tetrahedron::childEdge[el_type][1][3])));
      
      /****************************************************************************/
      /*  adjust pointers to the dof's in the refinement edge                     */
      /****************************************************************************/
      
      if (el->getDof(0) == edge[0]) {
	child[0]->setDof(node + 2, dof[1]);
	child[1]->setDof(node + 2, dof[2]);
      } else {
	child[0]->setDof(node + 2, dof[2]);
	child[1]->setDof(node + 2, dof[1]);
      }
    }
    
    if (mesh->getNumberOfDofs(FACE)) {
      int node = mesh->getNode(FACE);    
      
      /****************************************************************************/
      /*  set pointers to those dof's that are handed on from the parant          */
      /****************************************************************************/
      
      child[0]->setDof(node + 3, const_cast<DegreeOfFreedom*>(el->getDof(node + 1)));
      child[1]->setDof(node + 3, const_cast<DegreeOfFreedom*>(el->getDof(node + 0)));
      
      /****************************************************************************/
      /*  get new dof for the common face of child0 and child1                    */
      /****************************************************************************/
      
      DegreeOfFreedom *newDOF = mesh->getDof(FACE);
      child[0]->setDof(node, static_cast<DegreeOfFreedom*>(newDOF));
      child[1]->setDof(node, static_cast<DegreeOfFreedom*>(newDOF));
    }
    
    if (mesh->getNumberOfDofs(CENTER)) {
      int node = mesh->getNode(CENTER);
      child[0]->setDof(node, const_cast<DegreeOfFreedom*>(mesh->getDof(CENTER)));
      child[1]->setDof(node, const_cast<DegreeOfFreedom*>(mesh->getDof(CENTER)));
    }

    if (mesh->getNumberOfDofs(EDGE) || mesh->getNumberOfDofs(FACE))
      fillPatchConnectivity(refineList, index);
  }


  void RefinementManager3d::fillPatchConnectivity(RCNeighbourList &refineList,
						  int index)
  {
    FUNCNAME_DBG("RefinementManager3d::fillPatchConnectivity");

    Element *el = refineList.getElement(index);
    int el_type = refineList.getType(index);
    int n_type = 0;
    int dir, adjc, i_neigh, j_neigh;
    int node0, node1, oppVertex = 0;

    for (dir = 0; dir < 2; dir++) {
      Element *neigh = refineList.getNeighbourElement(index, dir);
      if (neigh) {
	n_type = refineList.getType(refineList.getNeighbourNr(index, dir));
	oppVertex = refineList.getOppVertex(index, dir);
      }

      if (!neigh || neigh->isLeaf()) {
	/****************************************************************************/
	/*  get new dof's in the midedge of the face of el and for the two midpoints*/
	/*  of the sub-faces. If face is an interior face those pointers have to be */
	/*  adjusted by the neighbour element also (see below)                      */
	/****************************************************************************/

	if (mesh->getNumberOfDofs(EDGE)) {
	  node0 = node1 = mesh->getNode(EDGE);
	  node0 += Tetrahedron::nChildEdge[el_type][0][dir];
	  node1 += Tetrahedron::nChildEdge[el_type][1][dir];
	  DegreeOfFreedom *newDOF = mesh->getDof(EDGE);
	  (const_cast<Element*>(el->getFirstChild()))->setDof(node0, newDOF);
	  (const_cast<Element*>(el->getSecondChild()))->setDof(node1, newDOF);
	}
	if (mesh->getNumberOfDofs(FACE)) {
	  node0 = mesh->getNode(FACE) + Tetrahedron::nChildFace[el_type][0][dir];
	  (const_cast<Element*>(el->getFirstChild()))->setDof(node0, mesh->getDof(FACE));
	  node1 = mesh->getNode(FACE) + Tetrahedron::nChildFace[el_type][1][dir];
	  (const_cast<Element*>(el->getSecondChild()))->setDof(node1, mesh->getDof(FACE));
	}
      } else {
	/****************************************************************************/
	/*  interior face and neighbour has been refined, look for position at the  */
	/*  refinement edge                                                         */
	/****************************************************************************/
      
	if (el->getDof(0) == neigh->getDof(0)) {
	  // Same position at refinement edge.
	  adjc = 0;
	} else {
	  // Different position at refinement edge.
	  adjc = 1;
	}

	for (int i = 0; i < 2; i++) {
	  int j = Tetrahedron::adjacentChild[adjc][i];

	  i_neigh = Tetrahedron::nChildFace[el_type][i][dir];
	  j_neigh = Tetrahedron::nChildFace[n_type][j][oppVertex - 2];

	  /****************************************************************************/
	  /*  adjust dof pointer in the edge in the common face of el and neigh and   */
	  /*  the dof pointer in the sub-face child_i-child_j (allocated by neigh!)   */
	  /****************************************************************************/

	  if (mesh->getNumberOfDofs(EDGE)) {
	    node0 = mesh->getNode(EDGE) + Tetrahedron::nChildEdge[el_type][i][dir];
	    node1 = mesh->getNode(EDGE) + Tetrahedron::nChildEdge[n_type][j][oppVertex - 2];

	    TEST_EXIT_DBG(neigh->getChild(j)->getDof(node1))
	      ("no dof on neighbour %d at node %d\n", 
	       neigh->getChild(j)->getIndex(), node1);

	    (const_cast<Element*>(el->getChild(i)))->
	      setDof(node0, const_cast<DegreeOfFreedom*>(neigh->getChild(j)->getDof(node1)));
	  }
	  if (mesh->getNumberOfDofs(FACE)) {
	    node0 = mesh->getNode(FACE) + i_neigh;
	    node1 = mesh->getNode(FACE) + j_neigh;

	    TEST_EXIT_DBG(neigh->getChild(j)->getDof(node1))
	      ("No DOF on neighbour %d at node %d!\n",
	       neigh->getChild(j)->getIndex(), node1);

	    (const_cast<Element*>(el->getChild(i)))->
	      setDof(node0, const_cast<DegreeOfFreedom*>(neigh->getChild(j)->getDof(node1)));
	  }

	}  /*   for (i = 0; i < 2; i++)                                       */
      }    /*   else of   if (!neigh  ||  !neigh->child[0])                   */
    }      /*   for (dir = 0; dir < 2; dir++)                                 */
  }


  void RefinementManager3d::newCoordsFct(ElInfo *elInfo, RCNeighbourList &refineList)
  {
    FUNCNAME("RefinementManager3d::newCoordsFct()");

    Element *el = elInfo->getElement();
    DegreeOfFreedom *edge[2];
    ElInfo *elinfo = elInfo;
    int dow = Global::getGeo(WORLD);
    Projection *projector = elInfo->getProjection(0);

    if (!projector || projector->getType() != VOLUME_PROJECTION)
      projector = elInfo->getProjection(4);    

    if (el->getFirstChild() && projector && (!el->isNewCoordSet())) {
      WorldVector<double> *new_coord = new WorldVector<double>;

      for (int j = 0; j < dow; j++)
	(*new_coord)[j] = (elInfo->getCoord(0)[j] + elInfo->getCoord(1)[j]) * 0.5;

      projector->project(*new_coord);

      el->setNewCoord(new_coord);
      /****************************************************************************/
      /*  now, information should be passed on to patch neighbours...             */
      /*  get the refinement patch                                                */
      /****************************************************************************/
      refineList.setElement(0, el);
      refineList.setElType(0, elInfo->getType());
      int n_neigh = 1;

      for (int i = 0; i < 2; i++)
	edge[i] = const_cast<DegreeOfFreedom*>(elInfo->getElement()->getDof(i));

      if (getRefinePatch(&elinfo, edge, 0, refineList, &n_neigh)) {
	// Domain's boundary was reached while looping around the refinement edge.

	getRefinePatch(&elinfo, edge, 1, refineList, &n_neigh);
      }

      for (int i = 1; i < n_neigh; i++) {            /* start with 1, as list[0]=el */
	TEST(!refineList.getElement(i)->isNewCoordSet())
	  ("non-nil new_coord in el %d refineList[%d] el %d (n_neigh=%d)\n",
	   el->getIndex(), i, refineList.getElement(i)->getIndex(), n_neigh);
	
	refineList.getElement(i)->setNewCoord(el->getNewCoord());
      }
    }
  }


  void RefinementManager3d::setNewCoords(int macroEl)
  {
    RCNeighbourList refineList(mesh->getMaxEdgeNeigh());
    ElInfo *elInfo;

    if (macroEl == -1)
      elInfo = stack->traverseFirst(mesh, -1, 
				    Mesh::CALL_EVERY_EL_PREORDER | 
				    Mesh::FILL_BOUND | Mesh::FILL_COORDS | 
				    Mesh::FILL_NEIGH);
    else
      elInfo = stack->traverseFirstOneMacro(mesh, macroEl, -1,
					    Mesh::CALL_EVERY_EL_PREORDER | 
					    Mesh::FILL_BOUND | Mesh::FILL_COORDS | 
					    Mesh::FILL_NEIGH);
    

    while (elInfo) {
      newCoordsFct(elInfo, refineList);
      elInfo = stack->traverseNext(elInfo);
    }
  }


  DegreeOfFreedom RefinementManager3d::refinePatch(DegreeOfFreedom *edge[2], 
						   RCNeighbourList &refineList,
						   int n_neigh, bool bound)
  {
    Tetrahedron *el = 
      dynamic_cast<Tetrahedron*>(const_cast<Element*>(refineList.getElement(0)));
    /* first element in the list */
    DegreeOfFreedom *dof[3] = {NULL, NULL, NULL};

    /****************************************************************************/
    /*  get new dof's in the refinement edge                                    */
    /****************************************************************************/

    dof[0] = mesh->getDof(VERTEX);
    mesh->incrementNumberOfVertices(1);
  
    if (mesh->getNumberOfDofs(EDGE)) {
      dof[1] = mesh->getDof(EDGE);
      dof[2] = mesh->getDof(EDGE);
    }

    for (int i = 0; i < n_neigh; i++)
      bisectTetrahedron(refineList, i, dof, edge);

    /****************************************************************************/
    /*  if there are functions to interpolate data to the finer grid, do so     */
    /****************************************************************************/
    for (int iadmin = 0; iadmin < mesh->getNumberOfDOFAdmin(); iadmin++) {
      std::list<DOFIndexedBase*>::iterator it;
      DOFAdmin* admin = const_cast<DOFAdmin*>(&mesh->getDofAdmin(iadmin));
      std::list<DOFIndexedBase*>::iterator end = admin->endDOFIndexed();
      for (it = admin->beginDOFIndexed(); it != end; it++)
	(*it)->refineInterpol(refineList, n_neigh);
    }

    if (!mesh->queryCoarseDOFs()) {
      /****************************************************************************/
      /*  if there should be no dof information on interior leaf elements remove  */
      /*  dofs from edges, faces and the centers of parents                       */
      /****************************************************************************/
      if (mesh->getNumberOfDofs(EDGE)) {
	/****************************************************************************/
	/*  remove dof of the midpoint of the common refinement edge                */
	/****************************************************************************/

	el = dynamic_cast<Tetrahedron*>(const_cast<Element*>(refineList.getElement(0)));
	mesh->freeDof(const_cast<DegreeOfFreedom*>(el->getDof(mesh->getNode(EDGE))), EDGE);
      }
      
      if (mesh->getNumberOfDofs(EDGE) || 
	  mesh->getNumberOfDofs(FACE) ||
	  mesh->getNumberOfDofs(CENTER)) {
	for (int i = 0; i < n_neigh; i++)
	  refineList.removeDOFParent(i);
      }
    }


    /****************************************************************************/
    /*  update the number of edges and faces; depends whether refinement edge   */
    /*  is a boundary edge or not                                               */
    /****************************************************************************/

    if (bound) {
      mesh->incrementNumberOfEdges(n_neigh + 2);
      mesh->incrementNumberOfFaces(2 * n_neigh + 1);
      newCoords = true; // added to allow BOUNDARY_PROJECTION
    } else {
      mesh->incrementNumberOfEdges(n_neigh + 1);
      mesh->incrementNumberOfFaces(2 * n_neigh);
    }
    
    return dof[0][0];
  }


  bool RefinementManager3d::getRefinePatch(ElInfo **elInfo, 
					   DegreeOfFreedom *edge[2], 
					   int direction,
					   RCNeighbourList &refineList, 
					   int *n_neigh)
  {
    FUNCNAME_DBG("RefinementManager3d::getRefinePatch()");

    int localNeighbour = 3 - direction;
    Tetrahedron *el = 
      dynamic_cast<Tetrahedron*>(const_cast<Element*>((*elInfo)->getElement()));

    if ((*elInfo)->getNeighbour(localNeighbour) == NULL)
      return true;    
  
    int oppVertex = (*elInfo)->getOppVertex(localNeighbour);
#if DEBUG
    int testIndex = (*elInfo)->getNeighbour(localNeighbour)->getIndex();
#endif
    ElInfo *neighInfo = stack->traverseNeighbour3d((*elInfo), localNeighbour);
    int neighElType = neighInfo->getType();

    TEST_EXIT_DBG(neighInfo->getElement()->getIndex() == testIndex)
      ("Should not happen!\n");

    Tetrahedron *neigh = 
      dynamic_cast<Tetrahedron*>(const_cast<Element*>(neighInfo->getElement())); 
    int vertices = mesh->getGeo(VERTEX);

    while (neigh != el) {
      // === Determine the common edge of the element and its neighbour. ===

      int edgeDof0, edgeDof1;

      // get local/elementwise DOF indices of Start and End Vertices of EDGE
      // on Neighbour Element
      for (edgeDof0 = 0; edgeDof0 < vertices; edgeDof0++)
	if (neigh->getDof(edgeDof0) == edge[0])
	  break;
      for (edgeDof1 = 0; edgeDof1 < vertices; edgeDof1++)
	if (neigh->getDof(edgeDof1) == edge[1])
	  break;

      if (edgeDof0 > 3 || edgeDof1 > 3) {
	for (edgeDof0 = 0; edgeDof0 < vertices; edgeDof0++)
	  if (mesh->associated(neigh->getDof(edgeDof0, 0), edge[0][0]))  
	    break;
	for (edgeDof1 = 0; edgeDof1 < vertices; edgeDof1++)
	  if (mesh->associated(neigh->getDof(edgeDof1, 0), edge[1][0]))  
	    break;
	    
	if (edgeDof0 > 3 || edgeDof1 > 3) {
	  for (edgeDof0 = 0; edgeDof0 < vertices; edgeDof0++)
	    if (mesh->indirectlyAssociated(neigh->getDof(edgeDof0, 0), edge[0][0]))  
	      break;
	  for (edgeDof1 = 0; edgeDof1 < vertices; edgeDof1++)
	    if (mesh->indirectlyAssociated(neigh->getDof(edgeDof1, 0), edge[1][0]))  
	      break;
	    
	  TEST_EXIT_DBG(edgeDof0 < vertices)
 	    ("DOF %d not found on element %d with nodes (%d %d %d %d)\n", 
 	     edge[0][0], neigh->getIndex(), neigh->getDof(0, 0),
 	     neigh->getDof(1, 0), neigh->getDof(2, 0), neigh->getDof(3, 0));

	  TEST_EXIT_DBG(edgeDof1 < vertices)
 	    ("DOF %d not found on element %d with nodes (%d %d %d %d)\n", 
 	     edge[1][0], neigh->getIndex(), neigh->getDof(0, 0),
 	     neigh->getDof(1, 0), neigh->getDof(2, 0), neigh->getDof(3, 0));
	}
      }

      TEST_EXIT_DBG(edgeDof0 < vertices && edgeDof1 < vertices)
	("DOF %d or DOF %d not found on element %d with nodes (%d %d %d %d)\n", 
	 edge[0][0], edge[1][0], neigh->getIndex(), neigh->getDof(0, 0),
	 neigh->getDof(1, 0), neigh->getDof(2, 0), neigh->getDof(3, 0));

      int edgeNo = Tetrahedron::edgeOfDofs[edgeDof0][edgeDof1];

      if (edgeNo) {
	// Only 0 can be a compatible commen refinement edge. Thus, neigh has not 
	// a compatible refinement edge. Refine it first.

	neigh->setMark(std::max(neigh->getMark(), 1));

	neighInfo = refineFunction(neighInfo);

	// === Now, go to a child at the edge and determine the opposite vertex ===
	// === for  this child; continue the looping around the edge with this  ===
	// === element.                                                         ===
	
	neighInfo = stack->traverseNext(neighInfo);
	neighElType = neighInfo->getType();
	bool reverseMode = stack->getTraverseFlag().isSet(Mesh::CALL_REVERSE_MODE);
	
	switch (edgeNo) {
	case 1: 
	  if (reverseMode) {
	    neighInfo = stack->traverseNext(neighInfo);
	    neighElType = neighInfo->getType();
	  }
	    
	  oppVertex = (oppVertex == 1 ? 3 : 2);
	  break;
	case 2: 
	  if (reverseMode) {
	    neighInfo = stack->traverseNext(neighInfo);
	    neighElType = neighInfo->getType();
	  }

	  oppVertex = (oppVertex == 2 ? 1 : 3);
	  break;
	case 3: 
	  if (!reverseMode) {
	    neighInfo = stack->traverseNext(neighInfo);
	    neighElType = neighInfo->getType();
	  }

	  if (neighElType != 1)
	    oppVertex = (oppVertex == 0 ? 3 : 2);
	  else
	    oppVertex = (oppVertex == 0 ? 3 : 1);
	  break;
	case 4:
	  if (!reverseMode) {
	    neighInfo = stack->traverseNext(neighInfo);
	    neighElType = neighInfo->getType();
	  }

	  if (neighElType != 1)
	    oppVertex = (oppVertex == 0 ? 3 : 1);
	  else
	    oppVertex = (oppVertex == 0 ? 3 : 2);
	  break;
	case 5:
	  if (neighElType != 1) {
	    if (!reverseMode) {
	      neighInfo = stack->traverseNext(neighInfo);
	      neighElType = neighInfo->getType();
	    }
	  }
	  oppVertex = 3;
	  break;
	}

	neigh = 
	  dynamic_cast<Tetrahedron*>(const_cast<Element*>(neighInfo->getElement()));
      } else {
	// Neigh is compatible devisible. Put neigh to the list of patch elements 
	// and go to next neighbour.

	TEST_EXIT_DBG(*n_neigh < mesh->getMaxEdgeNeigh())
	  ("too many neighbours %d in refine patch\n", *n_neigh);
		
	refineList.setElement(*n_neigh, neigh);
	refineList.setElType(*n_neigh, neighElType);
	
	// We have to go back to the starting element via oppVertex values.

	refineList.setOppVertex(*n_neigh, 0, oppVertex); 
	
	(*n_neigh)++;

	int i = (oppVertex != 3 ? 3 : 2);

	if (neighInfo->getNeighbour(i)) {
	  oppVertex = neighInfo->getOppVertex(i);
    
#if DEBUG
	  int testIndex = neighInfo->getNeighbour(i)->getIndex();
#endif

	  neighInfo = stack->traverseNeighbour3d(neighInfo, i);

	  TEST_EXIT_DBG(neighInfo->getElement()->getIndex() == testIndex)
	    ("Should not happen!\n");

	  neighElType = neighInfo->getType();
	  neigh = 
	    dynamic_cast<Tetrahedron*>(const_cast<Element*>(neighInfo->getElement()));
	} else {
	  break;
	}
      }
    }
   
    if (neigh == el) {
      (*elInfo) = neighInfo;

      return false;
    }

   
    // The domain's boundary is reached. Loop back to the starting el.
    
    int i = *n_neigh - 1;
    oppVertex = refineList.getOppVertex(i, 0);
    do {
      TEST_EXIT_DBG(neighInfo->getNeighbour(oppVertex) && i > 0)
	("While looping back domains boundary was reached or i == 0\n");
      oppVertex = refineList.getOppVertex(i--, 0);

#if DEBUG
      int testIndex = neighInfo->getNeighbour(oppVertex)->getIndex();
#endif
      neighInfo = stack->traverseNeighbour3d(neighInfo, oppVertex);

      int edgeDof0, edgeDof1;
      for (edgeDof0 = 0; edgeDof0 < vertices; edgeDof0++)
	if (neigh->getDof(edgeDof0) == edge[0])
	  break;
      for (edgeDof1 = 0; edgeDof1 < vertices; edgeDof1++)
	if (neigh->getDof(edgeDof1) == edge[1])
	  break;

      TEST_EXIT_DBG(neighInfo->getElement()->getIndex() == testIndex)
	("Should not happen!\n");
    } while (neighInfo->getElement() != el);

    (*elInfo) = neighInfo;
    
    return true;
  }


  ElInfo* RefinementManager3d::refineFunction(ElInfo* elInfo)
  {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    FUNCNAME_DBG("RefinementManager3d::refineFunction()");
#endif

    Element *el = elInfo->getElement();

    if (el->getMark() <= 0)  
      return elInfo;

    int bound = false;
    DegreeOfFreedom *edge[2];

    // === Get memory for a list of all elements at the refinement edge. ===

    RCNeighbourList refineList(mesh->getMaxEdgeNeigh());
    refineList.setElType(0, elInfo->getType());
    refineList.setElement(0, el);
    int n_neigh = 1;

    if (elInfo->getProjection(0) && 
	elInfo->getProjection(0)->getType() == VOLUME_PROJECTION)
      newCoords = true;


    // === Give the refinement edge the right orientation. ===

    if (el->getDof(0, 0) < el->getDof(1, 0)) {
      edge[0] = const_cast<DegreeOfFreedom*>(el->getDof(0));
      edge[1] = const_cast<DegreeOfFreedom*>(el->getDof(1));
    } else {
      edge[1] = const_cast<DegreeOfFreedom*>(el->getDof(0));
      edge[0] = const_cast<DegreeOfFreedom*>(el->getDof(1));
    }


#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    vector<FixRefinementPatch::EdgeInEl> refineEdges;
    FixRefinementPatch::getOtherEl(mesh, stack, refineEdges);
//     if (refineEdges.size()) {
//       MSG("FIX REFINEMENT PATH ON ELEMENT %d %d: %d additional edges\n", 
// 	  elInfo->getElement()->getIndex(),
// 	  elInfo->getMacroElement()->getIndex(),
// 	  refineEdges.size());
//     }
#endif

    // === Traverse and refine the refinement patch. ====

    if (getRefinePatch(&elInfo, edge, 0, refineList, &n_neigh)) {      
      // Domain's boundary was reached while looping around the refinement edge
      getRefinePatch(&elInfo, edge, 1, refineList, &n_neigh);
      bound = true;
    }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    // === If the refinement edge must be fixed, add also the other part  ===
    // === of this edge to the refinement patch.                          ===

    for (int edgeIndex = 0; 
 	 edgeIndex < static_cast<int>(refineEdges.size()); edgeIndex++) {
      Element *otherEl = refineEdges[edgeIndex].first;   
      TraverseStack stack2;
      ElInfo *elInfo2 = 
	stack2.traverseFirstOneMacro(mesh, otherEl->getIndex(), -1, 
				     Mesh::CALL_LEAF_EL | 
				     Mesh::FILL_NEIGH | 
				     Mesh::FILL_BOUND);

      bool foundEdge = false;

      while (elInfo2) {
	for (int i = 0; i < 6; i++) {
	  DofEdge edge2 = elInfo2->getElement()->getEdge(i);
	  if (edge2.first == *(edge[0]) &&
	      edge2.second == *(edge[1])) {

	    // We found the edge on the other element side on leaf level. Two
	    // cases may occure: either it is already a refinement edge, i.e,
	    // it has local edge number 0, or it is not a refinement edge. In
	    // the first case, we are fine and this leaf element may be set
	    // to all elements on the refinement edge. Otherwise, the element
	    // must be refined at least once to get a refinement edge.

	    if (i == 0) {
	      // Edge is refinement edge, so add it to refine list.

	      refineList.setElType(n_neigh, elInfo2->getType());
	      refineList.setElement(n_neigh, elInfo2->getElement());
	      n_neigh++;
	      
	      TraverseStack *tmpStack = stack;
	      stack = &stack2;
	      
	      if (getRefinePatch(&elInfo2, edge, 0, refineList, &n_neigh)) {
		getRefinePatch(&elInfo2, edge, 1, refineList, &n_neigh);
		bound = true;
	      }
	      
	      stack = tmpStack;
	      foundEdge = true;
	      break;
	    } else {
	      // Edge i not refinement edge, so refine the element further.

	      Element *el2 = elInfo2->getElement();
	      el2->setMark(std::max(el2->getMark(), 1));

	      TraverseStack *tmpStack = stack;
	      stack = &stack2;

	      elInfo2 = refineFunction(elInfo2);

	      stack = tmpStack;
	      break;
	    }
	  }
	}

	if (foundEdge)
	  break;

	elInfo2 = stack2.traverseNext(elInfo2);
      }

      TEST_EXIT_DBG(foundEdge)("Should not happen!\n");
    }
#endif

    // fill neighbour information inside the patch in the refinement list
    refineList.fillNeighbourRelations(n_neigh, bound);

    // ============ Check for periodic boundary ============

    DegreeOfFreedom *next_edge[2] = {NULL, NULL};
    DegreeOfFreedom *first_edge[2] = {edge[0], edge[1]};
    DegreeOfFreedom *last_edge[2] = {NULL, NULL};
    int n_neigh_periodic = 0;

    DegreeOfFreedom lastNewDof = -1;
    DegreeOfFreedom firstNewDof = -1;

    RCNeighbourList periodicList;

    while (edge[0] != NULL) {
      refineList.periodicSplit(edge, next_edge, 
			       &n_neigh, &n_neigh_periodic, 
			       periodicList);

      DegreeOfFreedom newDof = 
	refinePatch(edge, periodicList, n_neigh_periodic, bound);

      if (firstNewDof == -1)
	firstNewDof = newDof;

      if (lastNewDof != -1) {
	for (std::map<int, VertexVector*>::iterator it = mesh->getPeriodicAssociations().begin();
	     it != mesh->getPeriodicAssociations().end(); ++it) {
	  if (it->second) {
	    if (((*(it->second))[edge[0][0]] == last_edge[0][0] &&
		 (*(it->second))[edge[1][0]] == last_edge[1][0]) || 
		((*(it->second))[edge[0][0]] == last_edge[1][0] &&
		 (*(it->second))[edge[1][0]] == last_edge[0][0])) {
	      (*(it->second))[lastNewDof] = newDof;
	      (*(it->second))[newDof] = lastNewDof;
	    } 
	  }
	}
      }
      lastNewDof = newDof;

      last_edge[0] = edge[0];
      last_edge[1] = edge[1];

      edge[0] = next_edge[0];
      edge[1] = next_edge[1];
    }

    if (lastNewDof != firstNewDof) {
      for (std::map<int, VertexVector*>::iterator it = mesh->getPeriodicAssociations().begin();
	   it != mesh->getPeriodicAssociations().end(); ++it) {
	if (it->second) {
	  if (((*(it->second))[first_edge[0][0]] == last_edge[0][0] &&
	       (*(it->second))[first_edge[1][0]] == last_edge[1][0]) || 
	      ((*(it->second))[first_edge[0][0]] == last_edge[1][0] &&
	       (*(it->second))[first_edge[1][0]] == last_edge[0][0])) {	    
	    (*(it->second))[lastNewDof] = firstNewDof;
	    (*(it->second))[firstNewDof] = lastNewDof;
	  }
	}
      }
    }

    stack->update();

    return elInfo;
  }


  void FixRefinementPatch::getOtherEl(Mesh* mesh,
				      TraverseStack *stack, 
				      vector<EdgeInEl>& refineEdges)
  {
    FUNCNAME_DBG("FixRefinementPatch::getOtherEl()");

    if (!FixRefinementPatch::connectedEdges[mesh].empty()) {
      // === Get stack of current traverse. ===
      vector<ElInfo*> elInfos;
      vector<int> infos;      
      int stackUsed = stack->getStackData(elInfos, infos);         
      int checkIndex = stackUsed;
      int localEdgeNo = 0;
      ElInfo *elInfo = elInfos[stackUsed];
      
      TEST_EXIT_DBG(elInfo->getMesh() == mesh)("Something is wrong.\n");

      // === Calculate the refinement edge number on the macro element level. ===
      
      for (int i = 0; i < elInfo->getLevel(); i++) {
	TEST_EXIT_DBG(checkIndex >= 1)("Should not happen!\n");
	

	int parentType = elInfos[checkIndex - 1]->getType();
	int ithParentChild = 0;
	if (elInfos[checkIndex - 1]->getElement()->getChild(1) ==
	    elInfos[checkIndex]->getElement())
	  ithParentChild = 1;

	// If local edge number if equal or greater to 4, its a new edge and
	// thus cannot be part of the parent edge.
 	if (localEdgeNo >= 4) {
 	  localEdgeNo = -1;
 	  break;
 	}
	
	localEdgeNo = 
	  Tetrahedron::childEdge[parentType][ithParentChild][localEdgeNo];

	checkIndex--;
      }
      
     
      // If the refinement edge is part of an edge on the macro level, we must 
      // check if the refinement edge is part of an edge that must be fixed.
      if (localEdgeNo >= 0) {
	int macroElIndex = elInfos[checkIndex]->getElement()->getIndex();

	TEST_EXIT_DBG(elInfos[checkIndex]->getLevel() == 0)
	  ("Should not happen!\n");
	TEST_EXIT_DBG(macroElIndex == elInfo->getMacroElement()->getIndex())
	  ("Should not happen!\n");
	TEST_EXIT_DBG(localEdgeNo <= 5)("Should not happen!\n");

	refineEdges.clear();

	for (int i = 0; i < static_cast<int>(connectedEdges[mesh].size()); i++) {
	  TEST_EXIT_DBG(connectedEdges[mesh][i].first.first->getMesh() == mesh)("Something is wrong.\n");
	  if (connectedEdges[mesh][i].first.first->getIndex() == macroElIndex &&
	      connectedEdges[mesh][i].first.second == localEdgeNo) {
	    // We have found that this edge must be fixed.
	    refineEdges.push_back(connectedEdges[mesh][i].second);
	  }
	}
      }
    }
  }
  
} // end namespace AMDiS
