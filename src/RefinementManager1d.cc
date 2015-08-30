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

namespace AMDiS
{
  void RefinementManager1d::recursiveRefineFunction(ElInfo* elInfo)
  {
    Line* el =
      dynamic_cast<Line*>(const_cast<Element*>(elInfo->getElement())), *child[2];

    Mesh* mesh = elInfo->getMesh();

    if (elInfo->getProjection(0))
      newCoord(true);

    if (el->getMark() <= 0)
      return;

    child[0] = dynamic_cast<Line*>(mesh->createNewElement(el));
    child[1] = dynamic_cast<Line*>(mesh->createNewElement(el));

    int mark = std::max(0, el->getMark() - 1);
    child[0]->setMark(mark);
    child[1]->setMark(mark);
    el->setMark(0);

    /*--------------------------------------------------------------------------*/
    /*  transfer hided data from parent to children                             */
    /*--------------------------------------------------------------------------*/

    el->refineElementData(child[0], child[1]);

    el->setFirstChild(child[0]);
    el->setSecondChild(child[1]);

    if (child[0]->getMark() > 0)
      doMoreRecursiveRefine = true;

    DegreeOfFreedom* newVertexDOFs = mesh->getDof(VERTEX);
    child[0]->setDof(1, newVertexDOFs);
    child[1]->setDof(0, newVertexDOFs);

    /*--------------------------------------------------------------------------*/
    /*  the other vertices are handed on from the parent                        */
    /*--------------------------------------------------------------------------*/
    child[0]->setDof(0, const_cast<DegreeOfFreedom*>(el->getDof(0)));
    child[1]->setDof(1, const_cast<DegreeOfFreedom*>(el->getDof(1)));

    /*--------------------------------------------------------------------------*/
    /*  there is one more leaf element, two hierachical elements,               */
    /*  one more vertex                                                         */
    /*--------------------------------------------------------------------------*/

    mesh->incrementNumberOfLeaves(1);
    mesh->incrementNumberOfVertices(1);
    mesh->incrementNumberOfElements(2);

    if (mesh->getNumberOfDofs(CENTER))
    {
      /*--------------------------------------------------------------------------*/
      /* there are dofs at the barycenter of the triangles                        */
      /*--------------------------------------------------------------------------*/
      child[0]->setDof(mesh->getNode(CENTER), const_cast<DegreeOfFreedom*>(mesh->getDof(CENTER)));
      child[1]->setDof(mesh->getNode(CENTER), const_cast<DegreeOfFreedom*>(mesh->getDof(CENTER)));
    }

    /*--------------------------------------------------------------------------*/
    /*  if there are functions to interpolate data to the finer grid, do so     */
    /*--------------------------------------------------------------------------*/

    RCNeighbourList  ref_list(1);  // = {{nil, 0, 0}};
    ref_list.setElement(0, el);

    int nrAdmin = mesh->getNumberOfDOFAdmin();
    for (int iadmin = 0; iadmin < nrAdmin; iadmin++)
    {
      std::list<DOFIndexedBase*>::iterator it;
      DOFAdmin* admin = const_cast<DOFAdmin*>(&mesh->getDofAdmin(iadmin));
      std::list<DOFIndexedBase*>::iterator end = admin->endDOFIndexed();
      for (it = admin->beginDOFIndexed(); it != end; it++)
        (*it)->refineInterpol(ref_list, 1);
    }

    if (!mesh->queryCoarseDOFs() && mesh->getNumberOfDofs(CENTER))
    {
      mesh->freeDof(const_cast<DegreeOfFreedom*>( el->getDof(mesh->getNode(CENTER))), CENTER);
      el->setDof(mesh->getNode(CENTER), NULL);
    }
  }

  Flag RefinementManager1d::refineMesh(Mesh* aMesh)
  {
    mesh = aMesh;
    int nElements = mesh->getNumberOfLeaves();

    doMoreRecursiveRefine = true;
    while (doMoreRecursiveRefine)
    {
      doMoreRecursiveRefine = false;

      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL | Mesh::FILL_BOUND | Mesh::FILL_COORDS);
      while (elInfo)
      {
        recursiveRefineFunction(elInfo);
        elInfo = stack.traverseNext(elInfo);
      }
    }

    nElements = mesh->getNumberOfLeaves() - nElements;

    if (newCoords)
      setNewCoords(); // call of sub-class method

    return (nElements ? MESH_REFINED : Flag(0));
  }

  void RefinementManager1d::newCoordsFct(ElInfo* elInfo)
  {
    Element* el = elInfo->getElement();
    int dow = Global::getGeo(WORLD);

    Projection* projector = elInfo->getProjection(0);

    if (el->getFirstChild() && projector && (!el->isNewCoordSet()))
    {
      WorldVector<double>* new_coord = new WorldVector<double>;

      for (int j = 0; j < dow; j++)
        (*new_coord)[j] = (elInfo->getCoord(0)[j] + elInfo->getCoord(1)[j]) * 0.5;

      projector->project(*new_coord);
      el->setNewCoord(new_coord);
    }
  }

  void RefinementManager1d::setNewCoords(int)
  {
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(mesh, -1,
                                         Mesh::CALL_EVERY_EL_PREORDER | Mesh::FILL_BOUND | Mesh::FILL_COORDS);
    while (elInfo)
    {
      newCoordsFct(elInfo);
      elInfo = stack.traverseNext(elInfo);
    }
  }

} // end namespace AMDiS
