#include "CoarseningManager.hpp"

// AMDiS includes
#include "AdaptInstationary.hpp"
#include "AdaptStationary.hpp"
#include "FixVec.hpp"
#include "MacroElement.hpp"
#include "Mesh.hpp"
#include "ProblemStatBase.hpp"   // => MESH_COARSENED
#include "RCNeighbourList.hpp"
#include "Traverse.hpp"

namespace AMDiS
{

  Flag CoarseningManager::globalCoarsen(Mesh* aMesh, int mark)
  {
    if (mark >= 0)
      return 0;

    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(aMesh, -1, Mesh::CALL_LEAF_EL);
    while (elInfo)
    {
      elInfo->getElement()->setMark(mark);
      elInfo = stack.traverseNext(elInfo);
    }

    return coarsenMesh(aMesh);
  }

  void CoarseningManager::spreadCoarsenMark()
  {
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_POSTORDER);
    while (elInfo)
    {
      Element* el = elInfo->getElement();

      if (el->getChild(0))
      {
        // interior node of the tree
        int mark = std::max(el->getChild(0)->getMark(), el->getChild(1)->getMark());
        el->setMark(std::min(mark + 1, 0));
      }
      else
      {
        // leaf node of the tree
        if (el->getMark() < 0)
          el->setMark(el->getMark() - 1);
      }

      elInfo = stack.traverseNext(elInfo);
    }
  }

  void CoarseningManager::cleanUpAfterCoarsen()
  {
    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
    while (elInfo)
    {
      Element* el = elInfo->getElement();
      el->setMark(std::max(el->getMark(), 0));
      elInfo = stack.traverseNext(elInfo);
    }
  }

  Flag CoarseningManager::coarsenMesh(Mesh* aMesh)
  {
    mesh = aMesh;
    int nElements = mesh->getNumberOfLeaves();

    spreadCoarsenMark();

    stack = new TraverseStack;

    do
    {
      doMore = false;
      ElInfo* el_info =
        stack->traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_POSTORDER | Mesh::FILL_NEIGH);
      while (el_info)
      {
        coarsenFunction(el_info);
        el_info = stack->traverseNext(el_info);
      }
    }
    while (doMore);

    delete stack;

    cleanUpAfterCoarsen();

    nElements -= mesh->getNumberOfLeaves();

    if (nElements != 0)
    {
      aMesh->incChangeIndex();
      return MESH_COARSENED;
    }
    else
    {
      return Flag(0);
    }
  }

} // end namespace AMDiS
