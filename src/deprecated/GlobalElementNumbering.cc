#include "GlobalElementNumbering.h"
#include "MeshStructure.h"
#include "Mesh.h"
#include "Traverse.h"
#include "ElInfo.h"
#include "Element.h"

namespace AMDiS 
{
  GlobalElementNumbering::GlobalElementNumbering(MeshStructure *compositeStructure,
						 Mesh *localMesh)
  {
    int localIndex, globalIndex;
    compositeStructure->reset();
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(localMesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
    while(elInfo) {
      Element *element = elInfo->getElement();
      localIndex = element->getIndex();
      globalIndex = compositeStructure->getCurrentElement();
      localToGlobal_[localIndex] = globalIndex + 1;
      globalToLocal_[globalIndex] = localIndex + 1;
      if(element->isLeaf()) {
	compositeStructure->skipBranch();
      } else {
	compositeStructure->nextElement();
      }
      elInfo = stack.traverseNext(elInfo);
    }
  }

  int GlobalElementNumbering::getLocalIndex(int globalIndex)
  {
    return (globalToLocal_[globalIndex] - 1);
  }

  int GlobalElementNumbering::getGlobalIndex(int localIndex)
  {
    return (localToGlobal_[localIndex] - 1);
  }
  
} // end namespace AMDiS
