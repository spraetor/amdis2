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
#include "PeriodicBC.h"

namespace AMDiS 
{
  bool RefinementManager::doMoreRecursiveRefine = false;
  int RefinementManager::callRefineInterpol = 0;


  Flag RefinementManager::globalRefine(Mesh *aMesh, int mark)
  {
    if (mark <= 0) 
      return static_cast<Flag>(0);

    TraverseStack stack;
    ElInfo *elInfo = 
      stack.traverseFirst(aMesh, -1, 
			  Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS | Mesh::FILL_BOUND);
    while (elInfo) {
      elInfo->getElement()->setMark(mark);
      elInfo = stack.traverseNext(elInfo);
    }

    return refineMesh(aMesh);
  }


  Flag RefinementManager::refineMesh(Mesh *aMesh)
  {
    mesh = aMesh;
    int nElements = mesh->getNumberOfLeaves();
    newCoords = false;
    stack = new TraverseStack;
    doMoreRecursiveRefine = true;

    while (doMoreRecursiveRefine) {
      doMoreRecursiveRefine = false;
      ElInfo *elInfo = 
	stack->traverseFirst(mesh, -1, 
			     Mesh::CALL_LEAF_EL | Mesh::FILL_NEIGH | Mesh::FILL_BOUND);

      while (elInfo) {
	if (elInfo->getElement()->getMark() > 0) {	  
	  doMoreRecursiveRefine = 
	    doMoreRecursiveRefine || (elInfo->getElement()->getMark() > 1);
	  elInfo = refineFunction(elInfo);
	}	

	elInfo = stack->traverseNext(elInfo);
      }
    }
  
    if (newCoords)
      setNewCoords(); // call of sub-class method      

    delete stack;

    nElements -= mesh->getNumberOfLeaves();
 
    if (nElements != 0) {
      aMesh->incChangeIndex();
      return MESH_REFINED;
    } else {
      return Flag(0);
    }
  }


  void RefinementManager::refineMacroElement(Mesh *aMesh, int macroElIndex)
  {    
    mesh = aMesh;
    int nElements = mesh->getNumberOfLeaves();
    newCoords = false;
    doMoreRecursiveRefine = true;
    stack = new TraverseStack;
    
    while (doMoreRecursiveRefine) {
      doMoreRecursiveRefine = false;

      ElInfo *elInfo = 
	stack->traverseFirstOneMacro(mesh, macroElIndex, -1, 
				     Mesh::CALL_LEAF_EL | Mesh::FILL_NEIGH | Mesh::FILL_BOUND);
      
      while (elInfo) {
	if (elInfo->getElement()->getMark() > 0) {	  
	  doMoreRecursiveRefine = 
	    doMoreRecursiveRefine || (elInfo->getElement()->getMark() > 1);
	  elInfo = refineFunction(elInfo);
	}	
	elInfo = stack->traverseNext(elInfo);
      }
    }


    if (newCoords)
      setNewCoords(macroElIndex); // call of sub-class method

    delete stack;

    nElements -= mesh->getNumberOfLeaves();
    if (nElements != 0)
      aMesh->incChangeIndex();    
  }

} // end namespace AMDiS
