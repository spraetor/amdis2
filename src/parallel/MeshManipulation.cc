/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/


#include "parallel/MeshManipulation.h"
#include "DOFVector.h"
#include "Mesh.h"
#include "MeshStructure.h"
#include "BasisFunction.h"
#include "Traverse.h"
#include "Debug.h"
#include "FiniteElemSpace.h"

using namespace std;

namespace AMDiS { namespace Parallel {

  MeshManipulation::MeshManipulation(Mesh *m)
    : mesh(m)
  {
    switch (mesh->getDim()) {
    case 2:
      refineManager = new RefinementManager2d();
      break;
    case 3:
      refineManager = new RefinementManager3d();
      break;
    default:
      ERROR_EXIT("invalid dim!\n");
    }
  }


  MeshManipulation::~MeshManipulation()
  {
    delete refineManager;
  }

  void MeshManipulation::deleteDoubleDofs(std::vector<const FiniteElemSpace*>& feSpaces,
			  std::vector<MacroElement*>& newMacroEl,
			  ElementObjectDatabase &elObjDb)
  {
    std::set<MacroElement*> newMacroElSet(newMacroEl.begin(), newMacroEl.end());
    deleteDoubleDofs(feSpaces, newMacroElSet, elObjDb);
  }


  void MeshManipulation::deleteDoubleDofs(vector<const FiniteElemSpace*>& feSpaces,
					  std::set<MacroElement*>& newMacroEl, 
					  ElementObjectDatabase &elObjDb)
  {
    FUNCNAME("MeshManipulation::deleteDoubleDofs()");

    // Search for the FE space with the highest degree of polynomials. Using this
    // FE space ensures that deleting DOFs defined on it, also DOFs of lower 
    // order FE spaces will be deleted correctly.
    const FiniteElemSpace *feSpace = FiniteElemSpace::getHighest(feSpaces);

    // Create data structure that maps macro element indices to the 
    // corresponding pointers.
    map<int, MacroElement*> macroIndexMap;
    for (std::set<MacroElement*>::iterator it = newMacroEl.begin();
	 it != newMacroEl.end(); ++it)
      macroIndexMap[(*it)->getIndex()] = *it;
   
    // === Traverse mesh and put all "old" macro element to macrosProcessed  ===
    // === that stores all macro elements which are really connected to the  ===
    // === overall mesh structure.                                           ===

    std::set<int> macrosProcessed;   
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL);
    while (elInfo) {
      if (newMacroEl.count(elInfo->getMacroElement()) == 0) {
	int index = elInfo->getMacroElement()->getIndex();

	macrosProcessed.insert(index);
	macroIndexMap[index] = elInfo->getMacroElement();
      }

      elInfo = stack.traverseNext(elInfo);
    }

#if (DEBUG != 0)
    DOFVector<WorldVector<double> > coords(feSpace, "dofCorrds");
    feSpace->getMesh()->getDofIndexCoords(coords);
#endif

    // === Traverse all new elements and connect them to the overall mesh     ===
    // === structure by deleting all double DOFs on macro element's vertices, ===
    // === edges and faces.                                                   ===

    map<const DegreeOfFreedom*, const DegreeOfFreedom*> mapDelDofs;
    map<const DegreeOfFreedom*, GeoIndex> dofPosIndex;

    for (std::set<MacroElement*>::iterator it = newMacroEl.begin(); 
	 it != newMacroEl.end(); ++it) {

      // === Traverse all vertices of the new element. ===

      for (int i = 0; i < mesh->getGeo(VERTEX); i++) {
	vector<ElementObjectData> &vertexEl = 
	  elObjDb.getElementsVertex((*it)->getIndex(), i);

	for (vector<ElementObjectData>::iterator elIt = vertexEl.begin();
	     elIt != vertexEl.end(); ++elIt) {
	  if (elIt->elIndex == (*it)->getIndex())
	    continue;

	  if (macrosProcessed.count(elIt->elIndex)) {
	    TEST_EXIT_DBG(macroIndexMap.count(elIt->elIndex))
	      ("Should not happen!\n");

	    Element *el0 = (*it)->getElement();
	    Element *el1 = macroIndexMap[elIt->elIndex]->getElement();

	    const DegreeOfFreedom *dof0 = el0->getDof(i);
	    const DegreeOfFreedom *dof1 = el1->getDof(elIt->ithObject);

	    mapDelDofs[dof0] = dof1;
	    dofPosIndex[dof0] = VERTEX;

	    break;
	  } 
	}
      }

      for (int i = 0; i < mesh->getGeo(EDGE); i++) {
	ElementObjectData elObj((*it)->getIndex(), i);

	vector<ElementObjectData> &edgeEl = 
	  elObjDb.getElementsEdge((*it)->getIndex(), i);

      	for (vector<ElementObjectData>::iterator elIt = edgeEl.begin();
	     elIt != edgeEl.end(); ++elIt) {
	  if (elIt->elIndex == (*it)->getIndex())
	    continue;

	  if (macrosProcessed.count(elIt->elIndex)) {
	    TEST_EXIT_DBG(macroIndexMap.count(elIt->elIndex))
	      ("Should not happen!\n");
	    
	    Element *el0 = (*it)->getElement();	    
	    Element *el1 = macroIndexMap[elIt->elIndex]->getElement();

	    TEST_EXIT_DBG(el0->getMesh() == el1->getMesh())("Mesh is different.\n");

	    bool reverseMode = elObjDb.getEdgeReverseMode(elObj, *elIt);

	    BoundaryObject b0(el0, 0, EDGE, i, reverseMode);
	    BoundaryObject b1(el1, 0, EDGE, elIt->ithObject, false);

	    DofContainer dofs0, dofs1;
	    vector<GeoIndex> dofGeoIndex0, dofGeoIndex1;
	    el0->getAllDofs(feSpace, b0, dofs0, true, &dofGeoIndex0);
	    el1->getAllDofs(feSpace, b1, dofs1, true, &dofGeoIndex1);
	    
	    
#if (DEBUG != 0)
	    if (feSpaces.size())
	      debug::testDofsByCoords(coords, dofs0, dofs1);
	    else 
	      TEST_EXIT_DBG(dofs0.size() == dofs1.size())
		("Should not happen!\n");
#endif
 	    for (unsigned int i = 0; i < dofs0.size(); i++) {
 	      mapDelDofs[dofs0[i]] = dofs1[i];
	      dofPosIndex[dofs0[i]] = dofGeoIndex0[i];

	      TEST_EXIT_DBG(dofGeoIndex0[i] == dofGeoIndex1[i])
		("Should not happen: %d %d\n", dofGeoIndex0[i], dofGeoIndex1[i]);
	    }

	    break;
	  }
	}
      }


      for (int i = 0; i < mesh->getGeo(FACE); i++) {
	ElementObjectData elObj((*it)->getIndex(), i);

	vector<ElementObjectData> &faceEl = 
	  elObjDb.getElementsFace((*it)->getIndex(), i);

	for (vector<ElementObjectData>::iterator elIt = faceEl.begin();
	     elIt != faceEl.end(); ++elIt) {
	  if (elIt->elIndex == (*it)->getIndex())
	    continue;

	  if (macrosProcessed.count(elIt->elIndex)) {
	    TEST_EXIT_DBG(macroIndexMap.count(elIt->elIndex))
	      ("Should not happen!\n");

	    Element *el0 = (*it)->getElement();	    
	    Element *el1 = macroIndexMap[elIt->elIndex]->getElement();

	    bool reverseMode = elObjDb.getFaceReverseMode(elObj, *elIt);

	    BoundaryObject b0(el0, 0, FACE, i, reverseMode);
	    BoundaryObject b1(el1, 0, FACE, elIt->ithObject, false);

	    DofContainer dofs0, dofs1;
	    vector<GeoIndex> dofGeoIndex0, dofGeoIndex1;
	    el0->getAllDofs(feSpace, b0, dofs0, true, &dofGeoIndex0);
	    el1->getAllDofs(feSpace, b1, dofs1, true, &dofGeoIndex1);

#if (DEBUG != 0)
	    if (feSpaces.size())
	      debug::testDofsByCoords(coords, dofs0, dofs1);
	    else 
	      TEST_EXIT_DBG(dofs0.size() == dofs1.size())
		("Should not happen!\n");
#endif

 	    for (unsigned int i = 0; i < dofs0.size(); i++) {
 	      mapDelDofs[dofs0[i]] = dofs1[i];
	      dofPosIndex[dofs0[i]] = dofGeoIndex1[i];

	      TEST_EXIT_DBG(dofGeoIndex0[i] == dofGeoIndex1[i])
		("Should not happen: %d %d\n", 
		 dofGeoIndex0[i], dofGeoIndex1[i]);
	    }

	    break;
	  }
	}
      }
      macrosProcessed.insert((*it)->getIndex());
    }


    // === Remove all DOF replacments of the form A -> B, B -> C by A -> C. ===

    bool changed = false;
    do {
      changed = false;
      for (map<const DegreeOfFreedom*, const DegreeOfFreedom*>::iterator it = mapDelDofs.begin();
	   it != mapDelDofs.end(); ++it) {
	map<const DegreeOfFreedom*, const DegreeOfFreedom*>::iterator findIt = mapDelDofs.find(it->second);
	if (findIt != mapDelDofs.end()) {
	  TEST_EXIT_DBG(it->first != findIt->second)
	    ("Cycle %d -> %d and %d -> %d! Should not happen!\n",
	     *(it->first), *(it->second), *(findIt->first), *(findIt->second));
	  it->second = findIt->second;   
	  changed = true;
	}
      } 
    } while (changed);


    // === Set new DOF pointers in all elements of the mesh. ===

    elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
    while (elInfo) {
      for (int i = 0; i < mesh->getNumberOfNodes(); i++)
	if (mapDelDofs.count(elInfo->getElement()->getDof(i)) == 1)
	  elInfo->getElement()->setDof(i, const_cast<DegreeOfFreedom*>(mapDelDofs[elInfo->getElement()->getDof(i)]));

      elInfo = stack.traverseNext(elInfo);
    }


    // === And delete all double DOFs. ===

    for (map<const DegreeOfFreedom*, const DegreeOfFreedom*>::iterator it = mapDelDofs.begin();
	 it != mapDelDofs.end(); ++it)
      mesh->freeDof(const_cast<DegreeOfFreedom*>(it->first), 
		    dofPosIndex[it->first]);
  }


  bool MeshManipulation::fitElementToMeshCode(MeshStructure &code,
					      BoundaryObject &boundEl)
  {
    FUNCNAME("MeshManipulation::fitElementToMeshCode()");

    TEST_EXIT_DBG(boundEl.el)("No element given!\n");

    // If the code is empty, the element does not matter and the function can
    // return without chaning the mesh.
    if (code.empty())
      return false;

    // s0 and s1 are the number of the edge/face in both child of the element,
    // which contain the edge/face the function has to traverse through. If the
    // edge/face is not contained in one of the children, s0 or s1 is -1.
    int s0 = boundEl.el->getSubObjOfChild(0, boundEl.subObj, 
					  boundEl.ithObj, boundEl.elType);
    int s1 = boundEl.el->getSubObjOfChild(1, boundEl.subObj, 
					  boundEl.ithObj, boundEl.elType);

    TEST_EXIT_DBG(s0 != -1 || s1 != -1)("This should not happen!\n");

    bool meshChanged = false;
    Flag traverseFlag = 
      Mesh::CALL_EVERY_EL_PREORDER | Mesh::FILL_NEIGH | Mesh::FILL_BOUND;

    // Test for reverse mode, in which the left and right children of elements
    // are flipped.
    if (boundEl.reverseMode)
      traverseFlag |= Mesh::CALL_REVERSE_MODE;    


    // === If the edge/face is contained in both children. ===

    Mesh *mesh = boundEl.el->getMesh();

    if (s0 != -1 && s1 != -1) {
      // Create traverse stack and traverse within the mesh until the element,
      // which should be fitted to the mesh structure code, is reached.
      TraverseStack stack;
#if (DEBUG != 0)
      ElInfo *elInfo = 
	stack.traverseFirstOneMacro(mesh, boundEl.elIndex, -1, traverseFlag);
	
      TEST_EXIT(elInfo->getElement() == boundEl.el)
	("This should not happen!\n");
#else
      // stack must be initialized before passed to fitElementToMeshCode()
      stack.traverseFirstOneMacro(mesh, boundEl.elIndex, -1, traverseFlag);
#endif



      pcode = &code;
      pstack = &stack;
      rMode = boundEl.reverseMode;
      return fitElementToMeshCode(boundEl.subObj, boundEl.ithObj);
    }


    // === The edge/face is contained in only on of the both children. ===

    if (boundEl.el->isLeaf()) {

      // Create traverse stack and traverse the mesh to the element.
      TraverseStack stack;
      ElInfo *elInfo = 
	stack.traverseFirstOneMacro(mesh, boundEl.elIndex, -1, traverseFlag);

      TEST_EXIT_DBG(elInfo)("This should not happen!\n");

      // Code is not leaf, therefore refine the element.
      boundEl.el->setMark(1);
      refineManager->setMesh(boundEl.el->getMesh());
      refineManager->setStack(&stack);
      refineManager->refineFunction(elInfo);
      meshChanged = true;
      // If element is leaf and code contains only one leaf element, we are finished.
      if (code.getNumElements() == 1 && code.isLeafElement())
	return true;     


    }

    Element *child0 = boundEl.el->getFirstChild();
    Element *child1 = boundEl.el->getSecondChild();
    if (boundEl.reverseMode) {
      using std::swap;
      swap(s0, s1);
      swap(child0, child1);    
    }

    // === We know that the edge/face is contained in only one of the children. ===
    // === Therefore, traverse the mesh to this children and fit this element   ===
    // === To the mesh structure code.                                          ===

    TraverseStack stack;
    ElInfo *elInfo = 
      stack.traverseFirstOneMacro(mesh, boundEl.elIndex, -1, traverseFlag);

    if (s0 != -1) {
      while (elInfo && elInfo->getElement() != child0)
	elInfo = stack.traverseNext(elInfo);     

      pcode = &code;
      pstack = &stack;
      rMode = boundEl.reverseMode;
      meshChanged |= fitElementToMeshCode(boundEl.subObj, s0);
    } else {
      while (elInfo && elInfo->getElement() != child1) 
	elInfo = stack.traverseNext(elInfo);      

      pcode = &code;
      pstack = &stack;
      rMode = boundEl.reverseMode;
      meshChanged |= fitElementToMeshCode(boundEl.subObj, s1);
    }


    return meshChanged;
  }


  bool MeshManipulation::fitElementToMeshCode(GeoIndex subObj, int ithObj)
  {
    FUNCNAME("MeshManipulation::fitElementToMeshCode()");

    // === Test if there are more elements in stack to check with the code. ===

    ElInfo *elInfo = pstack->getElInfo();
    if (!elInfo)
      return false;


    // === Test if code contains a leaf element. If this is the case, the ===
    // === current element is finished. Traverse the mesh to the next     ===
    // === coarser element.                                               ===

    if (pcode->isLeafElement()) {
      int level = elInfo->getLevel();

      do {
	elInfo = pstack->traverseNext(elInfo);
      } while (elInfo && elInfo->getLevel() > level);

      return false;
    }


    bool meshChanged = false;
    Element *el = elInfo->getElement();
    int s0 = el->getSubObjOfChild(0, subObj, ithObj, elInfo->getType());
    int s1 = el->getSubObjOfChild(1, subObj, ithObj, elInfo->getType());

    // === If element is leaf (and the code is not), refine the element. ===

    bool elementRefined = true;

    if (el->isLeaf()) {
      TEST_EXIT_DBG(elInfo->getLevel() < 255)("This should not happen!\n");

      // In some situations refinement of an element can be omitted, altough 
      // it is defined in the mesh structure code. One case is when the
      // following holds:
      //   - the code is used to adapt a face of a tetrahedron
      //   - only one of the children of the current element contains this face
      //   - the face to be adapted of the current element is either the
      //     local face 0 or 1
      //   - the next structure code is 0, i.e., we would be finished after
      //     the refinement of this element
      // In this scenario, the element must be refined due to the structure
      // code, but the refinement does not introduce new DOFs on the face,
      // that should be adapted. Thus, we can ommit the refinement.
      if (subObj == FACE) {
  	if ((s0 != -1 && s1 == -1) || (s0 == -1 && s1 != -1)) {
  	  if (ithObj <= 1 && pcode->lookAhead() == 0) {
  	    elementRefined = false;
 	    pcode->nextElement();
 	    pstack->traverseNext(elInfo);
  	  }
  	}
      }

      if (elementRefined) {
	el->setMark(1);
	refineManager->setMesh(el->getMesh());
	refineManager->setStack(pstack);
	refineManager->refineFunction(elInfo);
	meshChanged = true;
      } 
    }


    if (elementRefined) {
      // === Continue fitting the mesh structure code to the children of the ===
      // === current element.                                                ===
      
      Element *child0 = el->getFirstChild();
      Element *child1 = el->getSecondChild();
      if (rMode) {
	using std::swap;
	swap(s0, s1);
	swap(child0, child1);
      }

    
      // === Traverse left child. ===
      
      if (s0 != -1) {
	// The edge/face is contained in the left child, therefore fit this
	// child to the mesh structure code.
	
	pstack->traverseNext(elInfo);
	pcode->nextElement();
	meshChanged |= fitElementToMeshCode(subObj, s0);
	elInfo = pstack->getElInfo();
      } else {
	// The edge/face is not contained in the left child. Hence we need
	// to traverse through all subelements of the left child until we get
	// the second child of the current element.
	
	do {
	  elInfo = pstack->traverseNext(elInfo);
	} while (elInfo && elInfo->getElement() != child1); 
	
	TEST_EXIT_DBG(elInfo != NULL)("This should not happen!\n");
      }  
      
      TEST_EXIT_DBG(elInfo->getElement() == child1)
	("Got wrong child with idx = %d! Searched for child idx = %d\n",
	 elInfo->getElement()->getIndex(), child1->getIndex());
      
      
      // === Traverse right child. ===
      
      if (s1 != -1) {
	// The edge/face is contained in the right child, therefore fit this
	// child to the mesh structure code.
	
	pcode->nextElement();
	meshChanged |= fitElementToMeshCode(subObj, s1);
      } else {
	// The edge/face is not contained in the right child. Hence we need
	// to traverse through all subelements of the right child until we have
	// finished traversing the current element with all its subelements.
	
	int level = elInfo->getLevel();
	
	do {
	  elInfo = pstack->traverseNext(elInfo);
	} while (elInfo && elInfo->getLevel() > level);
      }
      
    }
      
    return meshChanged;
  }
  
} }
