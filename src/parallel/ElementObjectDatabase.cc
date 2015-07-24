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


#include "VertexVector.h"
#include "parallel/ElementObjectDatabase.h"
#include "parallel/MeshLevelData.h"
#include <fstream>

#include "Serializer.h"

using namespace std;

namespace AMDiS { namespace Parallel {
  
  void ElementObjectDatabase::init(vector<Mesh*>& meshes_, 
				   map<Mesh*, vector<const FiniteElemSpace*> >& meshToFeSpaces)
  {
    meshes = meshes_;
    for(size_t i = 0; i < meshes.size(); i++) 
      feSpaces.push_back(meshToFeSpaces[meshes[i]][0]);
    
    macroMesh = meshes[0];
    feSpace = meshToFeSpaces[macroMesh][0];
  }

  void ElementObjectData::serialize(std::ostream &out) const
  {
    SerUtil::serialize(out, elIndex);
    SerUtil::serialize(out, ithObject);
  }

  /// Read this element object from disk.
  void ElementObjectData::deserialize(std::istream &in)
  {
    SerUtil::deserialize(in, elIndex);
    SerUtil::deserialize(in, ithObject);
  }


  void ElementObjectDatabase::create(map<int, int>& rankMap, MeshLevelData& ld)
  {
    FUNCNAME_DBG("ElementObjectDatabase::create()");

    macroElementRankMap = &rankMap;
    levelData = &ld;

    // === Reset temporary data ===
    
    tmpVertexElements.clear();
    tmpEdgeElements.clear();
    tmpFaceElements.clear();

    // === Fills macro element data structures. ===
    
    TraverseStack stack;
    ElInfo* elInfo = NULL;
    
    for(size_t i = 0; i < meshes.size(); i++) {
      
      elInfo = stack.traverseFirst(meshes[i], -1,
				   Mesh::CALL_LEAF_EL | 
				   Mesh::FILL_NEIGH | 
				   Mesh::FILL_BOUND);
      while (elInfo) {
	TEST_EXIT_DBG(elInfo->getLevel() == 0)("Should not happen!\n");
	
	Element *el = elInfo->getElement();
	
	// Macro data info is stored once.
	if (i == 0) {
	  macroElIndexTypeMap.insert(make_pair(el->getIndex(), elInfo->getType()));
	
	  // Add all sub object of the element to the variable elObjDb.
	  addElement(el);
	  
	  addElementPeriodicBoundary(elInfo);
	}

        // Element pointer has to be stored for each mesh one copy.
	macroElIndexMap[el->getIndex()][meshes[i]] = el;
	
	elInfo = stack.traverseNext(elInfo);
      }
    }
    
    // === Move temporary data to original one ===
    
    for (map<DegreeOfFreedom, vector<ElementObjectData> >::iterator it = tmpVertexElements.begin();
	 it != tmpVertexElements.end(); ++it)
      vertexElements[it->first] = it->second;
    tmpVertexElements.clear();

    for (map<DofEdge, vector<ElementObjectData> >::iterator it = tmpEdgeElements.begin();
	 it != tmpEdgeElements.end(); ++it)
      edgeElements[it->first] = it->second;
    tmpEdgeElements.clear();

    for (map<DofFace, vector<ElementObjectData> >::iterator it = tmpFaceElements.begin();
	 it != tmpFaceElements.end(); ++it)
      faceElements[it->first] = it->second;
    tmpFaceElements.clear();

    // Handle periodic boudaries
    createPeriodicData();

    // Create data about the reverse modes of neighbouring elements.
    createReverseModeData();
  }
  
  void ElementObjectDatabase::createMacroElementInfo(std::map<Mesh*, std::vector<MacroElement*> >& mel)
  {
    macroElIndexMap.clear();
    macroElIndexTypeMap.clear();
    
    std::map<Mesh*, std::vector<MacroElement*> >::iterator iter = mel.begin();
    while(iter != mel.end()) {
      Mesh* mesh = iter->first;
      for (vector<MacroElement*>::iterator it = iter->second.begin();
	   it != iter->second.end(); it++) {
	macroElIndexMap[(*it)->getIndex()][mesh] = (*it)->getElement();
        if(iter == mel.begin())
          macroElIndexTypeMap.insert(make_pair((*it)->getIndex(), (*it)->getElType()));
      }
      iter++;
    }
  }
      
  void ElementObjectDatabase::addElement(Element *el)
  {
    // === First, add all element objects to the database. ===

    for (int i = 0; i < el->getGeo(VERTEX); i++)
      addVertex(el, i);
    
    for (int i = 0; i < el->getGeo(EDGE); i++)
      addEdge(el, i);      
    
    for (int i = 0; i < el->getGeo(FACE); i++)
      addFace(el, i);
  }
  
  void ElementObjectDatabase::addElementPeriodicBoundary(ElInfo* elInfo)
  {
    FUNCNAME("ElementObjectDatabase::addElementPeriodicBoundary()");

    TEST_EXIT_DBG(macroMesh)("Mesh not set!\n");

    Element *el = elInfo->getElement();
    
    switch (macroMesh->getDim()) {
    case 2:
      for (int i = 0; i < el->getGeo(EDGE); i++) {
	if (macroMesh->isPeriodicAssociation(elInfo->getBoundary(EDGE, i))) {
	  // The current element's i-th edge is periodic.
	  Element *neigh = elInfo->getNeighbour(i);

	  TEST_EXIT_DBG(neigh)("Should not happen!\n");
	  
	  DofEdge edge0 = el->getEdge(i);
	  DofEdge edge1 = neigh->getEdge(elInfo->getOppVertex(i));
	  BoundaryType boundaryType = elInfo->getBoundary(EDGE, i);

	  // Add the periodic edge.
	  periodicEdges[make_pair(edge0, edge1)] = boundaryType;
	  periodicEdgeAssoc[edge0].insert(edge1);
	  
	  // Add both vertices of the edge to be periodic.
	  DegreeOfFreedom mappedDof0 = 
	    macroMesh->getPeriodicAssociations(boundaryType)[edge0.first];
	  DegreeOfFreedom mappedDof1 = 
	    macroMesh->getPeriodicAssociations(boundaryType)[edge0.second];
	  periodicVertices[make_pair(edge0.first, mappedDof0)] = boundaryType;
	  periodicVertices[make_pair(edge0.second, mappedDof1)] = boundaryType;

	  periodicDofAssoc[edge0.first].insert(boundaryType);
	  periodicDofAssoc[edge0.second].insert(boundaryType); 
	}
      }
      break;
    case 3:
      for (int i = 0; i < el->getGeo(FACE); i++) {
	if (macroMesh->isPeriodicAssociation(elInfo->getBoundary(FACE, i))) {
	  // The current element's i-th face is periodic.
	  Element *neigh = elInfo->getNeighbour(i);

	  TEST_EXIT_DBG(neigh)("Should not happen!\n");
	  
	  DofFace face0 = el->getFace(i);
	  DofFace face1 = neigh->getFace(elInfo->getOppVertex(i));
	  BoundaryType boundaryType = elInfo->getBoundary(FACE, i);

	  // Add the periodic face.
	  periodicFaces[make_pair(face0, face1)] = elInfo->getBoundary(i);
	  
	  /// Add all three vertices of the face to be periodic.
	  DegreeOfFreedom mappedDof0 = 
	    macroMesh->getPeriodicAssociations(boundaryType)[std::get<0>(face0)];
	  DegreeOfFreedom mappedDof1 = 
	    macroMesh->getPeriodicAssociations(boundaryType)[std::get<1>(face0)];
	  DegreeOfFreedom mappedDof2 = 
	    macroMesh->getPeriodicAssociations(boundaryType)[std::get<2>(face0)];

	  periodicVertices[make_pair(std::get<0>(face0), mappedDof0)] = boundaryType;
	  periodicVertices[make_pair(std::get<1>(face0), mappedDof1)] = boundaryType;
	  periodicVertices[make_pair(std::get<2>(face0), mappedDof2)] = boundaryType;
	  
	  periodicDofAssoc[std::get<0>(face0)].insert(boundaryType);
	  periodicDofAssoc[std::get<1>(face0)].insert(boundaryType);
	  periodicDofAssoc[std::get<2>(face0)].insert(boundaryType);
	  
	  TEST_EXIT_DBG(std::get<0>(face0) == 
			macroMesh->getPeriodicAssociations(boundaryType)[std::get<0>(face1)] &&
			std::get<1>(face0) == 
			macroMesh->getPeriodicAssociations(boundaryType)[std::get<1>(face1)] &&
			std::get<2>(face0) == 
			macroMesh->getPeriodicAssociations(boundaryType)[std::get<2>(face1)])
	    ("Should not happen!\n");
	  
	  // Create all three edges of the element and add them to be periodic.
	  DofEdge elEdge0 = make_pair(std::get<0>(face0), std::get<1>(face0));
	  DofEdge elEdge1 = make_pair(std::get<0>(face0), std::get<2>(face0));
	  DofEdge elEdge2 = make_pair(std::get<1>(face0), std::get<2>(face0));
	  DofEdge neighEdge0 = make_pair(std::get<0>(face1), std::get<1>(face1));
	  DofEdge neighEdge1 = make_pair(std::get<0>(face1), std::get<2>(face1));
	  DofEdge neighEdge2 = make_pair(std::get<1>(face1), std::get<2>(face1));
	  	  
	  periodicEdges[make_pair(elEdge0, neighEdge0)] = boundaryType;
	  periodicEdges[make_pair(elEdge1, neighEdge1)] = boundaryType;
	  periodicEdges[make_pair(elEdge2, neighEdge2)] = boundaryType;
	  
	  periodicEdgeAssoc[elEdge0].insert(neighEdge0);
	  periodicEdgeAssoc[elEdge1].insert(neighEdge1);
	  periodicEdgeAssoc[elEdge2].insert(neighEdge2);
	}
      }
      break;
    default:
      ERROR_EXIT("Should not happen!\n");
    }  
    
  }


  void ElementObjectDatabase::addVertex(Element *el, int ith)
  {
    DegreeOfFreedom vertex = el->getDof(ith, 0);
    int elIndex = el->getIndex();
    ElementObjectData elObj(elIndex, ith);
   
    tmpVertexElements[vertex].push_back(elObj);
    vertexLocalMap[elObj] = vertex;
  }


  void ElementObjectDatabase::addEdge(Element *el, int ith)
  {
    DofEdge edge = el->getEdge(ith);
    int elIndex = el->getIndex();
    ElementObjectData elObj(elIndex, ith);
    
    tmpEdgeElements[edge].push_back(elObj);
    edgeLocalMap[elObj] = edge;
  }

  
  void ElementObjectDatabase::addFace(Element *el, int ith)
  {
    DofFace face = el->getFace(ith);
    int elIndex = el->getIndex();
    ElementObjectData elObj(elIndex, ith);
    
    tmpFaceElements[face].push_back(elObj);
    faceLocalMap[elObj] = face;
  }


  void ElementObjectDatabase::createPeriodicData()
  {
    FUNCNAME_DBG("ElementObjectDatabase::createPeriodicData()");

    TEST_EXIT_DBG(macroMesh)("Mesh not set!\n");

    // === Return, if there are no periodic vertices, i.e., there are no  ===
    // === periodic boundaries in the mesh.                               ===
    
    if (periodicVertices.size() == 0)
      return;
 
    // === Calculate smallest periodic boundary ID in mesh. ===

    smallestPeriodicBcType = 0;
    for (map<BoundaryType, VertexVector*>::iterator it = macroMesh->getPeriodicAssociations().begin();
	 it != macroMesh->getPeriodicAssociations().end(); ++it)
      smallestPeriodicBcType = std::min(smallestPeriodicBcType, it->first);

    // === Get all vertex DOFs that have multiple periodic associations. ===
    
    // We group all vertices together, that have either two or three periodic
    // associations. For rectangular domains in 2D, the four corner vertices have
    // all two periodic associations. For box domains in 3D, the eight corner
    // vertices have all three periodic associations.

    vector<DegreeOfFreedom> multPeriodicDof2, multPeriodicDof3;
    for (map<DegreeOfFreedom, std::set<BoundaryType> >::iterator it = periodicDofAssoc.begin();
	 it != periodicDofAssoc.end(); ++it) {
      TEST_EXIT_DBG((macroMesh->getDim() == 2 && it->second.size() <= 2) ||
		    (macroMesh->getDim() == 3 && it->second.size() <= 3))
	("Should not happen!\n");
      
      if (it->second.size() == 2)
	multPeriodicDof2.push_back(it->first);
      if (it->second.size() == 3)
	multPeriodicDof3.push_back(it->first);
    }

    if (macroMesh->getDim() == 2) {
      TEST_EXIT_DBG(multPeriodicDof2.size() == 0 || 
		    multPeriodicDof2.size() == 4)
	("Should not happen (%d)!\n", multPeriodicDof2.size());
      TEST_EXIT_DBG(multPeriodicDof3.size() == 0)("Should not happen!\n");
    }

    if (macroMesh->getDim() == 3) {
      TEST_EXIT_DBG(multPeriodicDof3.size() == 0 || 
		    multPeriodicDof3.size() == 8)
	("Should not happen (%d)!\n", multPeriodicDof3.size());
    }

    if (multPeriodicDof2.size() > 0) {
      for (unsigned int i = 0; i < multPeriodicDof2.size(); i++) {
	DegreeOfFreedom dof0 = multPeriodicDof2[i];
	if (dof0 == -1)
	  continue;
	
	DegreeOfFreedom dof1 = -1;
	DegreeOfFreedom dof2 = -1;
	BoundaryType type0 = *(periodicDofAssoc[dof0].begin());
	BoundaryType type1 = *(++(periodicDofAssoc[dof0].begin()));
	
	for (PerBoundMap<DegreeOfFreedom>::iterator it = periodicVertices.begin();
	     it != periodicVertices.end(); ++it) {
	  if (it->first.first == dof0 && it->second == type0)
	    dof1 = it->first.second;
	  if (it->first.first == dof0 && it->second == type1)
	    dof2 = it->first.second;
	  
	  if (dof1 != -1 && dof2 != -1)
	    break;
	}
	
	TEST_EXIT_DBG(dof1 != -1 && dof2 != -1)("Should not happen!\n");
	
	DegreeOfFreedom dof3 = -1;
	for (PerBoundMap<DegreeOfFreedom>::iterator it = periodicVertices.begin();
	     it != periodicVertices.end(); ++it) {
	  if (it->first.first == dof1 && it->second == type1) {
	    dof3 = it->first.second;
	    
	    TEST_EXIT_DBG(periodicVertices[make_pair(dof2, dof3)] == type0)
	      ("Should not happen!\n");
	    
	    break;
	  }
	}
	  
	TEST_EXIT_DBG(dof3 != -1)("Should not happen!\n");
	TEST_EXIT_DBG(periodicVertices.count(make_pair(dof0, dof3)) == 0)
	  ("Should not happen!\n");
	TEST_EXIT_DBG(periodicVertices.count(make_pair(dof3, dof0)) == 0)
	  ("Should not happen!\n");

 	periodicVertices[make_pair(dof0, dof3)] = 
	  provideConnectedPeriodicBoundary(type0, type1);
 	periodicVertices[make_pair(dof3, dof0)] = 
 	  provideConnectedPeriodicBoundary(type0, type1);

	for (unsigned int j = i + 1; j < multPeriodicDof2.size(); j++)
	  if (multPeriodicDof2[j] == dof3)
	    multPeriodicDof2[j] = -1;
      }
    }

    if (multPeriodicDof3.size() > 0) {
      int nMultPeriodicDofs = multPeriodicDof3.size();

      for (int i = 0; i < nMultPeriodicDofs; i++) {
	for (int j = i + 1; j < nMultPeriodicDofs; j++) {
	  pair<DegreeOfFreedom, DegreeOfFreedom> perDofs0 = 
	    make_pair(multPeriodicDof3[i], multPeriodicDof3[j]);
	  pair<DegreeOfFreedom, DegreeOfFreedom> perDofs1 = 
	    make_pair(multPeriodicDof3[j], multPeriodicDof3[i]);
	  
	  if (periodicVertices.count(perDofs0) == 0) {
	    BoundaryType b = getNewBoundaryType();
	    periodicVertices[perDofs0] = b;
	    periodicVertices[perDofs1] = b;
	  }
	}
      }
    }
  
   
    // === Get all edges that have multiple periodic associations (3D only!). ===
    
    for (map<DofEdge, std::set<DofEdge> >::iterator it = periodicEdgeAssoc.begin();
	 it != periodicEdgeAssoc.end(); ++it) {
      if (it->second.size() > 1) {
	TEST_EXIT_DBG(macroMesh->getDim() == 3)("Should not happen!\n");
	TEST_EXIT_DBG(it->second.size() == 2)("Should not happen!\n");
	
	DofEdge edge0 = it->first;
	DofEdge edge1 = *(it->second.begin());
	DofEdge edge2 = *(++(it->second.begin()));
	
	pair<DofEdge, DofEdge> perEdge0 = make_pair(edge1, edge2);
	pair<DofEdge, DofEdge> perEdge1 = make_pair(edge2, edge1);

	TEST_EXIT_DBG(periodicEdges.count(make_pair(edge0, edge1)) == 1)
	  ("Should not happen!\n");
	TEST_EXIT_DBG(periodicEdges.count(make_pair(edge1, edge0)) == 1)
	  ("Should not happen!\n");

	BoundaryType type0 = periodicEdges[make_pair(edge0, edge1)];
	BoundaryType type1 = periodicEdges[make_pair(edge0, edge2)];
	BoundaryType type2 = provideConnectedPeriodicBoundary(type0, type1);

 	periodicEdges[perEdge0] = type2;
 	periodicEdges[perEdge1] = type2;
      }
    }


    // === In debug mode we make some tests, if the periodic structures are set ===
    // === in a symmetric way, i.e., if A -> B for a specific boundary type,    ===
    // === there must be a mapping B -> A with the same boundary type.          ===
    
#if (DEBUG != 0)
    for (PerBoundMap<DegreeOfFreedom>::iterator it = periodicVertices.begin();
	 it != periodicVertices.end(); ++it) {
      pair<DegreeOfFreedom, DegreeOfFreedom> testVertex = 
	make_pair(it->first.second, it->first.first);
      
      TEST_EXIT_DBG(periodicVertices.count(testVertex) == 1)
	("Should not happen!\n");
      TEST_EXIT_DBG(periodicVertices[testVertex] == it->second)
	("Should not happen!\n");
    }

    for (PerBoundMap<DofEdge>::iterator it = periodicEdges.begin();
	 it != periodicEdges.end(); ++it) {
      pair<DofEdge, DofEdge> testEdge = 
	make_pair(it->first.second, it->first.first);
      
      TEST_EXIT_DBG(periodicEdges.count(testEdge) == 1)("Should not happen!\n");
      TEST_EXIT_DBG(periodicEdges[testEdge] == it->second)("Should not happen!\n");
    }
    
    for (PerBoundMap<DofFace>::iterator it = periodicFaces.begin();
	 it != periodicFaces.end(); ++it) {
      pair<DofFace, DofFace> testFace = 
	make_pair(it->first.second, it->first.first);
      
      TEST_EXIT_DBG(periodicFaces.count(testFace) == 1)("Should not happen!\n");
      TEST_EXIT_DBG(periodicFaces[testFace] == it->second)("Should not happen!\n");
    }
#endif       
  }

  BoundaryType ElementObjectDatabase::getNewBoundaryType()
  {
    FUNCNAME_DBG("ElementObjectDatabase::getNewBoundaryType()");

    BoundaryType newPeriodicBoundaryType = 0;
    for (map<BoundaryType, VertexVector*>::iterator it = macroMesh->getPeriodicAssociations().begin();
	 it != macroMesh->getPeriodicAssociations().end(); ++it)
      newPeriodicBoundaryType = std::min(newPeriodicBoundaryType, it->first);
    
    TEST_EXIT_DBG(newPeriodicBoundaryType < 0)("Should not happen!\n");
    newPeriodicBoundaryType--;
    
    for(size_t i = 0; i < meshes.size(); i++) {
      meshes[i]->getPeriodicAssociations()[newPeriodicBoundaryType] = 
        new VertexVector(feSpaces[i]->getAdmin(), "");
    }
        
    return newPeriodicBoundaryType;
  }
  
  BoundaryType 
  ElementObjectDatabase::provideConnectedPeriodicBoundary(BoundaryType b0, 
							  BoundaryType b1)
  {
    FUNCNAME_DBG("ElementObjectDatabase::provideConnectedPeriodicBoundary()");

    std::pair<BoundaryType, BoundaryType> bConn = 
      (b0 <= b1 ? make_pair(b0, b1) : make_pair(b1, b0));

    if (bConnMap.count(bConn) == 0) {
      BoundaryType newPeriodicBoundaryType = getNewBoundaryType();

      for(size_t i = 0; i < meshes.size(); i++) {
      
        VertexVector &vecB0 = meshes[i]->getPeriodicAssociations(b0);
        VertexVector &vecB1 = meshes[i]->getPeriodicAssociations(b1);
        VertexVector &vecC = meshes[i]->getPeriodicAssociations(newPeriodicBoundaryType);

        DOFIteratorBase it(const_cast<DOFAdmin*>(feSpaces[i]->getAdmin()), USED_DOFS);

	for (it.reset(); !it.end(); ++it) {
	  if (!it.isDofFree()) {
	    TEST_EXIT_DBG(vecB1[vecB0[it.getDOFIndex()]] == vecB0[vecB1[it.getDOFIndex()]])
	      ("Should not happen!\n");

	    vecC[it.getDOFIndex()] = vecB1[vecB0[it.getDOFIndex()]];
	  }
	}
      }
      
      bConnMap[bConn] = newPeriodicBoundaryType;
    }
    
    return bConnMap[bConn];
  }


  void ElementObjectDatabase::updateRankData()
  {
    FUNCNAME("ElementObjectDatabase::updateRankData()");
    
    TEST_EXIT(macroElementRankMap)("Should not happen!\n");

    vertexInRank.clear();
    for (flat_map<DegreeOfFreedom, vector<ElementObjectData> >::iterator it = vertexElements.begin();
	 it != vertexElements.end(); ++it) {
      for (vector<ElementObjectData>::iterator it2 = it->second.begin(); 
	   it2 != it->second.end(); ++it2) {
	int elementInRank = (*macroElementRankMap)[it2->elIndex];	
	if (it2->elIndex > vertexInRank[it->first][elementInRank].elIndex)
	  vertexInRank[it->first][elementInRank] = *it2;
      }
    }
    
    edgeInRank.clear();
    for (flat_map<DofEdge, vector<ElementObjectData> >::iterator it = edgeElements.begin();
	 it != edgeElements.end(); ++it) {
      for (vector<ElementObjectData>::iterator it2 = it->second.begin(); 
	   it2 != it->second.end(); ++it2) {
	int elementInRank = (*macroElementRankMap)[it2->elIndex];
	if (it2->elIndex > edgeInRank[it->first][elementInRank].elIndex)
	  edgeInRank[it->first][elementInRank] = *it2;
      }
    }
    
    faceInRank.clear();
    for (flat_map<DofFace, vector<ElementObjectData> >::iterator it = faceElements.begin();
	 it != faceElements.end(); ++it) {
      for (vector<ElementObjectData>::iterator it2 = it->second.begin(); 
	   it2 != it->second.end(); ++it2) {
	int elementInRank = (*macroElementRankMap)[it2->elIndex];
	if (it2->elIndex > faceInRank[it->first][elementInRank].elIndex)
	  faceInRank[it->first][elementInRank] = *it2;
      }
    }
  }


  void ElementObjectDatabase::clear()
  {
    vertexElements.clear();
    edgeElements.clear();
    faceElements.clear();

    vertexLocalMap.clear();
    edgeLocalMap.clear();
    faceLocalMap.clear();

    vertexInRank.clear();
    edgeInRank.clear();
    faceInRank.clear();

    bConnMap.clear();
    periodicVertices.clear();
    periodicEdges.clear();
    periodicFaces.clear();
    periodicDofAssoc.clear();
    periodicEdgeAssoc.clear();

    edgeReverseMode.clear();
    faceReverseMode.clear();

    macroElementRankMap = NULL;
    macroElIndexMap.clear();
    macroElIndexTypeMap.clear();
  }


  void ElementObjectDatabase::createReverseModeData()
  {
    FUNCNAME_DBG("ElementObjectDatabase::createReverseModeData()");

    // === In 2D, all reverse modes are always true! ===

    if (macroMesh->getDim() == 2)
      return;


    // === First, create reverse modes for all "directly" neighbouring elements. ===

    for (flat_map<DofEdge, vector<ElementObjectData> >::iterator edgeIt = edgeElements.begin();
	 edgeIt != edgeElements.end(); ++edgeIt) {

      vector<ElementObjectData>& els = edgeIt->second;

      for (unsigned int i = 0; i < els.size(); i++) {
	BoundaryObject obj0(macroElIndexMap[els[i].elIndex][macroMesh], 
			    macroElIndexTypeMap[els[i].elIndex], 
			    EDGE, els[i].ithObject);

	for (unsigned int j = i + 1; j < els.size(); j++) {
	  BoundaryObject obj1(macroElIndexMap[els[j].elIndex][macroMesh], 
			      macroElIndexTypeMap[els[j].elIndex], 
			      EDGE, els[j].ithObject);

	  bool reverseMode = 
	    BoundaryObject::computeReverseMode(obj0, obj1, feSpace, INTERIOR);
	  if (reverseMode) {
	    edgeReverseMode.insert(make_pair(els[i], els[j]));
	    edgeReverseMode.insert(make_pair(els[j], els[i]));
	  }
	}
      }
    }

    for (flat_map<DofFace, vector<ElementObjectData> >::iterator faceIt = faceElements.begin();
	 faceIt != faceElements.end(); ++faceIt) {

      vector<ElementObjectData>& els = faceIt->second;

      for (unsigned int i = 0; i < els.size(); i++) {
	BoundaryObject obj0(macroElIndexMap[els[i].elIndex][macroMesh], 
			    macroElIndexTypeMap[els[i].elIndex], 
			    FACE, els[i].ithObject);

	for (unsigned int j = i + 1; j < els.size(); j++) {
	  BoundaryObject obj1(macroElIndexMap[els[j].elIndex][macroMesh], 
			      macroElIndexTypeMap[els[j].elIndex], 
			      FACE, els[j].ithObject);

	  bool reverseMode = 
	    BoundaryObject::computeReverseMode(obj0, obj1, feSpace, INTERIOR);
	  if (reverseMode) {
	    faceReverseMode.insert(make_pair(els[i], els[j]));
	    faceReverseMode.insert(make_pair(els[j], els[i]));
	  }
	}
      }
    }


    // === And create reverse modes for periodic neighbouring elements. ===

    for (PerBoundMap<DofEdge>::iterator edgeIt = periodicEdges.begin();
	 edgeIt != periodicEdges.end(); ++edgeIt) {
      vector<ElementObjectData> &edges0 = edgeElements[edgeIt->first.first];
      vector<ElementObjectData> &edges1 = edgeElements[edgeIt->first.second];

      for (unsigned int i = 0; i < edges0.size(); i++) {
	BoundaryObject obj0(macroElIndexMap[edges0[i].elIndex][macroMesh], 
			    macroElIndexTypeMap[edges0[i].elIndex], 
			    EDGE, edges0[i].ithObject);

	for (unsigned int j = 0; j < edges1.size(); j++) {
	  BoundaryObject obj1(macroElIndexMap[edges1[j].elIndex][macroMesh], 
			      macroElIndexTypeMap[edges1[j].elIndex], 
			      EDGE, edges1[j].ithObject);

	  bool reverseMode = 
	    BoundaryObject::computeReverseMode(obj0, obj1, feSpace, 
					       edgeIt->second);
	  if (reverseMode) {
	    edgeReverseMode.insert(make_pair(edges0[i], edges1[j]));
	    edgeReverseMode.insert(make_pair(edges1[j], edges0[i]));
	  }
	}
      }
    }

    for (PerBoundMap<DofFace>::iterator faceIt = periodicFaces.begin();
	 faceIt != periodicFaces.end(); ++faceIt) {
      vector<ElementObjectData> &faces0 = faceElements[faceIt->first.first];
      vector<ElementObjectData> &faces1 = faceElements[faceIt->first.second];

      TEST_EXIT_DBG((faces0.size() == 1) && (faces1.size() == 1))("Should not happen!\n");

      BoundaryObject obj0(macroElIndexMap[faces0[0].elIndex][macroMesh], 
			  macroElIndexTypeMap[faces0[0].elIndex], 
			  FACE, faces0[0].ithObject);      
      BoundaryObject obj1(macroElIndexMap[faces1[0].elIndex][macroMesh], 
			  macroElIndexTypeMap[faces1[0].elIndex], 
			  FACE, faces1[0].ithObject);
      
      bool reverseMode = 
	BoundaryObject::computeReverseMode(obj0, obj1, feSpace, faceIt->second);
      if (reverseMode) {
	faceReverseMode.insert(make_pair(faces0[0], faces1[0]));
	faceReverseMode.insert(make_pair(faces1[0], faces0[0]));
      }
    }  
  }


  int ElementObjectDatabase::getIterateOwner(int level)
  {
    FUNCNAME("ElementObjectDatabase::getIterateOwner()");

    TEST_EXIT_DBG(macroElementRankMap)("Should not happen!\n");

    switch (iterGeoPos) {
    case VERTEX:
      return getOwner(vertexElements[vertexIter->first], level);
      break;
    case EDGE:
      return getOwner(edgeElements[edgeIter->first], level);
      break;
    case FACE:
      return getOwner(faceElements[faceIter->first], level);
      break;
    default:
      ERROR_EXIT("There is something reallllly wrong!\n");
      return -1;
    }
    
    return -1;
  }


  int ElementObjectDatabase::getOwner(DegreeOfFreedom vertex, int level)
  {
    return getOwner(vertexElements[vertex], level);
  }


  int ElementObjectDatabase::getOwner(DofEdge edge, int level)
  {
    return getOwner(edgeElements[edge], level);
  }


  int ElementObjectDatabase::getOwner(DofFace face, int level)
  {
    return getOwner(faceElements[face], level);
  }


  int ElementObjectDatabase::getOwner(vector<ElementObjectData>& objData, 
				      int level)
  {
    FUNCNAME_DBG("ElementObjectDatabase::getOwner()");

    int owner = -1;

    std::set<int> &levelRanks = levelData->getLevelRanks(level);
    bool allRanks = (level == 0);
   
    for (vector<ElementObjectData>::iterator it = objData.begin();
	 it != objData.end(); ++it) {
      int elRank = (*macroElementRankMap)[it->elIndex];
      if (allRanks || levelRanks.count(elRank))
	owner = std::max(owner, elRank);
    }
    
    TEST_EXIT_DBG(owner >= 0)("Cannot find owner on level %d\n", level);
    
    return owner;
  }


  int ElementObjectDatabase::getIterateMaxLevel()
  {
    FUNCNAME("ElementObjectDatabase::getIterateMaxLevel()");

    int nLevel = levelData->getNumberOfLevels();
    TEST_EXIT_DBG(nLevel > 0)("Should not happen!\n");

    vector<std::set<int> > ranksInLevel;
    ranksInLevel.resize(nLevel);

    switch (iterGeoPos) {
    case VERTEX:
      {
	vector<ElementObjectData>& vertexData = vertexElements[vertexIter->first];
	for (vector<ElementObjectData>::iterator it = vertexData.begin();
	     it != vertexData.end(); ++it)
	  ranksInLevel[0].insert((*macroElementRankMap)[it->elIndex]);
      }
      break;
    case EDGE:
      {
	vector<ElementObjectData>& edgeData = edgeElements[edgeIter->first];
	for (vector<ElementObjectData>::iterator it = edgeData.begin();
	     it != edgeData.end(); ++it)
	  ranksInLevel[0].insert((*macroElementRankMap)[it->elIndex]);
      }
      break;
    case FACE:
      {
	vector<ElementObjectData>& faceData = faceElements[faceIter->first];
	for (vector<ElementObjectData>::iterator it = faceData.begin();
	     it != faceData.end(); ++it)
	  ranksInLevel[0].insert((*macroElementRankMap)[it->elIndex]);
      }
      break;
    default:
      ERROR_EXIT("Should not happen!\n");
    }

    for (std::set<int>::iterator it = ranksInLevel[0].begin(); 
	 it != ranksInLevel[0].end(); ++it) 
      for (int level = 1; level < nLevel; level++)
	ranksInLevel[level].insert(levelData->getLevelId(level, *it));

    int maxLevel = 0;
    for (int level = 1; level < nLevel; level++)
      if (ranksInLevel[level].size() > 1)
	maxLevel = level;
    
    return maxLevel;
  }


  void ElementObjectDatabase::serialize(ostream &out)
  {
    int nSize = vertexElements.size();
    SerUtil::serialize(out, nSize);
    for (flat_map<DegreeOfFreedom, vector<ElementObjectData> >::iterator it = vertexElements.begin();
	 it != vertexElements.end(); ++it) {
      SerUtil::serialize(out, it->first);
      serialize(out, it->second);
    }

    nSize = edgeElements.size();
    SerUtil::serialize(out, nSize);
    for (flat_map<DofEdge, vector<ElementObjectData> >::iterator it = edgeElements.begin();
	 it != edgeElements.end(); ++it) {
      SerUtil::serialize(out, it->first);
      serialize(out, it->second);
    }

    nSize = faceElements.size();
    SerUtil::serialize(out, nSize);
    for (flat_map<DofFace, vector<ElementObjectData> >::iterator it = faceElements.begin();
	 it != faceElements.end(); ++it) {
      SerUtil::serialize(out, it->first);
      serialize(out, it->second);
    }



    nSize = vertexLocalMap.size();
    SerUtil::serialize(out, nSize);
    for (boost::container::flat_map<ElementObjectData, DegreeOfFreedom>::iterator it = vertexLocalMap.begin();
	 it != vertexLocalMap.end(); ++it) {
      it->first.serialize(out);
      SerUtil::serialize(out, it->second);
    }

    nSize = edgeLocalMap.size();
    SerUtil::serialize(out, nSize);
    for (boost::container::flat_map<ElementObjectData, DofEdge>::iterator it = edgeLocalMap.begin();
	 it != edgeLocalMap.end(); ++it) {
      it->first.serialize(out);
      SerUtil::serialize(out, it->second);
    }

    nSize = faceLocalMap.size();
    SerUtil::serialize(out, nSize);
    for (boost::container::flat_map<ElementObjectData, DofFace>::iterator it = faceLocalMap.begin();
	 it != faceLocalMap.end(); ++it) {
      it->first.serialize(out);
      SerUtil::serialize(out, it->second);
    }


    nSize = vertexInRank.size();
    SerUtil::serialize(out, nSize);
    for (flat_map<DegreeOfFreedom, flat_map<int, ElementObjectData> >::iterator it = vertexInRank.begin();
	 it != vertexInRank.end(); ++it) {
      SerUtil::serialize(out, it->first);
      serialize(out, it->second);
    }

    nSize = edgeInRank.size();
    SerUtil::serialize(out, nSize);
    for (flat_map<DofEdge, flat_map<int, ElementObjectData> >::iterator it = edgeInRank.begin();
	 it != edgeInRank.end(); ++it) {
      SerUtil::serialize(out, it->first);
      serialize(out, it->second);
    }

    nSize = faceInRank.size();
    SerUtil::serialize(out, nSize);
    for (flat_map<DofFace, flat_map<int, ElementObjectData> >::iterator it = faceInRank.begin();
	 it != faceInRank.end(); ++it) {
      SerUtil::serialize(out, it->first);
      serialize(out, it->second);
    }

    SerUtil::serialize(out, periodicVertices);
    SerUtil::serialize(out, periodicEdges);
    SerUtil::serialize(out, periodicFaces);


    nSize = periodicDofAssoc.size();
    SerUtil::serialize(out, nSize);
    for (map<DegreeOfFreedom, std::set<BoundaryType> >::iterator it = periodicDofAssoc.begin();
	 it != periodicDofAssoc.end(); ++it) {
      SerUtil::serialize(out, it->first);
      SerUtil::serialize(out, it->second);
    }

    nSize = periodicEdgeAssoc.size();
    SerUtil::serialize(out, nSize);
    for (map<DofEdge, std::set<DofEdge> >::iterator it = periodicEdgeAssoc.begin();
	 it != periodicEdgeAssoc.end(); ++it) {
      SerUtil::serialize(out, it->first);
      SerUtil::serialize(out, it->second);
    }


    nSize = edgeReverseMode.size();
    SerUtil::serialize(out, nSize);
    for (std::set<pair<ElementObjectData, ElementObjectData> >::iterator it = edgeReverseMode.begin();
	 it != edgeReverseMode.end(); ++it) {
      it->first.serialize(out);
      it->second.serialize(out);
    }

    nSize = faceReverseMode.size();
    SerUtil::serialize(out, nSize);
    for (std::set<pair<ElementObjectData, ElementObjectData> >::iterator it = faceReverseMode.begin();
	 it != faceReverseMode.end(); ++it) {
      it->first.serialize(out);
      it->second.serialize(out);
      SerUtil::serialize(out, it->second);
    }
  }


  void ElementObjectDatabase::deserialize(istream &in)
  {
    int nSize;
    SerUtil::deserialize(in, nSize);
    vertexElements.clear();
    for (int i = 0; i < nSize; i++) {
      DegreeOfFreedom dof;
      vector<ElementObjectData> data;
      SerUtil::deserialize(in, dof);
      deserialize(in, data);
      vertexElements[dof] = data;
    }

    SerUtil::deserialize(in, nSize);
    edgeElements.clear();
    for (int i = 0; i < nSize; i++) {
      DofEdge edge;
      vector<ElementObjectData> data;
      SerUtil::deserialize(in, edge);
      deserialize(in, data);
      edgeElements[edge] = data;
    }

    SerUtil::deserialize(in, nSize);
    faceElements.clear();
    for (int i = 0; i < nSize; i++) {
      DofFace face;
      vector<ElementObjectData> data;
      SerUtil::deserialize(in, face);
      deserialize(in, data);
      faceElements[face] = data;
    }
   


    SerUtil::deserialize(in, nSize);
    vertexLocalMap.clear();
    for (int i = 0; i < nSize; i++) {
      ElementObjectData data;
      DegreeOfFreedom dof;
      data.deserialize(in);
      SerUtil::deserialize(in, dof);
      vertexLocalMap[data] = dof;
    }

    SerUtil::deserialize(in, nSize);
    edgeLocalMap.clear();
    for (int i = 0; i < nSize; i++) {
      ElementObjectData data;
      DofEdge edge;
      data.deserialize(in);
      SerUtil::deserialize(in, edge);
      edgeLocalMap[data] = edge;
    }

    SerUtil::deserialize(in, nSize);
    faceLocalMap.clear();
    for (int i = 0; i < nSize; i++) {
      ElementObjectData data;
      DofFace face;
      data.deserialize(in);
      SerUtil::deserialize(in, face);
      faceLocalMap[data] = face;
    }



    SerUtil::deserialize(in, nSize);
    vertexInRank.clear();
    for (int i = 0; i < nSize; i++) {
      DegreeOfFreedom dof;
      flat_map<int, ElementObjectData> data;
      SerUtil::deserialize(in, dof);
      deserialize(in, data);
      vertexInRank[dof] = data;
    }

    SerUtil::deserialize(in, nSize);
    edgeInRank.clear();
    for (int i = 0; i < nSize; i++) {
      DofEdge edge;
      flat_map<int, ElementObjectData> data;
      SerUtil::deserialize(in, edge);
      deserialize(in, data);
      edgeInRank[edge] = data;
    }

    SerUtil::deserialize(in, nSize);
    faceInRank.clear();
    for (int i = 0; i < nSize; i++) {
      DofFace face;
      flat_map<int, ElementObjectData> data;
      SerUtil::deserialize(in, face);
      deserialize(in, data);
      faceInRank[face] = data;
    }

    SerUtil::deserialize(in, periodicVertices);
    SerUtil::deserialize(in, periodicEdges);
    SerUtil::deserialize(in, periodicFaces);

    

    SerUtil::deserialize(in, nSize);
    periodicDofAssoc.clear();
    for (int i = 0; i < nSize; i++) {
      DegreeOfFreedom dof;
      std::set<DegreeOfFreedom> dofs;
      SerUtil::deserialize(in, dof);
      SerUtil::deserialize(in, dofs);
      periodicDofAssoc[dof] = dofs;
    }

    SerUtil::deserialize(in, nSize);
    periodicEdgeAssoc.clear();
    for (int i = 0; i < nSize; i++) {
      DofEdge edge;
      std::set<DofEdge> edges;
      SerUtil::deserialize(in, edge);
      SerUtil::deserialize(in, edges);
      periodicEdgeAssoc[edge] = edges;
    }



    SerUtil::deserialize(in, nSize);
    edgeReverseMode.clear();
    for (int i = 0; i < nSize; i++) {
      ElementObjectData obj0, obj1;
      obj0.deserialize(in);
      obj1.deserialize(in);

      edgeReverseMode.insert(make_pair(obj0, obj1));
    }

    SerUtil::deserialize(in, nSize);
    faceReverseMode.clear();
    for (int i = 0; i < nSize; i++) {
      ElementObjectData obj0, obj1;
      obj0.deserialize(in);
      obj1.deserialize(in);

      faceReverseMode.insert(make_pair(obj0, obj1));
    }
  }


  void ElementObjectDatabase::serialize(ostream &out, 
					vector<ElementObjectData>& elVec)
  {
    int nSize = elVec.size();
    SerUtil::serialize(out, nSize);
    for (int i = 0; i < nSize; i++)
      elVec[i].serialize(out);
  }


  void ElementObjectDatabase::deserialize(istream &in, 
					  vector<ElementObjectData>& elVec)
  {
    int nSize;
    SerUtil::deserialize(in, nSize);
    elVec.resize(nSize);
    for (int i = 0; i < nSize; i++)
      elVec[i].deserialize(in);
  }

  
  void ElementObjectDatabase::serialize(ostream &out, 
					flat_map<int, ElementObjectData>& data)
  {
    int nSize = data.size();
    SerUtil::serialize(out, nSize);
    for (flat_map<int, ElementObjectData>::iterator it = data.begin();
	 it != data.end(); ++it) {
      SerUtil::serialize(out, it->first);
      it->second.serialize(out);
    }
  }

  
  void ElementObjectDatabase::deserialize(istream &in, 
					  flat_map<int, ElementObjectData>& data)
  {
    int nSize;
    SerUtil::deserialize(in, nSize);
    for (int i = 0; i < nSize; i++) {
      int index;
      ElementObjectData elObj;
      SerUtil::deserialize(in, index);
      elObj.deserialize(in);

      data[index] = elObj;
    }
  }

  //TODO: 047 not yet reimplemented.
  unsigned long ElementObjectDatabase::calculateMemoryUsage()
  {
    FUNCNAME("ElementObjectDatabase::calculateMemoryUsage()");

    const unsigned int structElObjDataSize = sizeof(ElementObjectData);
    const unsigned int dofSize = sizeof(DegreeOfFreedom);
    const unsigned int edgeSize = sizeof(DofEdge);
    const unsigned int faceSize = sizeof(DofFace);
    const unsigned int vectorOverhead = sizeof(vector<int>);
    const unsigned int mapOverhead = 48; //sizeof(_Rb_tree<int, int>);
    const unsigned int flatMapOverhead = 24;
//     const unsigned int mapEntryOverhead = 40; // sizeof(_Rb_tree_node_base);
    const unsigned int setOverhead = 48; 
    const unsigned int setEntryOverhead = 40;


    unsigned long value = 0;
    unsigned long tmp = 0;

    // vertexElements
    tmp = flatMapOverhead;
    for (flat_map<DegreeOfFreedom, vector<ElementObjectData> >::iterator mapIt =
	   vertexElements.begin(); mapIt != vertexElements.end(); ++mapIt)      
      tmp += dofSize + vectorOverhead + mapIt->second.size() * structElObjDataSize;
    MSG("EL-OBJ-DB MEM 01: %d\n", tmp);
    value += tmp;

    // edgeElements
    tmp = flatMapOverhead;
    for (flat_map<DofEdge, vector<ElementObjectData> >::iterator mapIt =
	   edgeElements.begin(); mapIt != edgeElements.end(); ++mapIt)
      tmp += dofSize + vectorOverhead + mapIt->second.size() * structElObjDataSize;
    MSG("EL-OBJ-DB MEM 02: %d\n", tmp);
    value += tmp;

    // faceElements
    tmp = flatMapOverhead;
    for (flat_map<DofFace, vector<ElementObjectData> >::iterator mapIt =
	   faceElements.begin(); mapIt != faceElements.end(); ++mapIt)
      tmp += dofSize + vectorOverhead + mapIt->second.size() * structElObjDataSize;
    MSG("EL-OBJ-DB MEM 03: %d\n", tmp);
    value += tmp;
    tmp = 0;

    // vertexLocalMap
    tmp += flatMapOverhead + vertexLocalMap.size() * (structElObjDataSize + dofSize);

    // edgeLocalMap
    tmp += flatMapOverhead + edgeLocalMap.size() * (structElObjDataSize + edgeSize);

    // faceLocalMap
    tmp += flatMapOverhead + faceLocalMap.size() * (structElObjDataSize + faceSize);

    MSG("EL-OBJ-DB MEM 04: %d\n", tmp);
    value += tmp;

    // vertexInRank
    tmp = flatMapOverhead;
    for (flat_map<DegreeOfFreedom, flat_map<int, ElementObjectData> >::iterator mapIt = 
	   vertexInRank.begin(); mapIt != vertexInRank.end(); ++mapIt)
      tmp += dofSize + flatMapOverhead + mapIt->second.size() * (sizeof(int) + structElObjDataSize);
    MSG("EL-OBJ-DB MEM 05: %d\n", tmp);
    value += tmp;

    // edgeInRank
    tmp = mapOverhead;
    for (flat_map<DofEdge, flat_map<int, ElementObjectData> >::iterator mapIt = 
	   edgeInRank.begin(); mapIt != edgeInRank.end(); ++mapIt)
      tmp += edgeSize + flatMapOverhead + mapIt->second.size() * (sizeof(int) + structElObjDataSize);
    MSG("EL-OBJ-DB MEM 06: %d\n", tmp);
    value += tmp;

    // faceInRank
    tmp = mapOverhead;
    for (flat_map<DofFace, flat_map<int, ElementObjectData> >::iterator mapIt = 
	   faceInRank.begin(); mapIt != faceInRank.end(); ++mapIt)
      tmp += faceSize + flatMapOverhead + mapIt->second.size() * (sizeof(int) + structElObjDataSize);
    MSG("EL-OBJ-DB MEM 07: %d\n", tmp);
    value += tmp;
    tmp = 0;

    if (bConnMap.size() || periodicVertices.size() || periodicDofAssoc.size()) {
      ERROR_EXIT("Not yet implemented for periodic meshes!\n");
    }

    // edgeReverseMode
    tmp += setOverhead + edgeReverseMode.size() * (setEntryOverhead + 2 * structElObjDataSize);

    // faceReverseMode
    tmp += setOverhead + faceReverseMode.size() * (setEntryOverhead + 2 * structElObjDataSize);

    // macroElIndexMap
    tmp += flatMapOverhead + macroElIndexMap.size() * (sizeof(int) + sizeof(int*));

    // macroElIndexTypeMap
    tmp += flatMapOverhead + macroElIndexTypeMap.size() * (sizeof(int) + sizeof(int));

    MSG("EL-OBJ-DB MEM 08: %d\n", tmp);
    value += tmp;

    return value;
  }

} }
