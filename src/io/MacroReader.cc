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


#include <string.h>
#include <map>
#include <iostream>
#include <fstream>

#include "MacroReader.h"
#include "MacroWriter.h"
#include "MacroElement.h"
#include "MacroInfo.h"
#include "Boundary.h"
#include "FiniteElemSpace.h"
#include "Mesh.h"
#include "FixVec.h"
#include "ElInfo.h"
#include "Initfile.h"
#include "DOFIterator.h"
#include "LeafData.h"
#include "VertexVector.h"

namespace AMDiS { namespace io {
  
  void MacroReader::restoreMacroDofs(MacroInfo& macroInfo, int check)
  {
    FUNCNAME("MacroReader::restoreMacroDofs()");
    
    TEST_EXIT(macroInfo.isInitialized())
      ("Cannot set dof since macro info is not initialized.\n");
      
    Mesh* mesh = macroInfo.mesh;
    int nVertices = macroInfo.nVertices;
    int nElements = macroInfo.nElements;
    int **mel_vertex = macroInfo.mel_vertex;
    DegreeOfFreedom** dof = macroInfo.dof;
    
    TEST_EXIT(mesh->getNumberOfMacros() == nElements) 
      ("Mesh's current macro doesn't match the macro file.\n");
    
    mesh->setNumberOfElements(nElements);
    mesh->setNumberOfLeaves(nElements);
    mesh->setNumberOfVertices(nVertices);
    
    // Vertices dofs
    for (int i = 0; i < nVertices; i++)
      dof[i] = mesh->getDof(VERTEX);
    
    std::deque<MacroElement*>::iterator it = mesh->firstMacroElement();
    
    for (int i = 0; i < mesh->getNumberOfMacros(); i++) 
      for (int k = 0; k < mesh->getGeo(VERTEX); k++) 
	const_cast<Element*>((*(it + i))->getElement())->
	  setDof(k, dof[mel_vertex[i][k]]);
	  
    // Edge and face dofs
    if (mesh->getDim() > 1)
      boundaryDOFs(mesh);
    
    // Center dofs
    if (mesh->getNumberOfDofs(CENTER)) {
      it = mesh->firstMacroElement();
      
      for (int i = 0; i < mesh->getNumberOfMacros(); i++)
	const_cast<Element*>(it[i]->getElement())->
	  setDof(mesh->getNode(CENTER), mesh->getDof(CENTER));
    }
    
    if (check)
      checkMesh(mesh);
  }

  MacroInfo* MacroReader::readMacro(std::string filename, 
				    Mesh* mesh,
				    std::string periodicFile,
				    int check)
  {
    FUNCNAME("MacroReader::readMacro()");

    TEST_EXIT(filename != "")("no file specified; filename NULL pointer\n");  

    MacroInfo *macroInfo = new MacroInfo();
    macroInfo->readAMDiSMacro(filename, mesh);

    std::deque<MacroElement*>::iterator mel = macroInfo->mel.begin();
    int **melVertex = macroInfo->mel_vertex;
    WorldVector<double> *coords = macroInfo->coords;
    DegreeOfFreedom **dof = macroInfo->dof;

    // === read periodic data =================================
    if (periodicFile != "") {   
      FILE *file = fopen(periodicFile.c_str(), "r");
      TEST_EXIT(file)("can't open file %s\n", periodicFile.c_str());

      int n;
      int dim = mesh->getDim();

      int el1, el2;
      int *verticesEl1 = new int[dim];
      int *verticesEl2 = new int[dim];
      int mode = -1; // 0: drop dofs, 1: associate dofs
      int result = 0;
      BoundaryType boundaryType;
      MacroReader::PeriodicMap periodicMap;

      result = fscanf(file, "%*s %d", &n);
      result = fscanf(file, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
   
      for (int i = 0; i < n; i++) {
	std::map<int, int> vertexMapEl1;
	std::map<int, int> vertexMapEl2;

	result = fscanf(file, "%d", &mode);
	TEST_EXIT(result == 1)("mode?\n");
      
	result = fscanf(file, "%d", &boundaryType);
	TEST_EXIT(result == 1)("boundaryType?\n");
      
	result = fscanf(file, "%d", &el1);
	TEST_EXIT(result == 1)("el1?\n");

	for (int j = 0; j < dim; j++) {
	  result = fscanf(file, "%d", &verticesEl1[j]);
	  TEST_EXIT(result == 1)("vertEl1[%d]\n", j);
	}
	result = fscanf(file, "%d", &el2);
	TEST_EXIT(result == 1)("el2?\n");
	for (int j = 0; j < dim; j++) {
	  result = fscanf(file, "%d", &verticesEl2[j]);
	  TEST_EXIT(result == 1)("vertEl2[%d]\n", j);
	}
	for (int j = 0; j < dim; j++) {
	  if (mode == 0)
	    periodicMap.setEntry(melVertex[el1][verticesEl1[j]], 
				 melVertex[el2][verticesEl2[j]]);	  
	  vertexMapEl1[verticesEl1[j]] = verticesEl2[j];
	  vertexMapEl2[verticesEl2[j]] = verticesEl1[j];
	}

	// calculate sides of periodic vertices
	int sideEl1 = 0, sideEl2 = 0;
	if (dim == 1) {
	  sideEl1 = verticesEl1[0];
	  sideEl2 = verticesEl2[0];
	} else {
	  for (int j = 0; j < dim + 1; j++) {
	    sideEl1 += j;
	    sideEl2 += j;
	  }
	  for (int j = 0; j < dim; j++) {
	    sideEl1 -= verticesEl1[j];
	    sideEl2 -= verticesEl2[j];
	  }
	}

	// create periodic info
	DimVec<WorldVector<double> > periodicCoordsEl1(dim - 1, NO_INIT);
	DimVec<WorldVector<double> > periodicCoordsEl2(dim - 1, NO_INIT);

	Element *element1 = const_cast<Element*>((*(mel + el1))->getElement());
	Element *element2 = const_cast<Element*>((*(mel + el2))->getElement());
     
	// for all vertices of this side
	for (int j = 0; j < dim; j++) {
	  periodicCoordsEl1[element1->getPositionOfVertex(sideEl1, verticesEl1[j])] = 
	    coords[melVertex[el2][vertexMapEl1[verticesEl1[j]]]];

	  periodicCoordsEl2[element2->getPositionOfVertex(sideEl2, verticesEl2[j])] =
	    coords[melVertex[el1][vertexMapEl2[verticesEl2[j]]]];
	}
	     
	// decorate leaf data
	ElementData *ld1 = element1->getElementData();
	ElementData *ld2 = element2->getElementData();

	TEST_EXIT_DBG(ld1)
	  ("Should not happen: no element data pointer in macro element %d!\n",
	   element1->getIndex());

	TEST_EXIT_DBG(ld2)
	  ("Should not happen: no element data pointer in macro element %d!\n",
	   element2->getIndex());

	LeafDataPeriodic *ldp1 = 
	  dynamic_cast<LeafDataPeriodic*>(ld1->getElementData(PERIODIC));
	LeafDataPeriodic *ldp2 = 
	  dynamic_cast<LeafDataPeriodic*>(ld2->getElementData(PERIODIC));

	if (!ldp1) {
	  ldp1 = new LeafDataPeriodic(ld1);
	  element1->setElementData(ldp1);
	}

	if (!ldp2) {
	  ldp2 = new LeafDataPeriodic(ld2);
	  element2->setElementData(ldp2);
	}

	ldp1->addPeriodicInfo(mode, boundaryType, sideEl1, &periodicCoordsEl1);
	ldp2->addPeriodicInfo(mode, boundaryType, sideEl2, &periodicCoordsEl2);

	if (mode != 0) {
	  VertexVector *associated = mesh->periodicAssociations[boundaryType];

	  if (!associated) {
	    associated = new VertexVector(mesh->getVertexAdmin(), "vertex vector");
	    mesh->periodicAssociations[boundaryType] = associated;
	    VertexVector::Iterator it(associated, ALL_DOFS);
	    for (it.reset2(); !it.end(); ++it)
	      *it = it.getDOFIndex();
	  }

	  for (int j = 0; j < dim; j++) {
#ifdef DEBUG
	  {
		  unsigned initData(melVertex[el1][verticesEl1[j]]);
		  unsigned oldData((*associated)[melVertex[el1][verticesEl1[j]]]);
		  unsigned newData(melVertex[el2][vertexMapEl1[verticesEl1[j]]]);
		  if( initData != oldData && newData != oldData ) {
			  MSG("warning: element %d overwrites assoc index %d: %d -> %d\n", 
					  el1, initData, oldData, newData);
		  }
	  }
#endif
	    (*associated)[melVertex[el1][verticesEl1[j]]] =
	      melVertex[el2][vertexMapEl1[verticesEl1[j]]];
	    (*associated)[melVertex[el2][verticesEl2[j]]] =
	      melVertex[el1][vertexMapEl2[verticesEl2[j]]];
	  }
	}
      }    

      delete [] verticesEl1;
      delete [] verticesEl2;

      // change periodic vertex dofs
      for (int i = 0; i < mesh->getNumberOfVertices(); i++) {
	if (periodicMap.getEntry(i) != -1) {
	  mesh->freeDof(dof[i], VERTEX);
	  dof[i] = dof[periodicMap.getEntry(i)];

	  std::map<BoundaryType, VertexVector*>::iterator assoc;
	  std::map<BoundaryType, VertexVector*>::iterator assocEnd =
	    mesh->periodicAssociations.end();

	  for (assoc = mesh->periodicAssociations.begin(); assoc != assocEnd; ++assoc) {

	    DegreeOfFreedom a = (*(assoc->second))[i];
	    if (a != i) {
	      (*(assoc->second))[i] = i;
	      (*(assoc->second))[a] = periodicMap.getEntry(i);
	    }
	  }

	}
      }

#if (DEBUG != 0)
      std::map<BoundaryType, VertexVector*>::iterator assoc;
      std::map<BoundaryType, VertexVector*>::iterator assocEnd =
	mesh->periodicAssociations.end();
      for (assoc = mesh->periodicAssociations.begin(); 
	   assoc != assocEnd; ++assoc)
	for (int i = 0; i < mesh->getNumberOfVertices(); i++)
	  if (i != (*(assoc->second))[i])
	    MSG("association %d: vertex %d -> vertex %d\n", 
		assoc->first, i, (*(assoc->second))[i]);

      for (int i = 0; i < mesh->getNumberOfVertices(); i++)
	if (periodicMap.getEntry(i) != -1)
	  MSG("identification : vertex %d is now vertex %d\n", 
	      i, periodicMap.getEntry(i));
#endif

    } // periodicFile


    // =========================================================
    // Set coordinate
    for (int i = 0; i < mesh->getNumberOfMacros(); i++) {
      for (int k = 0; k < mesh->getGeo(VERTEX); k++) {
	(*(mel + i))->setCoord(k, coords[melVertex[i][k]]);

	const_cast<Element*>((*(mel + i))->getElement())->
	  setDof(k, dof[melVertex[i][k]]);
      }
    }

    if (!macroInfo->neigh_set) {
      TEST_EXIT(periodicFile == "")
	("periodic boundary condition => element neighbours must be set\n");
      computeNeighbours(mesh);
    } else {
      /****************************************************************************/
      /* fill MEL oppVertex values when reading neighbour information form file  */
      /****************************************************************************/

      for (int i = 0; i < mesh->getNumberOfMacros(); i++) {
	for (int k = 0; k < mesh->getGeo(NEIGH); k++) {
	  MacroElement *neigh = const_cast<MacroElement*>(mel[i]->getNeighbour(k));

	  if (neigh) {
	    if (!macroInfo->neigh_inverse_set) {
	      int j = 0;
	      for (; j < mesh->getGeo(NEIGH); j++)
		if (neigh->getNeighbour(j) == *(mel + i))  
		  break;
	  
	      TEST_EXIT(j < mesh->getGeo(NEIGH))("el %d no neighbour of neighbour %d\n", 
						mel[i]->getIndex(), neigh->getIndex());
	      mel[i]->setOppVertex(k, j);
	    } else {
	      int j = 0;
	      for (; j < mesh->getGeo(NEIGH); j++)
		if (neigh->getNeighbourInv(j) == *(mel + i))  
		  break;
	  
	      TEST_EXIT(j < mesh->getGeo(NEIGH))("el %d no neighbour-inv of neighbour %d\n", 
						mel[i]->getIndex(), neigh->getIndex());
	      mel[i]->setOppVertex(k, j);
	    }
	  } else {
	    mel[i]->setOppVertex(k, -1);
	  }
	}
      }
    }

    if (!macroInfo->bound_set)
      macroInfo->dirichletBoundary();    
  
    if (mesh->getDim() > 1)
      boundaryDOFs(mesh);

    // initial boundary projections
    int numFaces = mesh->getGeo(FACE);
    int dim = mesh->getDim();
    mel = mesh->firstMacroElement();
    for (int i = 0; i < mesh->getNumberOfLeaves(); i++) {
      MacroElement *macroEl = *(mel + i);
      Projection *projector = macroEl->getProjection(0);
      if (projector && projector->getType() == VOLUME_PROJECTION) {
	for (int j = 0; j <= dim; j++)
	  projector->project(macroEl->getCoord(j));	
      } else {
	for (int j = 0; j < mesh->getGeo(EDGE); j++) {
	  projector = macroEl->getProjection(numFaces + j);
	  if (projector) {
	    int vertex0 = Global::getReferenceElement(dim)->getVertexOfEdge(j, 0);
	    int vertex1 = Global::getReferenceElement(dim)->getVertexOfEdge(j, 1);
	    projector->project(macroEl->getCoord(vertex0));
	    projector->project(macroEl->getCoord(vertex1));
	  }
	}
      }
    }

    macroInfo->fillBoundaryInfo(mesh);

    if (mesh->getNumberOfDofs(CENTER)) {
      for (int i = 0; i < mesh->getNumberOfMacros(); i++)
	const_cast<Element*>(mel[i]->getElement())->
	  setDof(mesh->getNode(CENTER), mesh->getDof(CENTER));
    }

    // === Domain size ===

    WorldVector<double> x_min, x_max;

    for (int j = 0; j < Global::getGeo(WORLD); j++) {
      x_min[j] =  1.E30;
      x_max[j] = -1.E30;
    }

    for (int i = 0; i < mesh->getNumberOfVertices(); i++) {
      for (int j = 0; j < Global::getGeo(WORLD); j++) {
	x_min[j] = std::min(x_min[j], coords[i][j]);
	x_max[j] = std::max(x_max[j], coords[i][j]);
      }
    }

    for (int j = 0; j < Global::getGeo(WORLD); j++)
      mesh->setDiameter(j, x_max[j] - x_min[j]);

    if (check) {
      checkMesh(mesh);

      if (mesh->getDim() > 1)
	macroTest(mesh);
    }

    return macroInfo;
  }
  
  void MacroReader::computeNeighbours(Mesh *mesh)
  {
    FUNCNAME("MacroReader::computeNeighbours()");

    int dim = mesh->getDim();
    FixVec<DegreeOfFreedom*, DIMEN> dof(dim, NO_INIT);

    for (int i = 0; i < mesh->getNumberOfLeaves(); i++) {
      for (int k = 0; k < mesh->getGeo(NEIGH); k++) {
	mesh->getMacroElement(i)->setOppVertex(k, AMDIS_UNDEFINED);
	mesh->getMacroElement(i)->setNeighbour(k, NULL);
      }
    }

    for (int i = 0; i < mesh->getNumberOfLeaves(); i++) {
      for (int k = 0; k < mesh->getGeo(NEIGH); k++) {
	if (mesh->getMacroElement(i)->getBoundary(k) != INTERIOR) {
	  mesh->getMacroElement(i)->setNeighbour(k, NULL);
	  mesh->getMacroElement(i)->setOppVertex(k, -1);
	  continue;
	}

	if (mesh->getMacroElement(i)->getOppVertex(k) == AMDIS_UNDEFINED) {
	  if (dim == 1) {
	    dof[0] = const_cast<DegreeOfFreedom*>(mesh->getMacroElement(i)->
						  getElement()->getDof(k));
	  } else {
	    for (int l = 0; l < dim; l++)
	      dof[l] = const_cast<DegreeOfFreedom*>(mesh->getMacroElement(i)->
						    getElement()->
						    getDof((k + l + 1) % (dim + 1)));
	  }
	  
	  int j = 0;
	  for (j = i + 1; j < mesh->getNumberOfLeaves(); j++) {
	    int m = mesh->getMacroElement(j)->getElement()->oppVertex(dof);
	    if (m != -1) {
	      mesh->getMacroElement(i)->setNeighbour(k, mesh->getMacroElement(j));
	      mesh->getMacroElement(j)->setNeighbour(m, mesh->getMacroElement(i));
	      mesh->getMacroElement(i)->setOppVertex(k, m);
	      mesh->getMacroElement(j)->setOppVertex(m, k);
	      break;
	    }
	  }

	  if (j >= mesh->getNumberOfLeaves()) {
	    std::cout << "----------- ERROR ------------" << std::endl;
	    std::cout << "Cannot find neighbour " << k << " of element " << i << std::endl;
	    std::cout << "  dim = " << dim << std::endl;
	    std::cout << "  coords of element = ";
	    for (int l = 0; l <= dim; l++) {
	      std::cout << mesh->getMacroElement(i)->getCoord(l);
	      if (l < dim)
		std::cout << " / ";	      
	    }
	    std::cout << std::endl;
	    std::cout << "  dofs = ";
	    for (int l = 0; l < dim; l++)
	      std::cout << *(dof[l]) << " ";	    
	    std::cout << std::endl;

	    ERROR_EXIT("\n");
	  }    
	}
      }
    }
  }


  /****************************************************************************/
  /*  boundaryDOFs:                                                           */
  /*  adds dof's at the edges of a given macro triangulation and calculates   */
  /*  the number of edges                                                     */
  /****************************************************************************/

  void MacroReader::boundaryDOFs(Mesh *mesh)
  {
    FUNCNAME("Mesh::boundaryDOFs()");

    int lnode = mesh->getNode(EDGE);
    int k, lne = mesh->getNumberOfLeaves();
    int max_n_neigh = 0, n_neigh, ov;
    std::deque<MacroElement*>::iterator mel;
    const MacroElement* neigh;
    DegreeOfFreedom *dof;

    mesh->setNumberOfEdges(0);
    mesh->setNumberOfFaces(0);

    int dim = mesh->getDim();

    switch (dim) {
    case 2:
      for (mel = mesh->firstMacroElement(); mel != mesh->endOfMacroElements(); mel++) {
	// check for periodic boundary
	Element *el = const_cast<Element*>((*mel)->getElement());
	ElementData *ed = el->getElementData(PERIODIC);

	DimVec<bool> periodic(dim, DEFAULT_VALUE, false);

	if (ed) {
	  std::list<LeafDataPeriodic::PeriodicInfo> &periodicInfos = 
	    dynamic_cast<LeafDataPeriodic*>(ed)->getInfoList();
	  std::list<LeafDataPeriodic::PeriodicInfo>::iterator it;
	  std::list<LeafDataPeriodic::PeriodicInfo>::iterator end = periodicInfos.end();
	  for (it = periodicInfos.begin(); it != end; ++it)
	    if (it->type != 0)
	      periodic[it->elementSide] = true;
	}

	for (int i = 0; i < mesh->getGeo(NEIGH); i++) {
	  if (!(*mel)->getNeighbour(i) || 
	      ((*mel)->getNeighbour(i)->getIndex() < (*mel)->getIndex())) {

	    mesh->incrementNumberOfEdges(1);

	    if (mesh->getNumberOfDofs(EDGE)) {
	      dof = mesh->getDof(EDGE);
	      el->setDof(lnode + i, dof);
      
	      if ((*mel)->getNeighbour(i)) {
		Element *neigh = 
		  const_cast<Element*>((*mel)->getNeighbour(i)->getElement());

		if (periodic[i])
		  neigh->setDof(lnode + (*mel)->getOppVertex(i), mesh->getDof(EDGE));
		else
		  neigh->setDof(lnode + (*mel)->getOppVertex(i), dof);		
	      }
	    }
	  }  
	}
      }
      break;
    case 3:
      lnode = mesh->getNode(FACE);
      mel = mesh->firstMacroElement();
      for (int i = 0; i < lne; i++) {

	// check for periodic boundary
	Element *el = const_cast<Element*>((*(mel+i))->getElement());
	
	ElementData *ed = el->getElementData(PERIODIC);

	DimVec<bool> periodic(dim, DEFAULT_VALUE, false);
      
	if (ed) {
	  std::list<LeafDataPeriodic::PeriodicInfo> &periodicInfos = 
	    dynamic_cast<LeafDataPeriodic*>(ed)->getInfoList();
	  std::list<LeafDataPeriodic::PeriodicInfo>::iterator it;
	  std::list<LeafDataPeriodic::PeriodicInfo>::iterator end = periodicInfos.end();
	  for (it = periodicInfos.begin(); it != end; ++it)
	    if (it->type != 0)
	      periodic[it->elementSide] = true;
	}

	for (k = 0; k < mesh->getGeo(EDGE); k++) {
	  // === Check for not counted edges. ===
	  n_neigh = 1;

	  if (newEdge(mesh, (*(mel + i)), k, &n_neigh)) {
	    mesh->incrementNumberOfEdges(1);
	    max_n_neigh = std::max(max_n_neigh, n_neigh);
	  }
	}
      
	for (k = 0; k < mesh->getGeo(NEIGH); k++) {
	  neigh = (*(mel + i))->getNeighbour(k);
	  // === Face is counted and dof is added by the element with bigger index. ===
	  if (neigh && (neigh->getIndex() > (*(mel + i))->getIndex()))  
	    continue;
	
	  mesh->incrementNumberOfFaces(1);
	
	  if (mesh->getNumberOfDofs(FACE)) {
	    TEST_EXIT(!(*(mel + i))->getElement()->getDof(lnode + k))
	      ("dof %d on element %d already set\n", 
	       lnode + k, (*(mel + i))->getIndex());
	  
	    const_cast<Element*>((*(mel + i))->getElement())->setDof(lnode + k,
								     mesh->getDof(FACE));

	    if (neigh) {
	      ov = (*(mel + i))->getOppVertex(k);
	      TEST_EXIT(!neigh->getElement()->getDof(lnode + ov))
		("dof %d on neighbour %d already set\n", 
		 lnode + ov, neigh->getIndex());
	    
	      Element *neighEl = 
		const_cast<Element*>((*(mel + i))->getNeighbour(k)->getElement());

	      if (periodic[k])
		neighEl->setDof(lnode+ov, mesh->getDof(FACE));
	      else
		neighEl->setDof(lnode+ov, const_cast<DegreeOfFreedom*>((*(mel + i))->getElement()->
							   getDof(lnode + k)));	      
	    }
	  }
	}
      }
      break;
    default: 
      ERROR_EXIT("invalid dim\n");
    }
    
    if (3 == dim)
      mesh->setMaxEdgeNeigh(std::max(8, 2 * max_n_neigh));
    else
      mesh->setMaxEdgeNeigh(dim - 1);    
  }

  /* 
     testet mesh auf auftretende Zyklen
  
     wenn Zyklus auftritt:
     ordnet Eintraege in MacroElement-Struktur um, so dass kein Zyklus auftritt
     erzeugt neue Macro-Datei nameneu mit umgeordnetem Netz 
     (wenn nameneu=NULL wird keine MAcro-Datei erzeugt)
  */      

  void MacroReader::macroTest(Mesh *mesh)
  {
    FUNCNAME("MacroReader::macroTest()");
   
    int i = macrotest(mesh);

    if (i >= 0) {
      WARNING("There is a cycle beginning in macro element %d\n", i);
      WARNING("Entries in MacroElement structures get reordered\n");
      umb(NULL, mesh, umbVkantMacro);
    }
  }
  
  /****************************************************************************/
  /*  macro_test():                              Author: Thomas Kastl (1998)  */
  /****************************************************************************/
  /*
    testet mesh auf auftretende Zyklen
  
    wenn mesh zyklenfrei -> -1
    sonst ->  globaler Index des Macroelementes bei dem ein Zyklus beginnt 
  */

  int MacroReader::macrotest(Mesh *mesh)
  {
    std::deque<MacroElement*>::const_iterator macro, mac;
    int flg = 0;
    int dim = mesh->getDim();
    int *test = new int[mesh->getNumberOfMacros()];
    int *zykl = new int[mesh->getNumberOfMacros()];
 
    for (int i = 0; i < mesh->getNumberOfMacros(); i++)
      test[i] = 0;

    int zykstart = -1;
    std::deque<MacroElement*>::const_iterator macrolfd = mesh->firstMacroElement();

    while (macrolfd != mesh->endOfMacroElements()) {
      if (test[(*macrolfd)->getIndex()] == 1) {
	macrolfd++;
      } else {
	for (int i = 0; i < mesh->getNumberOfMacros(); i++)
	  zykl[i] = 0;	
    
	macro = macrolfd;
	flg = 2;
	do {
	  if (zykl[(*macro)->getIndex()] == 1) {
	    flg = 0;
	    zykstart = (*macro)->getIndex();
	  } else {
	    zykl[(*macro)->getIndex()] = 1;
      
	    if (test[(*macro)->getIndex()] == 1) {
	      flg = 1;
	    } else if ((*macro)->getNeighbour(dim) == NULL) {
	      flg = 1;
	      test[(*macro)->getIndex()] = 1;
	    } else if ((*macro) == (*macro)->getNeighbour(dim)->getNeighbour(dim)) {
	      flg = 1;
	      test[(*macro)->getIndex()] = 1;
	      test[(*macro)->getNeighbour(dim)->getIndex()] = 1;
	    } else {
	      for (mac = mesh->firstMacroElement(); 
		   (*mac) != (*macro)->getNeighbour(dim); mac++);
	      macro = mac;
	    } 
	  }	  
	} while(flg == 2);
 
	if (flg == 1)
	  macrolfd++;
	else 
	  macrolfd=mesh->endOfMacroElements();	
      }
    }
  
    delete [] zykl;
    delete [] test;
 
    return zykstart;
  }

  //   waehlt geeignete Verfeinerungskanten, so dass kein Zyklus auftritt (recumb)

  //   ele     Integer-Vektor der Dimension Anzahl der Macro-Elemente
  //           zur Speicherung der neuen Verfeinerungskanten
  //           (wird nur benoetigt, wenn umbvk=umb_vkant_macrodat) 
  
  //   umbvk   Fkt. zur Neuordnung der Verfeinerungskanten
  //           = umb_vkant_macro :
  //               Eintraege in MacroElement-Struktur und Eintraege in macro->el
  //               werden tatsaechlich umgeordnet
  //               -> ALBERT-Routine write_macro kann zum erzeugen einer
  //                  neuen Macro-Datei angewendet werden 
  //           = umb_vkant_macrodat :
  //               Eintraege werden nicht veraendert, es werden nur die lokalen
  //               Indices der Kanten, die zu Verfeinerungskanten werden im
  //               Integer-Vektor ele abgespeichert
  //               -> print_Macrodat zur Erzeugung einer zyklenfreien neuen
  //                  Macro-Datei kann angewendet werden

  void MacroReader::umb(int *ele, Mesh *mesh,
			void (*umbvk)(Mesh*, MacroElement*, int, int*))
  {
    int *test = new int[mesh->getNumberOfMacros()];
  
    for (int i = 0; i < static_cast<int>(mesh->getNumberOfMacros()); i++)
      test[i] = 0;

    recumb(mesh, (*mesh->firstMacroElement()), NULL, test, 0, 0, ele, umbvk);

    delete [] test;
  }

  bool MacroReader::newEdge(Mesh *mesh, MacroElement *mel,
			    int mel_edge_no, int *n_neigh)
  {
    FUNCNAME("MacroElement::newEdge()"); 
    MacroElement *nei;
    const DegreeOfFreedom *dof[2];
    DegreeOfFreedom *edge_dof = NULL;
    int j, k, opp_v, node = 0;
    BoundaryType lbound = INTERIOR;
    Projection *lproject = NULL;
    const int max_no_list_el = 100;
    BoundaryType *list_bound[100];
    Projection **list_project[100];
    Element *el = const_cast<Element*>(mel->getElement());
    int edge_no = mel_edge_no;
    static int next_el[6][2] = {{3,2},{1,3},{1,2},{0,3},{0,2},{0,1}};
    int vertices = mesh->getGeo(VERTEX);

    int mel_index = mel->getIndex();

    list_bound[0] = &(mel->boundary[mesh->getGeo(FACE)+edge_no]);
    list_project[0] = &(mel->projection[mesh->getGeo(FACE)+edge_no]);

    if (mesh->getNumberOfDofs(EDGE)) {
      node = mesh->getNode(EDGE);
      if (el->getDof(node+edge_no)) {
	/****************************************************************************/
	/*  edge was counted by another macro element and dof was added on the      */
	/*  complete patch                                                          */
	/****************************************************************************/
	return false;
      } else {
	edge_dof = mesh->getDof(EDGE);
	el->setDof(node+edge_no, edge_dof);
      }
    }

    for (j = 0; j < 2; j++)
      dof[j] = el->getDof(el->getVertexOfEdge(edge_no, j));
    
    /****************************************************************************/
    /*  first look for all neighbours in one direction until a boundary is      */
    /*  reached :-( or we are back to mel :-)                                   */
    /*  if the index of a neighbour is smaller than the element index, the edge */
    /*  is counted by this neighbour, return 0.                                 */
    /*  If we are back to element, return 1, to count the edge                  */
    /****************************************************************************/

    nei = mel->getNeighbour(next_el[edge_no][0]);
    opp_v = mel->getOppVertex(next_el[edge_no][0]);


    if (mel->getBoundary(next_el[edge_no][0])) {
      lbound = newBound(mel->getBoundary(next_el[edge_no][0]), lbound);
      lproject = mel->getProjection(next_el[edge_no][0]);
    }

    while (nei  &&  nei != mel) {
      for (j = 0; j < vertices; j++)
	if (nei->getElement()->getDof(j) == dof[0])  break;
      for (k = 0; k < vertices; k++)
	if (nei->getElement()->getDof(k) == dof[1])  break;

      // check for periodic boundary
      if (j == 4 || k == 4) {
	nei = NULL;
	break;
      }

      if (mesh->getNumberOfDofs(EDGE)) {
	TEST_EXIT(nei->index > mel_index)
	  ("neighbour index %d < element index %d\n", nei->getIndex(), mel_index);
      }

      if (!mesh->getNumberOfDofs(EDGE) && nei->getIndex() < mel_index) 
	return false;
    
      edge_no = Tetrahedron::edgeOfDofs[j][k];

      TEST_EXIT(*n_neigh < max_no_list_el)
	("too many neigbours for local list\n");

      list_bound[(*n_neigh)] = 
	&(nei->boundary[mesh->getGeo(FACE)+edge_no]);

      list_project[(*n_neigh)++] = 
	&(nei->projection[mesh->getGeo(FACE)+edge_no]);

      if (mesh->getNumberOfDofs(EDGE))
	nei->element->setDof(node+edge_no,edge_dof);      

      if (next_el[edge_no][0] != opp_v) {
	if (nei->getBoundary(next_el[edge_no][0])) {
	  lbound = newBound(nei->getBoundary(next_el[edge_no][0]), lbound);
	  Projection *neiProject = nei->getProjection(next_el[edge_no][0]);
	  if (!lproject) {
	    lproject = neiProject;
	  } else {
	    if (neiProject && (lproject->getID() < neiProject->getID()))
	      lproject = neiProject;
	  }
	}
	opp_v = nei->getOppVertex(next_el[edge_no][0]);
	nei = nei->getNeighbour(next_el[edge_no][0]);
      } else {
	if (nei->getBoundary(next_el[edge_no][1])) {
	  lbound = newBound(nei->getBoundary(next_el[edge_no][1]), lbound);
	  Projection *neiProject = nei->getProjection(next_el[edge_no][1]);	
	  if (!lproject) {
	    lproject = neiProject;
	  } else {
	    if (neiProject && (lproject->getID() < neiProject->getID()))
	      lproject = neiProject;
	  }
	}
	opp_v = nei->getOppVertex(next_el[edge_no][1]);
	nei = nei->getNeighbour(next_el[edge_no][1]);
      }
    }

    if (!nei) {
      /****************************************************************************/
      /*  while looping around the edge the domain's boundary was reached. Now,   */
      /*  loop in the other direction to the domain's boundary		    */
      /****************************************************************************/
      edge_no = mel_edge_no;
      
      nei = mel->getNeighbour(next_el[edge_no][1]);
      opp_v = mel->getOppVertex(next_el[edge_no][1]);
      if (mel->getBoundary(next_el[edge_no][1])) {
	lbound = newBound(mel->getBoundary(next_el[edge_no][1]), lbound); 
	Projection *neiProject =  mel->getProjection(next_el[edge_no][1]);
	if (!lproject) {
	  lproject = neiProject;
	} else {
	  if (neiProject && (lproject->getID() < neiProject->getID()))
	    lproject = neiProject;	  
	}
      }
      
      while (nei) {
	for (j = 0; j < vertices; j++)
	  if (nei->getElement()->getDof(j) == dof[0])  break;
	for (k = 0; k < vertices; k++)
	  if (nei->getElement()->getDof(k) == dof[1])  break;
	
	// check for periodic boundary
	if (j == 4 || k == 4)
	  return false;
	
	if (mesh->getNumberOfDofs(EDGE)) {
	  TEST_EXIT(nei->getIndex() > mel_index)
	    ("neighbour index %d < element index %d\n", nei->getIndex(),
	     mel_index);
	}
	
	if (nei->getIndex() < mel_index)  
	  return false;
	
	edge_no = Tetrahedron::edgeOfDofs[j][k];
	
	TEST_EXIT(*n_neigh < max_no_list_el)("too many neigbours for local list\n");
	
	list_bound[(*n_neigh)] = &(nei->boundary[mesh->getGeo(FACE) + edge_no]);
	list_project[(*n_neigh)++] = &(nei->projection[mesh->getGeo(FACE) + edge_no]);

	if (mesh->getNumberOfDofs(EDGE)) {
	  TEST_EXIT(!nei->getElement()->getDof(node+edge_no))
	    ("dof %d on element %d is already set, but not on element %d\n",
	     node + edge_no, nei->getIndex(), mel_index);
	  
	  nei->element->setDof(node+edge_no, edge_dof);
	}

	if (next_el[edge_no][0] != opp_v) {
	  if (nei->getBoundary(next_el[edge_no][0])) {
	    lbound = newBound(nei->getBoundary(next_el[edge_no][0]), lbound);
	    Projection *neiProject = nei->getProjection(next_el[edge_no][0]);
	    if (!lproject) {
	      lproject = neiProject;
	    } else {
	      if (neiProject &&( lproject->getID() < neiProject->getID()))
		lproject = neiProject;	      
	    }
	  }
	  
	  opp_v = nei->getOppVertex(next_el[edge_no][0]);
	  nei = nei->getNeighbour(next_el[edge_no][0]);
	} else {
	  if (nei->getBoundary(next_el[edge_no][1])) {
	    lbound = newBound(nei->getBoundary(next_el[edge_no][1]), lbound); 
	    Projection *neiProject = nei->getProjection(next_el[edge_no][1]);
	    if (!lproject) {
	      lproject = neiProject;
	    } else {
	      if (neiProject && (lproject->getID() < neiProject->getID()))
		lproject = neiProject;	      
	    }
	  }
	  
	  opp_v = nei->getOppVertex(next_el[edge_no][1]);
	  nei = nei->getNeighbour(next_el[edge_no][1]);
	}
      }
    }
    
    for (j = 0; j < *n_neigh; j++) {
      *(list_bound[j]) = lbound;
      *(list_project[j]) = lproject;
    }
  
    return true;
  }

  void MacroReader::fillMelBoundary(Mesh *mesh, MacroElement *mel, 
				    FixVec<BoundaryType,NEIGH> const& ind)
  {
    for (int i = 0; i < mesh->getGeo(NEIGH); i++)
      mel->boundary[i] = ind[i];        
  }


  void MacroReader::fillMelNeigh(MacroElement *mel,
				 std::deque<MacroElement*>& macro_elements, 
				 FixVec<int,NEIGH> const& ind)
  {
    int dim = mel->element->getMesh()->getDim();

    for (int k = 0; k < Global::getGeo(NEIGH, dim); k++) {
      if (ind[k] >= 0) 
	mel->neighbour[k] = macro_elements[ind[k]];
      else
	mel->neighbour[k] = NULL;
    }
  }


  void MacroReader::fillMelNeighInv(MacroElement *mel,
				 std::deque<MacroElement*>& macro_elements, 
				 FixVec<int,NEIGH> const& ind)
  {
    int dim = mel->element->getMesh()->getDim();

    for (int k = 0; k < Global::getGeo(NEIGH, dim); k++) {
      if (ind[k] >= 0) 
	mel->neighbour_inv[k] = macro_elements[ind[k]];
      else
	mel->neighbour_inv[k] = NULL;
    }
  }


  //   ordnet Eintraege in Macro-Element macro bzgl. Verfeinerungskante ka um
  //   (coord, bound, boundary, neigh, oppVertex)

  //   ordnet Eintraege in macro->el bzgl. Verfeinerungskante ka um
  //   (Element-DOF's (macro->el->dof) in Ecken und auf Kanten,
  //    wenn NEIGH_IN_EL macro->el->neigh, macro->el->oppVertex)
  //   (wird fuer ALBERT-Routine write_macro benoetigt)

  //   ele wird nicht benoetigt (es kann NULL uebergeben werden)    
  void MacroReader::umbVkantMacro(Mesh *mesh, MacroElement* me, int ka, int *)
  {
    MacroElement* macr=new MacroElement(mesh->getDim());
    int i;
    int n0;
    DegreeOfFreedom *d[7];
  
    int vertices = mesh->getGeo(VERTEX);
    int facesPlusEdges = mesh->getGeo(EDGE) + mesh->getGeo(FACE);

    if (ka == 2);
    else { 
      for (i = 0; i < 3; i++) {
	macr->coord[i]=me->coord[i];
	macr->setBoundary(facesPlusEdges + i, me->getBoundary(facesPlusEdges + i));
	macr->setBoundary(i, me->getBoundary(i));
	macr->setNeighbour(i, me->getNeighbour(i));
	macr->setOppVertex(i,me->getOppVertex(i));
      }    
  
      for (i = 0; i < 7; i++)
	d[i] = const_cast<DegreeOfFreedom*>(me->getElement()->getDof(i));      

      if (ka == 1) { 
	me->coord[0] = macr->coord[2];
	me->coord[1] = macr->coord[0];
	me->coord[2] = macr->coord[1];

	me->setBoundary(facesPlusEdges + 0,macr->getBoundary(facesPlusEdges + 2));
	me->setBoundary(facesPlusEdges + 1,macr->getBoundary(facesPlusEdges + 0));
	me->setBoundary(facesPlusEdges + 2,macr->getBoundary(facesPlusEdges + 1));

	me->setBoundary(0, macr->getBoundary(2));
	me->setBoundary(1, macr->getBoundary(0));
	me->setBoundary(2, macr->getBoundary(1));

	me->setNeighbour(0,const_cast<MacroElement*>(macr->getNeighbour(2)));
	me->setNeighbour(1,const_cast<MacroElement*>(macr->getNeighbour(0)));
	me->setNeighbour(2,const_cast<MacroElement*>(macr->getNeighbour(1)));

	me->setOppVertex(0,macr->getOppVertex(2));
	me->setOppVertex(1,macr->getOppVertex(0));
	me->setOppVertex(2,macr->getOppVertex(1));


	if (mesh->getNumberOfDofs(VERTEX)) {                /* Ecken */
	  n0 = mesh->getNode(VERTEX);              
        
	  const_cast<Element*>(me->getElement())->setDof(n0,d[n0+2]);     
	  const_cast<Element*>(me->getElement())->setDof(n0+1,d[n0]);  
	  const_cast<Element*>(me->getElement())->setDof(n0+2,d[n0+1]);   
	}
 
	if (mesh->getNumberOfDofs(EDGE)) {                  /* Kanten */
	  n0 = mesh->getNode(EDGE);    
       
	  const_cast<Element*>(me->getElement())->setDof(n0,d[n0+2]);  
	  const_cast<Element*>(me->getElement())->setDof(n0+1,d[n0]);  
	  const_cast<Element*>(me->getElement())->setDof(n0+2,d[n0+1]);
	} 
      } else {
	me->coord[0] = macr->coord[1];
	me->coord[1] = macr->coord[2];
	me->coord[2] = macr->coord[0];

	me->setBoundary(facesPlusEdges + 0,macr->getBoundary(facesPlusEdges + 1));
	me->setBoundary(facesPlusEdges + 1,macr->getBoundary(facesPlusEdges + 2));
	me->setBoundary(facesPlusEdges + 2,macr->getBoundary(facesPlusEdges + 0));

	me->setBoundary(0, macr->getBoundary(1));
	me->setBoundary(1, macr->getBoundary(2));
	me->setBoundary(2, macr->getBoundary(0));

	me->setNeighbour(0,const_cast<MacroElement*>(macr->getNeighbour(1)));
	me->setNeighbour(1,const_cast<MacroElement*>(macr->getNeighbour(2)));
	me->setNeighbour(2,const_cast<MacroElement*>(macr->getNeighbour(0)));

	me->setOppVertex(0,macr->getOppVertex(1));
	me->setOppVertex(1,macr->getOppVertex(2));
	me->setOppVertex(2,macr->getOppVertex(0));
    
	if (mesh->getNumberOfDofs(VERTEX)) {                /* Ecken */
	  n0 = mesh->getNode(VERTEX);              
        
	  const_cast<Element*>(me->getElement())->setDof(n0,d[n0+1]);     
	  const_cast<Element*>(me->getElement())->setDof(n0+1,d[n0+2]);  
	  const_cast<Element*>(me->getElement())->setDof(n0+2,d[n0]);   
	}
 
	if (mesh->getNumberOfDofs(EDGE)) {                  /* Kanten */
	  n0 = mesh->getNode(EDGE);    
       
	  const_cast<Element*>(me->getElement())->setDof(n0,d[n0+1]);  
	  const_cast<Element*>(me->getElement())->setDof(n0+1,d[n0+2]);  
	  const_cast<Element*>(me->getElement())->setDof(n0+2,d[n0]);
	} 
      }
  

      for (i = 0; i < vertices; i++)   /* oppVertex der Nachbarn umsetzen*/  
	if (me->getNeighbour(i))
	  const_cast<MacroElement*>(me->getNeighbour(i))->setOppVertex(me->getOppVertex(i), i);
    }
    delete macr;
  }


  //   durchlaeuft rek. Macro-Triangulierung mit Start auf Macro-Element macro und
  //   waehlt geignete Verfeinerungskanten, so dass kein Zyklus auftritt
  
  //   (Umbennenung der lokalen Indices der Macro-Elemente
  //    oder Speichern der neuen Verfeinerungskanten mit Fkt. umbvk)
  
  //   waehlt als neue Verfeinerungskante die laengste Kante eines Elementes
  //   fuehrt "fiktiv verlaengerte" Kanten ein, 
  //   wenn ein Element 2 laengste Kanten besitzt  
 
  //   macroalt  Zeiger auf "Rekursions-Vorgaenger"-Macroelement
  //   lg        Laenge der Verfeinerungskante von "Vorgaenger"
  //   ka = 1, wenn Verfeinerungskante von "Vorgaenger" fiktiv verlaengert
  //   ka = 0, sonst
  //   test      Vektor von Flags, die angeben, ob Macro-Element schon getestet (=1)
  //   ele       Integer-Vektor der Dimension Anzahl der Macro-Elemente
  //             zur Speicherung der neuen Verfeinerungskanten
  //             (wird nur benoetigt, wenn umbvk=umb_vkant_macrodat) 

  void MacroReader::recumb(Mesh *mesh,  
			   MacroElement *mel, MacroElement *macroalt,
			   int *test, double lg, int ka, int *ele, 
			   void (*umbvk)(Mesh *mesh, MacroElement*, int k, int *el))
  {
    double l[3];
    int v[3];
    int k = 0;
    MacroElement *n;

    int vertices = mesh->getGeo(VERTEX);

    if (!mel || test[mel->getIndex()]==1) {
      return;
    } else {
      test[mel->getIndex()]=1; 
      
      laengstekante(mel->coord,l,v);
      
      if (v[1] == mesh->getGeo(VERTEX)) {            /*nur eine laengste Kante*/
	umbvk(mesh, mel, v[0], ele);      
	  
	for (int i = 0; i < vertices; i++) {
	  recumb(mesh,
		 mel->neighbour[i],
		 mel,
		 test,
		 l[0],
		 0,
		 ele,
		 umbvk);
	}
	return;
      } else { 
	if (ka == 1) {
	  if (fabs((l[0]-lg)/lg) < 0.0000001) {
	    for (int i = 0; i < vertices; i++) 
	      if (mel->neighbour[i] == macroalt)
		k = i;	    
		  
	    umbvk(mesh,mel,k,ele);
		  
	    for (int i = 0; i < vertices; i++) {
	      recumb(mesh,
		     mel->neighbour[i],
		     mel,
		     test,
		     l[0],
		     0,
		     ele,
		     umbvk);
	    }
	    return;
	  }
	}
	  
	n = const_cast<MacroElement*>(mel->getNeighbour(v[0])); 
	/*Nachbar an fiktiv verlaengerter Kante*/
	umbvk(mesh, mel, v[0], ele);   
	  
	recumb(mesh, n, mel,test,l[0],1,ele,umbvk);
	for (int i = 0; i < vertices; i++) { 
	  recumb(mesh,
		 mel->neighbour[i],
		 mel,
		 test,
		 l[0],
		 0,
		 ele,
		 umbvk);
	}
	return;
      }
    }             
  }
  
  //   berechnet aus Koord. der Eckpunkte eines Elementes 
  //   die Laenge der laengsten Kante

  //   l[0] = Laenge der laengsten Kante
  //   v[0] = lokaler Index der laengsten Kante
  //   v[1] = anderer lokaler Index, wenn ex. 2 laengste Kanten
  //        = 3, wenn ex. nur eine laengste Kante
  //   v[2] = dritter lokaler Index, wenn ex. 3 laengste Kanten
  //        = 3, wenn ex. nur eine oder zwei laengste Kanten      

  void MacroReader::laengstekante(FixVec<WorldVector<double>,VERTEX> coord, double *l, int *v)
  {
    int dim = coord.getSize() - 1;

    int i;
    int k;
    double lg;
    double eps;
    double kz; 

    int vertices = Global::getGeo(VERTEX,dim);

    l[0]=absteukl(coord[1],coord[2]);
    l[1]=absteukl(coord[0],coord[2]);
    l[2]=absteukl(coord[0],coord[1]);

    lg=l[0];
    kz=l[0];
    for (i=1; i < vertices; i++)
      {
	if (l[i] > lg) lg=l[i];
	if (l[i] < kz) kz=l[i];
      }

    eps=std::min(0.000001,kz/10000);
    k=0;
    for (i=0; i < vertices; i++)
      { 
	if (fabs(l[i]-lg) < eps)
	  {
	    v[k]=i;
	    k++;
	  } 
      }
    for (i=k; i < vertices; i++)
      {
	v[i]=Global::getGeo(VERTEX,dim);
      }

    l[0]=lg;
  }

  void MacroReader::checkMesh(Mesh *mesh)
  {
    FUNCNAME("MacroReader::checkMesh()");
    
    int nused, nfree;
    DOFAdmin  *localAdmin=mesh->admin[0];
    Flag fill_flag;
    int error_detected = 0;

    MSG("Checking mesh ...\n");

    fill_flag = Mesh::CALL_EVERY_EL_PREORDER | Mesh::FILL_NEIGH | Mesh::FILL_BOUND;

    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1, fill_flag | Mesh::FILL_ADD_ALL);
    while (elInfo) {
      basicCheckFct(elInfo, mesh);
      elInfo = stack.traverseNext(elInfo);
    }

    if (mesh->preserveCoarseDOFs)
      fill_flag = Mesh::CALL_EVERY_EL_PREORDER;
    else 
      fill_flag = Mesh::CALL_LEAF_EL;
  
    fill_flag |= Mesh::FILL_NEIGH | Mesh::FILL_BOUND;

    for (unsigned int iadmin = 0; iadmin < mesh->admin.size(); iadmin++) {
      localAdmin = mesh->admin[iadmin];
   
      if (localAdmin->getSize() > 0) {
	if (static_cast<int>(mesh->dof_used.size()) < localAdmin->getSize()) 
	  mesh->dof_used.resize(localAdmin->getSize() + 1000);
	
	for (unsigned int i = 0; i < mesh->dof_used.size(); i++) 
	  mesh->dof_used[i] = 0;
	
	nused = nfree = 0;

	TraverseStack stack;
	ElInfo *elInfo = 
	       stack.traverseFirst(mesh, -1, fill_flag | Mesh::FILL_ADD_ALL);
	while (elInfo) {
	  basicDOFCheckFct(elInfo, mesh, iadmin);
	  elInfo = stack.traverseNext(elInfo);
	}
	
	DOFIteratorBase it(localAdmin, USED_DOFS);
	for (it.reset(); !it.end(); ++it) {
	  nused++;
	  if (!mesh->dof_used[it.getDOFIndex()])  {
	    error_detected++;
	    MSG("dof[%d] not used??\n",it.getDOFIndex());
	  }
	}
	
	DOFIteratorBase freeIt(localAdmin, FREE_DOFS);
	for (freeIt.reset(); !freeIt.end(); ++freeIt) {
	  nfree++;
	  if (mesh->dof_used[freeIt.getDOFIndex()]) {
	    error_detected++;
	    MSG("dof[%d] used??\n",freeIt.getDOFIndex());
	  }
	}
	
	TEST(nused + nfree == localAdmin->getSize())
	  ("nused = %d, nfree = %d, admin.size = %d ????\n",
	   nused, nfree, localAdmin->getSize());
	TEST(nused == localAdmin->getUsedDofs())
	  ("nused = %d, admin.used_count = %d ?????\n",
	   nused, localAdmin->getUsedDofs());
      }
    }

    if (!error_detected) {
      MSG("checking done; no error detected\n");
    } else {
      MSG("checking done; %d error%s detected\n", error_detected,
	  error_detected == 1 ? "" : "s");

      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_INORDER);
      while (elInfo) {
	basicNodeFct(elInfo, mesh);
	elInfo = stack.traverseNext(elInfo);
      }
      WAIT_REALLY;
    }
  }


  int MacroReader::basicCheckFct(ElInfo* elinfo, Mesh *mesh)
  {
    FUNCNAME("MacroReader::basicCheckFct()");

    int j, k, opp_v;
    Element *el = elinfo->getElement();
    int error_detected = 0;
    int dim = mesh->getDim();

    elinfo->testFlag(Mesh::FILL_NEIGH);

    for (int i = 0; i < mesh->getGeo(NEIGH); i++) {
      const Element *neig = elinfo->getNeighbour(i);

      if (neig) {
	// If element has neighbour but also a non periodic boundary, this is
	// an error.
	if (elinfo->getBoundary(i) > 0) { // < 0 => periodic boundary
	  if (!error_detected)
	    MSG("error detected!!!\n");
	  error_detected++;
	  MSG("interior (*boundary)[%d] non NULL on element = %d\n", 
	      i, el->getIndex());
	}
	
	opp_v = elinfo->getOppVertex(i);
	if (opp_v < 0 || opp_v >= mesh->getGeo(NEIGH)) {
	  if (!error_detected)
	    MSG("error detected!!!\n");
	  error_detected++;
	  MSG("opp_v = %d\n", opp_v);
	}
	
	if (elinfo->getBoundary(i) > 0) {  // < 0 => periodic boundary
	  if (dim == 1) {
	    if (el->getDof(i) != neig->getDof(opp_v)) {
	      if (!error_detected)
		MSG("error detected!!!\n");
	      error_detected++;
	      MSG("neighbour error\n");
	    }
	  } else {
	    for (j = 1; j < mesh->getGeo(VERTEX); j++) {
	      for (k = 1; k < mesh->getGeo(VERTEX); k++)
		if (el->getDof((i+j) % mesh->getGeo(VERTEX)) == 
		    neig->getDof((opp_v+k) % mesh->getGeo(VERTEX))) 
		  break;
	      
	      if (k >= mesh->getGeo(VERTEX)) {
		if (!error_detected)
		  MSG("error detected!!!\n");
		error_detected++;
		MSG("dof %d of el %d at face %d isn't dof of neigh %d at face %d\n",
		    el->getDof((i+j) % 3,0), el->getIndex(), i, neig->getIndex(), 
		    opp_v);
	      }	   
	    }
	  }
	}
      } else {
	if (elinfo->getBoundary(i) == INTERIOR) {
	  if (!error_detected)
	    MSG("error detected!!!\n");
	  error_detected++;
	  MSG("(*boundary)[%d] on domains boundary is NULL on element = %d\n",
	      i, el->getIndex());
	}
      }
    }

    return error_detected;
  }

  void MacroReader::basicDOFCheckFct(ElInfo* elinfo, Mesh *mesh, int iadmin)
  {
    FUNCNAME("MacroReader::basicDOFCheckFct()");
    
    Element* el = elinfo->getElement();
    const DOFAdmin& admin = mesh->getDofAdmin(iadmin);
    const Element *neig;
    const DegreeOfFreedom *dof;

    if (mesh->dof_used.size() == 0)
      return;

    int ndof = admin.getNumberOfDofs(VERTEX);
    if (ndof) {
      int j0 = admin.getNumberOfPreDofs(VERTEX);
      TEST_EXIT(j0 + ndof <= mesh->getNumberOfDofs(VERTEX))
	("admin.getNumberOfPreDofs(VERTEX) %d + nDOF %d > mesh->nDOF %d\n",
	 j0, ndof, mesh->getNumberOfDofs(VERTEX));
      int i0 = mesh->getNode(VERTEX);
      for (int i = 0; i < mesh->getGeo(VERTEX); i++) {
	if ((dof = el->getDof(i0 + i)) == NULL) {
	  ERROR("no vertex dof %d on element %d\n", i, el->getIndex());
	} else {
	  for (int j = 0; j < ndof; j++) {
	    int jdof = dof[j0 + j];
	    TEST(jdof >= 0 && jdof < static_cast<int>(mesh->dof_used.size()))
	      ("vertex dof = %d invalid? size = %d\n", jdof, mesh->dof_used.size());
	    mesh->dof_used[jdof]++;
	  }
	}
      }
      /* neighbour vertex dofs have been checked in check_fct() */
    }

  
    if (mesh->getDim() > 1) {
      ndof = admin.getNumberOfDofs(EDGE);
      
      if (ndof) {
	// === Check for higher order DOFs on edges. ===

	int j0 = admin.getNumberOfPreDofs(EDGE);

	TEST_EXIT(j0 + ndof <= mesh->getNumberOfDofs(EDGE))
	  ("admin.getNumberOfPreDofs(EDGE) %d + nDOF %d > mesh->nDOF %d\n",
	   j0, ndof, mesh->getNumberOfDofs(EDGE));

	int i0 = mesh->getNode(EDGE);
	
	for (int i = 0; i < mesh->getGeo(EDGE); i++) {
	  dof = el->getDof(i0 + i);

	  if (dof == NULL) {
	    ERROR("no edge dof %d on element %d\n", i, el->getIndex());
	  } else {
	    for (int j = 0; j < ndof; j++) {
	      int jdof = dof[j0 + j];
	      TEST(jdof >= 0 && jdof < static_cast<int>(mesh->dof_used.size()))
		("edge dof=%d invalid? size=%d\n",jdof, mesh->dof_used.size());
	      mesh->dof_used[jdof]++;
	    }
	  }

	  if (el->getFirstChild() == NULL) {
	    if (mesh->getDim() == 2) {
	      neig = elinfo->getNeighbour(i);

	      // Just check the edge if it is not part of a periodic boundary. In
	      // the case of periodic boundaries the DOFs on neighbouring edges
	      // cannot fit together!

	      if (neig && elinfo->getBoundary(i) >= 0) {
		int ov = elinfo->getOppVertex(i);

		TEST(neig->getDof(i0 + ov) == dof)
		  ("el %d edge %d dof %8X: wrong dof %8X in neighbour %d edge %d\n",
		   el->getIndex(), i, dof, neig->getDof(i0 + ov), 
		   neig->getIndex(), ov);
	      }
	    } else { // dim == 3
	      for (int in = 0; in < mesh->getGeo(NEIGH); in++) {
		if ((in != el->getVertexOfEdge(i,0)) && 
		    (in != el->getVertexOfEdge(i,1)) &&
		    (neig = elinfo->getNeighbour(in))) {
		  int found = 0;
		  for (int k = 0; k < mesh->getGeo(EDGE); k++) {
		    if (neig->getDof(i0 + k) == dof) found++;
		  }
		  TEST(found==1)("el %d edge %d dof found=%d in neighbour %d\n",
				 el->getIndex(), i, found, neig->getIndex());
		}
	      }
	    }
	  }
	}
      }
    }

    if (mesh->getDim() == 3) {
      ndof = admin.getNumberOfDofs(FACE);
      if (ndof) {
	int j0 = admin.getNumberOfPreDofs(FACE);
	TEST_EXIT(j0 + ndof <= mesh->getNumberOfDofs(FACE))
	  ("admin->n0_dof[FACE] %d + nDOF %d > mesh->nDOF %d\n",
	   j0, ndof, mesh->getNumberOfDofs(FACE));
	int i0 = mesh->getNode(FACE);
	for (int i = 0; i < mesh->getGeo(FACE); i++) {
	  TEST(dof = el->getDof(i0 + i))("no face dof %d ???\n", i);
	  for (int j = 0; j < ndof; j++) {
	    int jdof = dof[j0 + j];
	    TEST(jdof >= 0 && jdof < static_cast<int>(mesh->dof_used.size()))
	      ("face dof=%d invalid? size=%d\n",jdof, mesh->dof_used.size());
	    mesh->dof_used[jdof]++;
	  }

	  if (el->getChild(0) == NULL) {
	    if ((neig = elinfo->getNeighbour(i))) {
	      int ov = elinfo->getOppVertex(i);

	      TEST(neig->getDof(i0 + ov) == dof)
		("el %d face %d dof %8X: wrong dof %8X in neighbour %d face %d\n",
		 el->getIndex(), i, dof, neig->getDof(i0 + ov), neig->getIndex(),
		 ov);
	    }
	  }
	}
      }
    }

    ndof = admin.getNumberOfDofs(CENTER);
    if (ndof) {
      int i0 = mesh->getNode(CENTER);
      TEST(dof = el->getDof(i0))("no center dof???\n");
      int j0 = admin.getNumberOfPreDofs(CENTER);
      TEST_EXIT(j0 + ndof <= mesh->getNumberOfDofs(CENTER))
	("admin.getNumberOfPreDofs(CENTER) %d + nDOF %d > mesh->nDOF %d\n",
	 j0, ndof, mesh->getNumberOfDofs(CENTER));
      for (int j = 0; j < ndof; j++) {
	int jdof = dof[j0 + j];
	TEST(jdof >= 0 && jdof < static_cast<int>(mesh->dof_used.size()))
	  ("center dof=%d invalid? size=%d\n",jdof, mesh->dof_used.size());
	mesh->dof_used[jdof]++;
      }
    }
  }

  void MacroReader::basicNodeFct(ElInfo* elinfo, Mesh *mesh)
  {
    FUNCNAME("MacroReader::basicNodeFct()");

    Element *el = elinfo->getElement();
    MSG("el %4d: ", el->getIndex());
    for (int i = 0; i < mesh->getGeo(VERTEX); i++)
      Msg::print("%4d%s", el->getDof(i,0), 
		 i < mesh->getGeo(VERTEX)-1 ? ", " : "\n");
  }
  
} } // end namespace io, AMDiS
