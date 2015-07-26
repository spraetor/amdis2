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

#include <cstring>

#include "MacroInfo.h"
#include "Mesh.h"
#include "MacroReader.h"
#include "FixVecConvert.h"
#include "SurfaceRegion_ED.h"
#include "ElementRegion_ED.h"
#include "MacroElement.h"

namespace AMDiS 
{
  void MacroInfo::fill(Mesh *pmesh, int nelements, int nvertices)
  {
    FUNCNAME("MacroInfo::fill()");

    TEST_EXIT(pmesh)("no mesh\n");

    mesh = pmesh;
    nElements = nelements;
    nVertices = nvertices;
    
    mesh->setNumberOfElements(nElements);
    mesh->setNumberOfLeaves(nElements);
    mesh->setNumberOfVertices(nVertices);

    for (int i = 0; i < nElements; i++) {
      MacroElement *newMacro = new MacroElement(mesh->getDim());
      mel.push_back(newMacro);
      mesh->addMacroElement(mel[i]);
    }

    dof = new DegreeOfFreedom*[nVertices];
    coords = new WorldVector<double>[nVertices];
    mel_vertex = new int*[nElements];

    for (int i = 0; i < nElements; i++)
      mel_vertex[i] = new int[mesh->getGeo(VERTEX)];

    for (int i = 0; i < nVertices; i++)
      dof[i] = mesh->getDof(VERTEX);

    for (int i = 0; i < nElements; i++) {
      mel[i]->element = mesh->createNewElement();
      mel[i]->index = i;
      mel[i]->elType = 0;
    }
    neigh_set = false;
    neigh_inverse_set = false;
    bound_set = false;
  }


  void MacroInfo::clear()
  {
    for (int i = 0; i < mesh->getNumberOfMacros(); i++)
      delete [] mel_vertex[i];

    delete [] mel_vertex;
    delete [] coords;
    coords = NULL;  
    delete [] dof;
    dof = NULL;

    mesh = NULL;
    neigh_set = false;
    neigh_inverse_set = false;
    
    nElements = 0;
    nVertices = 0;
    initialized = false;
  }


  /********************************************************************************/
  /*  read_indices()  reads dim + 1 indices from  file  into  id[0 - dim],        */
  /*    returns true if dim + 1 inputs arguments could be read successfully by    */
  /*    fscanf(), else false                                                      */
  /********************************************************************************/

  int MacroInfo::read_indices(FILE *file, Vector<int> &id)
  {
    int dim = mesh->getDim();
    TEST_EXIT_DBG(id.getSize() == dim+1)("Vector has wrong dimension!");

    for (int i = 0; i <= dim; i++)
      if (fscanf(file, "%d", &id[i]) != 1)
	return false;

    return true;
  }

#define N_KEYS      15
#define N_MIN_KEYS  7
  static const char *keys[N_KEYS] = {
    "DIM",                   //  0 
    "DIM_OF_WORLD",          //  1
    "number of vertices",    //  2
    "number of elements",    //  3
    "vertex coordinates",    //  4
    "element vertices",      //  5
    "element boundaries",    //  6
    "element neighbours",    //  7
    "element type",          //  8
    "projections",           //  9
    "element region",        // 10
    "surface region",        // 11
    "mesh name",             // 12
    "time",                   // 13
    "element neighbours inverse"    //  14   (for graph-structured meshes with >1 neighbors)
  };


  static int get_key_no(const char *key)
  {
    for (int i = 0; i < N_KEYS; i++)
      if (!strcmp(keys[i], key))  
	return i;

    return -1;
  }

#include <ctype.h>

  static const char *read_key(const char *line)
  {
    static char key[100];
    char *k = key;

    while (isspace(*line)) 
      line++;
    while ((*k++ = *line++) != ':');
    *--k = '\0';
  
    return const_cast<const char *>(key);
  }

  void MacroInfo::readAMDiSMacro(std::string filename, Mesh* mesh)
  {
    FUNCNAME("MacroInfo::readAMDiSMacro()");
    
    if (initialized) {
      WARNING("The old data in macro info will be removed.\n");
      clear();
    }

    int dim, dow;
    int nElements, nVertices;
    int j, k;
    double dbl;
    char line[256];
    int line_no, n_keys, sort_key[N_KEYS], nv_key, ne_key;
    int key_def[N_KEYS] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
    const char *key;
    DimVec<int> *ind = NULL;
    FixVec<int, NEIGH> *ind_neigh = NULL;

    TEST_EXIT(filename != "")("No filename specified!\n");

    FILE *file = fopen(filename.c_str(), "r");
    TEST_EXIT(file)("cannot open file %s\n", filename.c_str());


    // === Looking for all keys in the macro file. ===

    line_no = n_keys = 0;
    while (fgets(line, 255, file)) {
      line_no++;
      if (!strchr(line, ':'))  
	continue;
      key = read_key(line);
      int i_key = get_key_no(key);
      TEST_EXIT(i_key >= 0)
	("macro file %s must not contain key %s on line %d\n",
	 filename.c_str(), key, line_no);
      TEST_EXIT(!key_def[i_key])
	("key %s defined second time on line %d in file %s\n");

      sort_key[n_keys++] = i_key;
      key_def[i_key] = true;
    }
    fclose(file);


    // === Test, if there is data for every key and if all is defined in === 
    // === right order.                                                  ===

    for (int i_key = 0; i_key < N_MIN_KEYS; i_key++) {
      for (j = 0; j < n_keys; j++)
	if (sort_key[j] == i_key)  
	  break;

      TEST_EXIT(j < n_keys)("You do not have specified data for %s in %s\n",
			    keys[i_key], filename.c_str());

      for (j = 0; j < n_keys; j++)
	if (sort_key[j] == 2)  break;
      nv_key = j;
      for (j = 0; j < n_keys; j++)
	if (sort_key[j] == 3)  break;
      ne_key = j;
    
      switch (i_key) {
      case 0:
      case 1:
	TEST_EXIT(sort_key[i_key] < 2)
	  ("You have to specify DIM or mesh->getGeo(WORLD) before all other data\n");
	break;
      case 4: 
	TEST_EXIT(nv_key < i_key)
	  ("Before reading data for %s, you have to specify the %s in file\n",
	   keys[4], keys[2], filename.c_str());
	break;
      case 5: 
	TEST_EXIT(nv_key < i_key  &&  ne_key < i_key)
	  ("Before reading data for %s, you have to specify the %s and %s in file %s\n",
	   keys[5], keys[3], keys[2], filename.c_str());
      case 6:
      case 7:
      case 8:
	TEST_EXIT(ne_key < i_key)
	  ("Before reading data for %s, you have to specify the %s in file %s\n",
	   keys[i_key], keys[3], filename.c_str());
      }
    }

    for (int i_key = 0; i_key < N_KEYS; i_key++)
      key_def[i_key] = false;


    // === And now, reading data. ===
	
    file = fopen(filename.c_str(), "r");
    TEST_EXIT(file)("cannot open file %s\n", filename.c_str());

    int result;

    for (int i_key = 0; i_key < n_keys; i_key++) {

      switch (sort_key[i_key]) {
	
      case 0:
	// line "DIM"
	result = fscanf(file, "%*s %d", &dim);
	TEST_EXIT(result == 1)("cannot read DIM correctly in file %s\n", filename.c_str());

	ind = new DimVec<int>(dim, NO_INIT);
	ind_neigh = new FixVec<int, NEIGH>(dim, NO_INIT);

	key_def[0] = true;
	break;

      case 1:
	// line "DIM_OF_WORLD"
	result = fscanf(file, "%*s %d", &dow);
	TEST_EXIT(result == 1)
	  ("cannot read Global::getGeo(WORLD) correctly in file %s\n", filename.c_str());
	TEST_EXIT(dow == Global::getGeo(WORLD))
	  ("dimension of world = %d != Global::getGeo(WORLD) = %d\n", 
	   dow, Global::getGeo(WORLD));

	key_def[1] = true;
	break;

      case 2:
	// line "number of vertices"
	result = fscanf(file, "%*s %*s %*s %d", &nVertices);
	TEST_EXIT(result == 1)
	  ("cannot read number of vertices correctly in file %s\n", filename.c_str());
	TEST_EXIT(nVertices > 0)
	  ("number of vertices = %d must be bigger than 0\n", nVertices);

	key_def[2] = true;
	if (key_def[3])
	  fill(mesh, nElements, nVertices);
	break;

      case 3:
	// line "number of elements"
	result = fscanf(file, "%*s %*s %*s %d", &nElements);
	TEST_EXIT(result == 1)
	  ("cannot read number of elements correctly in file %s\n", filename.c_str());
	TEST_EXIT(nElements > 0)
	  ("number of elements = %d must be bigger than 0\n", nElements);

	key_def[3] = true;
	if (key_def[2])
	  fill(mesh, nElements, nVertices);
	break;

      case 4:
	// block "vertex coordinates"
	result = fscanf(file, "%*s %*s");
	for (int i = 0; i < nVertices; i++) {
	  for (j = 0; j <Global::getGeo(WORLD) ; j++) {
	    result = fscanf(file, "%lf", &dbl);
	    TEST_EXIT(result == 1)
	      ("error while reading coordinates, check file %s\n", filename.c_str());
	    coords[i][j] = dbl;
	  }
	}
	key_def[4] = true;
	break;

      case 5:
	// block "element vertices"
	result = fscanf(file, "%*s %*s");

	// === Global number of vertices on a single element. ===

	for (int i = 0; i < nElements; i++) {
	  result = read_indices(file, *ind);
	  TEST_EXIT(result)
	    ("cannot read vertex indices of element %d in file %s\n",  i, filename.c_str());

	  for (k = 0; k < mesh->getGeo(VERTEX); k++)
	    mel_vertex[i][k] = (*ind)[k];
	}

	key_def[5] = true;
	break;

      case 6:
	// block "element boundaries"
	result = fscanf(file, "%*s %*s");


	// === MEL boundary pointers. ===

	for (int i = 0; i < nElements; i++) {  
	  // boundary information of ith element 

	  result = read_indices(file, *ind_neigh);
	  TEST_EXIT(result)
	    ("cannot read boundary type of element %d in file %s\n", i, filename.c_str());

	  // fill boundary of macro-element
	  io::MacroReader::fillMelBoundary(mesh, mel[i], *ind_neigh);
	}

	this->fillBoundaryInfo(mesh);
                   
	bound_set = true;
	key_def[6] = true;
	break;

      case 7:
	// block "element neighbours"
	result = fscanf(file, "%*s %*s");

	// ===  Fill MEL neighbour pointers:                               ===
	// ===    if they are specified in the file: read them from file,  ===
	// ===    else init them by a call of fill_mel_neighbour()         ===

	neigh_set = true;
	for (int i = 0; i < nElements; i++) {
	  //  neighbour information about ith element

	  if (read_indices(file, *ind_neigh)) {
	    io::MacroReader::fillMelNeigh(mel[i], mel, *ind_neigh);
	  } else {
	    // setting of neighbours fails

	    neigh_set = false; 
	    break;
	  }
	}

	key_def[7] = true;
	break;

      case 8:
	// block "element type"
	result = fscanf(file, "%*s %*s");

	// === MEL elType ===

	if (dim == 2 || dim == 1)
	  ERROR("there is no element type in 2d and 2d; ignoring data for elType\n");

	for (int i = 0; i < nElements; i++) {
	  result = fscanf(file, "%d", &j);
	  TEST_EXIT(result == 1)
	    ("cannot read elType of element %d in file %s\n", i, filename.c_str());
	  if (dim == 3)
	    (mel)[i]->elType = j;
	}

	key_def[8] = true;
	break;

      case 9:
	// block "projections"
	{
	  result = fscanf(file, "%*s");

	  int nFaces = mesh->getGeo(FACE);
	  int nEdgesAtBoundary = 0;

	  for (k = 1; k < dim; k++)
	    nEdgesAtBoundary += k;

	  for (int i = 0; i < nElements; i++) {
	    result = read_indices(file, *ind);
	    TEST_EXIT(result)
	      ("cannot read boundary projector of element %d in file %s\n", i, filename.c_str());
	
	    Projection *projector = Projection::getProjection((*ind)[0]);

	    if (projector && projector->getType() == VOLUME_PROJECTION) {
	      mel[i]->setProjection(0, projector);
	    } else { // boundary projection
	      for (j = 0; j < mesh->getGeo(NEIGH); j++) {
		projector = Projection::getProjection((*ind)[j]);
		if (projector) {
		  mel[i]->setProjection(j, projector);
		  if (dim > 2) {
		    for (k = 0; k < nEdgesAtBoundary; k++) {
		      int edgeNr = Global::getReferenceElement(dim)->getEdgeOfFace(j, k);
		      mel[i]->setProjection(nFaces + edgeNr, projector);
		    }
		  }
		}
	      }
	    }
	  }
	}
	key_def[9] = true;
	break;

      case 10:
	// block "element region"
	result = fscanf(file, "%*s %*s");

	// === MEL regions. ===

	for (int i = 0; i < nElements; i++) {
	  result = fscanf(file, "%d", &j);
	  TEST_EXIT(result == 1)
	    ("cannot read region of element %d in file %s\n", i, filename.c_str());
	  if (j >= 0) {
	    Element *el = mel[i]->getElement();
	    ElementRegion_ED *elementRegion = 
	      new ElementRegion_ED(el->getElementData());
	    elementRegion->setRegion(j);
	    el->setElementData(elementRegion);
	  }
	}
	key_def[10] = true;
	break;

      case 11:
	// block "surface region"
	result = fscanf(file, "%*s %*s");
	for (int i = 0; i < nElements; i++) {
	  result = read_indices(file, *ind);
	  TEST_EXIT(result)
	    ("cannot read surface regions of element %d in file %s\n", i, filename.c_str());

	  Element *el = mel[i]->getElement();

	  for (j = 0; j < mesh->getGeo(NEIGH); j++) {
	    if ((*ind)[j] >= 0) {
	      SurfaceRegion_ED *surfaceRegion = 
		new SurfaceRegion_ED(el->getElementData());
	      surfaceRegion->setSide(j);
	      surfaceRegion->setRegion((*ind)[j]);
	      el->setElementData(surfaceRegion);
	    }
	  }
	}
	key_def[11] = true;
	break;

      case 12:
	// line "mesh name"
	result = fscanf(file, "%*s %*s %*s");
	break;

      case 13:
	// line "time"
	result = fscanf(file, "%*s %*s");
	break;

      case 14:
	// block "element neighbours inverse"
	result = fscanf(file, "%*s %*s %*s");

	// ===  Fill MEL neighbour pointers:                               ===
	// ===    if they are specified in the file: read them from file,  ===
	// ===    else init them by a call of fill_mel_neighbour()         ===

	neigh_inverse_set = true;
	for (int i = 0; i < nElements; i++) {
	  //  neighbour information about ith element

	  if (read_indices(file, *ind_neigh)) {
	    io::MacroReader::fillMelNeighInv(mel[i], mel, *ind_neigh);
	  } else {
	    // setting of neighbours fails
	    ERROR_EXIT("Problem while reading 'element neighbours inverse' of element %d\n", i);
	    neigh_inverse_set = false; 
	    break;
	  }
	}

	key_def[7] = true;
	break;

      }
    }

    if (ind)
      delete ind;

    fclose(file);
    
    initialized = true;
  }


  void MacroInfo::dirichletBoundary()
  {
    // === Traverse all elements and set either interior boundaries, if the    ===
    // === element has a neighbour on this bound. Otherwise set some arbitrary ===
    // === Dirichlet boundary conditions.                                      ===

    for (int i = 0; i < static_cast<int>(mel.size()); i++)
      for (int k = 0; k < mesh->getGeo(NEIGH); k++)
	if (mel[i]->neighbour[k])
	  mel[i]->boundary[k] = INTERIOR;
	else if (neigh_inverse_set && mel[i]->neighbour_inv[k])
	  mel[i]->boundary[k] = INTERIOR;
	else
	  mel[i]->boundary[k] = 1;
  }


  void MacroInfo::fillBoundaryInfo(Mesh *mesh)
  {
    FUNCNAME("MacroInfo::fillBoundaryInfo()");

    int i;
    int nVertices = mesh->getNumberOfVertices();
    std::deque<MacroElement*>::iterator melIt;
    std::vector<BoundaryType> bound(nVertices);

    switch (mesh->getDim()) {
    case 1:
      break;
    case 2:
      for (i = 0; i < nVertices; i++)
	bound[i] = INTERIOR;

      for (i = 0, melIt = mesh->firstMacroElement(); 
	   melIt != mesh->endOfMacroElements(); 
	   ++melIt, i++) {
	for (int j = 0; j < mesh->getGeo(NEIGH); j++) {
	  if ((*melIt)->getBoundary(j) != INTERIOR) {
	    if ((*melIt)->getBoundary(j) > 0) {
	      // Here we have some Dirichlet boundary conditions.

	      int j1 = mel_vertex[i][(j + 1) % 3];
	      int j2 = mel_vertex[i][(j + 2) % 3];
	      
	      bound[j1] = std::max(bound[j1], (*melIt)->getBoundary(j));
	      bound[j2] = std::max(bound[j2], (*melIt)->getBoundary(j));
	    } else {
	      // Otherwise we have some Neumann boundary conditions.

	      int j1 = mel_vertex[i][(j + 1) % 3];
	      int j2 = mel_vertex[i][(j + 2) % 3];
	      
	      if (bound[j1] != INTERIOR)
		bound[j1] = std::max(bound[j1], (*melIt)->getBoundary(j));
	      else
		bound[j1] = (*melIt)->getBoundary(j);
	      
	      if (bound[j2] != INTERIOR)
		bound[j2] = std::max(bound[j2], (*melIt)->getBoundary(j));
	      else
		bound[j2] = (*melIt)->getBoundary(j);
	    }
	  }
	}
      }

      for (i = 0, melIt = mesh->firstMacroElement(); 
	   melIt != mesh->endOfMacroElements(); 
	   ++melIt, i++) 
	for (int j = 0; j < mesh->getGeo(VERTEX); j++)
	  (*melIt)->setBoundary(3 + j, bound[mel_vertex[i][j]]);
	
      break;
    case 3:
      for (i = 0; i < nVertices; i++)
	bound[i] = INTERIOR;

      for (i = 0, melIt = mesh->firstMacroElement(); 
	   melIt != mesh->endOfMacroElements(); 
	   ++melIt, i++) {
	for (int j = 0; j < mesh->getGeo(NEIGH); j++) {
	  for (int k = 1; k < 4; k++)
	    bound[mel_vertex[i][(j + k) % 4]] =
	      ((*melIt)->getBoundary(j) != INTERIOR) ?
	      newBound((*melIt)->getBoundary(j),
		       bound[mel_vertex[i][(j + k) % 4]]) :
	      bound[mel_vertex[i][(j + k) % 4]];
	}
      }

      for (i = 0, melIt = mesh->firstMacroElement(); 
	   melIt != mesh->endOfMacroElements(); 
	   ++melIt, i++) 	
	for (int j = 0; j < mesh->getGeo(VERTEX); j++)
	  (*melIt)->setBoundary(10 + j, bound[mel_vertex[i][j]]);
      
      break;
    default: 
      ERROR_EXIT("invalid dim\n");
    }
  }
  
}
