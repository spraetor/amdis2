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


#include "parallel/CheckerPartitioner.h"
#include "Traverse.h"

using namespace std;

namespace AMDiS { namespace Parallel {

  CheckerPartitioner::CheckerPartitioner(string name, MPI::Intracomm *comm)
    : MeshPartitioner(name, comm),
      mpiRank(mpiComm->Get_rank()),
      mpiSize(mpiComm->Get_size()),
      mode(0),
      multilevel(false)
  {
    string modestr = "";
    Parameters::get(initFileStr + "->mode", modestr);
    
    if (modestr == "x-stripes")
      mode = 1;
    else if (modestr == "y-stripes")
      mode = 2;
    else if (modestr == "z-stripes")
      mode = 3;
    else if (modestr == "tetrahedron-stripes")
      mode = 4;
    else if (modestr == "multilevel")
      multilevel = true;
    else {
      if (modestr != "") {
	ERROR_EXIT("No partitioner mode \"%s\"!\n", modestr.c_str());
      }
    }
  }
  

  bool CheckerPartitioner::createInitialPartitioning()
  {
    FUNCNAME("CheckerPartitioner::createInitialPartitioning()");

    // In one of the stripes mode, we have to check if the number of macro
    // elements with together with the number of nodes.
    if (mode == 1 || mode == 2 || mode == 3) {
      int elCounter = 0;
      TraverseStack stack;
      ElInfo *elInfo = stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL);
      while (elInfo) {
	elCounter++;
	elInfo = stack.traverseNext(elInfo);
      }
      
      if (mesh->getDim() == 2) {
	TEST_EXIT(elCounter == 2 * mpiSize * mpiSize)
	  ("The number of macro elements is %d, but must be %d for %d number of nodes!",
	   elCounter, 2 * mpiSize * mpiSize, mpiSize);
      }

      if (mesh->getDim() == 3) {
	TEST_EXIT(elCounter == 6 * static_cast<int>(pow(mpiSize, 1.5)))
	  ("The number of macro elements is %d, but must be %d for %d number of nodes!",
	   elCounter, 6 * static_cast<int>(pow(mpiSize, 1.5)), mpiSize);
      }

    }

    if (multilevel) {
      TEST_EXIT(MPI::COMM_WORLD.Get_size() == 16)
	("Multilevel partitioning is implemented for 16 nodes only!\n");
    }


    if (mode == 4) {
      TEST_EXIT(mesh->getDim() == 3)("Works only in 3D!\n");
      createTetrahedronStripes();
    }

    
    int dim = mesh->getDim();
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL);
    while (elInfo) {
      Element *el = elInfo->getElement();
      int elIndex = el->getIndex();
      int boxIndex = elIndex / (dim == 2 ? 2 : 6);
      int elInRank = -1;

      switch (mode) {
      case 0:
	if (!multilevel) {
	  elInRank = boxIndex;	
	} else {
	  TEST_EXIT_DBG(boxIndex >= 0 && boxIndex <= 15)("Wrong box index!\n");
	  int rs[16] = {0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};
	  elInRank = rs[boxIndex];
	}
	break;

      case 1:
	// x-slices

	{
	  if (dim == 2)
	    elInRank = elIndex / (2 * mpiSize);
	  
	  if (dim == 3) {
	    int boxSliceY = 
	      (boxIndex % mpiSize) / static_cast<int>(sqrt(mpiSize));
	    int boxSliceZ = boxIndex / mpiSize;
	    elInRank = boxSliceY * static_cast<int>(sqrt(mpiSize)) + boxSliceZ;
	  }
	}

	break;

      case 2:
	// y-slices

	{
	  if (dim == 2)
	    elInRank = (elIndex % (2 * mpiSize)) / 2;
	  
	  if (dim == 3) {
	    int boxSliceX = 
	      (boxIndex % mpiSize) % static_cast<int>(sqrt(mpiSize));
	    int boxSliceZ = boxIndex / mpiSize;
	    elInRank = boxSliceX * static_cast<int>(sqrt(mpiSize)) + boxSliceZ;
	  }
	}

	break;

      case 3:
	// z-slices

	{
	  int boxSliceX = (boxIndex % mpiSize) % static_cast<int>(sqrt(mpiSize));
	  int boxSliceY = (boxIndex % mpiSize) / static_cast<int>(sqrt(mpiSize));

	  elInRank = boxSliceX * static_cast<int>(sqrt(mpiSize)) + boxSliceY;
	}

	break;

      case 4:
	// tetrahedron-stripes

	elInRank = elStripeInRank[elIndex];

#if 0
	{

	  int nStripes = elInStripe.size();
	  int elsPerStripe = elInStripe[0].size();  
	  int procPerStripe = mpiSize / nStripes;
	  TEST_EXIT(mpiSize % nStripes == 0)("Should not happen!\n");
	  
	  int inStripe = -1;
	  int stripePos = -1;
	  for (int stripe = 0; stripe < nStripes; stripe++) {
	    for (int pos = 0; pos < elInStripe[stripe].size(); pos++) {
	      if (elInStripe[stripe][pos] == elIndex) {
		inStripe = stripe;
		stripePos = pos;
		break;
	      }
	    }
	    if (inStripe >= 0)
	      break;
	  }

	  TEST_EXIT(inStripe >= 0)("Should not happen!\n");

	  elInRank = inStripe;
	}
#endif

	break;

      default:
	ERROR_EXIT("Mode %d does not exists for checker based mesh partitioning!\n", 
		   mode);
      }

      TEST_EXIT_DBG(elInRank >= 0)("Should not happen!\n");
      TEST_EXIT_DBG(elInRank < mpiSize)("Should not happen!\n");
      
      elementInRank[elIndex] = (elInRank == mpiRank);
      partitionMap[elIndex] = elInRank;	
      
      elInfo = stack.traverseNext(elInfo);
    }

    return true;
  }


  void CheckerPartitioner::createTetrahedronStripes()
  {
    FUNCNAME("CheckerPartitioner::createTetrahedronStripes()");

    vector<vector<MacroElement*> > stripes;
    vector<vector<int> > elInStripe;

    int nElements = 0;
    TraverseStack stack;
    ElInfo *elInfo = 
      stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL | Mesh::FILL_COORDS);
    while (elInfo) {
      TEST_EXIT(elInfo->getLevel() == 0)("Should not happen!\n");

//       Element *el = elInfo->getElement();
//       int elIndex = el->getIndex();

      int zeroCoordCounter = 0;
      for (int i = 0; i < mesh->getGeo(VERTEX); i++)
	if (fabs(elInfo->getCoord(i)[2]) < 1e-10)
	  zeroCoordCounter++;      

      if (zeroCoordCounter == 3) {
	vector<MacroElement*> tmp;
	tmp.push_back(elInfo->getMacroElement());
	stripes.push_back(tmp);

	vector<int> tmpIndex;
	tmpIndex.push_back(elInfo->getMacroElement()->getIndex());
	elInStripe.push_back(tmpIndex);
      }

      nElements++;
      elInfo = stack.traverseNext(elInfo);
    }
    
    TEST_EXIT(mpiSize % stripes.size() == 0)
      ("Should not happen! mpiSize = %d but %d bottom elements found!\n",
       mpiSize, stripes.size());

    int testElementCounter = 0;
    for (size_t stripe = 0; stripe < stripes.size(); stripe++) {
      MacroElement *mel = stripes[stripe][0];

      set<int> localDofs;
      for (int i = 0; i < mesh->getGeo(VERTEX); i++)
	if (fabs(mel->getCoord(i)[2]) < 1e-10)
	  localDofs.insert(i);
      TEST_EXIT(localDofs.size() == 3)("Should not happen!\n");


      while (mel != NULL) {
	int replaceDof = -1;
	for (int i = 0; i < mesh->getGeo(VERTEX); i++)
	  if (localDofs.count(i) == 0)
	    replaceDof = i;
	
	bool found = false;
	for (std::set<int>::iterator dit = localDofs.begin(); 
	     dit != localDofs.end(); ++dit) {
	  WorldVector<double> c0 = mel->getCoord(*dit);
	  WorldVector<double> c1 = mel->getCoord(replaceDof);
	  
	  if (fabs(c0[0] - c1[0]) < 1e-10 &&
	      fabs(c0[1] - c1[1]) < 1e-10 &&
	      fabs(c0[2] - c1[2]) > 1e-10) {
	    found = true;
	    localDofs.erase(dit);
	    localDofs.insert(replaceDof);
	    break;
	  }
	}
	TEST_EXIT(found)("Should not happen!\n");
	
	
	set<DegreeOfFreedom> faceDofs;      
	for (std::set<int>::iterator dit = localDofs.begin(); 
	     dit != localDofs.end(); ++dit)
	  faceDofs.insert(mel->getElement()->getDof(*dit, 0));	
	TEST_EXIT(faceDofs.size() == 3)("Should not happen!\n");

	
	int localFace = -1;
	for (int i = 0; i < mesh->getGeo(FACE); i++) {
	  DofFace face = mel->getElement()->getFace(i);
	  bool allVertexInFace = 
	    faceDofs.count(std::get<0>(face)) &&
	    faceDofs.count(std::get<1>(face)) &&
	    faceDofs.count(std::get<2>(face));
	  if (allVertexInFace) {
	    localFace = i;
	    break;
	  }
	}
	TEST_EXIT(localFace >= 0)("Should not happen!\n");

	MacroElement *mel_neigh = mel->getNeighbour(localFace);
	if (mel_neigh) {
	  stripes[stripe].push_back(mel_neigh);
	  elInStripe[stripe].push_back(mel_neigh->getIndex());
	}	
	
	mel = mel_neigh;
      }

      testElementCounter += stripes[stripe].size();
    }

    TEST_EXIT(testElementCounter == nElements)("Should not happen!\n");

    size_t elsPerStripe = stripes[0].size();
    for (size_t i = 0; i < stripes.size(); i++) {
      TEST_EXIT(stripes[i].size() == elsPerStripe)
	("Should not happen!\n");
    }


    // === Computing mapping from macro element indices to ranks ===

    size_t nStripes = elInStripe.size();
    size_t procPerStripe = mpiSize / nStripes;
    size_t elsPerRank = elsPerStripe / procPerStripe;
    TEST_EXIT(mpiSize % nStripes == 0)("Should not happen!\n");
    TEST_EXIT(elsPerStripe % procPerStripe == 0)("Should not happen!\n");

    elStripeInRank.clear();

    int rankCount = 0;
    for (size_t i = 0; i < nStripes; i++) {
      for (size_t j = 0; j < procPerStripe; j++) {
	for (size_t k = 0; k < elsPerRank; k++)
	  elStripeInRank[elInStripe[i][j * elsPerRank + k]] = rankCount;
	rankCount++;
      }
    }
  }

} }
