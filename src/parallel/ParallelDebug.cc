#include "parallel/ParallelDebug.h"
#include "parallel/MeshDistributor.h"
#include "parallel/MpiHelper.h"
#include "ProblemStat.h"
#include "DOFVector.h"
#include "FixVec.h"
#include "StdMpi.h"
#include "Debug.h"
#include "Timer.h"
#include "io/VtkWriter.h"
#include "ElementDofIterator.h"

namespace AMDiS 
{ 
  namespace Parallel 
  {

  using namespace std;


  void ParallelDebug::testInteriorBoundary(MeshDistributor &pdb)
  {
    FUNCNAME("ParallelDebug::testInteriorBoundary()");

    vector<int*> sendBuffers, recvBuffers;

    MSG("WARNING: TEST ONLY FOR LEVEL 0!\n");
    MPI::Intracomm &mpiComm = pdb.levelData.getMpiComm(0);
    int mpiRank = mpiComm.Get_rank();

    MPI::Request *request = new MPI::Request[pdb.intBoundary[0].own.size() + 
			 pdb.intBoundary[0].other.size() +
                         pdb.intBoundary[0].periodic.size() * 2];
    int requestCounter = 0;


    // === Send rank's boundary information. ===

    for (RankToBoundMap::iterator rankIt = pdb.intBoundary[0].own.begin();
	 rankIt != pdb.intBoundary[0].own.end(); ++rankIt) {

      int nSendInt = rankIt->second.size();
      int* buffer = new int[nSendInt];
      for (int i = 0; i < nSendInt; i++)
	buffer[i] = (rankIt->second)[i].rankObj.elIndex;
      
      sendBuffers.push_back(buffer);
      
      request[requestCounter++] =
	mpiComm.Isend(buffer, nSendInt, MPI_INT, rankIt->first, 0);
    }


    // === Receive information from other ranks about the interior boundaries. ====

    for (RankToBoundMap::iterator rankIt = pdb.intBoundary[0].other.begin();
	 rankIt != pdb.intBoundary[0].other.end(); ++rankIt) {
      int nRecvInt = rankIt->second.size();
      int *buffer = new int[nRecvInt];
      recvBuffers.push_back(buffer);

      request[requestCounter++] = 
	mpiComm.Irecv(buffer, nRecvInt, MPI_INT, rankIt->first, 0);
    }


    // === To the last, do the same of periodic boundaries. ===

    for (RankToBoundMap::iterator rankIt = pdb.intBoundary[0].periodic.begin();       
	 rankIt != pdb.intBoundary[0].periodic.end(); ++rankIt) {
      if (rankIt->first == mpiRank)
	continue;

      int nValues = rankIt->second.size();
      int* sBuffer = new int[nValues];
      for (int i = 0; i < nValues; i++)
	sBuffer[i] = (rankIt->second)[i].rankObj.elIndex;

      sendBuffers.push_back(sBuffer);

      request[requestCounter++] =
	mpiComm.Isend(sBuffer, nValues, MPI_INT, rankIt->first, 0);

      int *rBuffer = new int[nValues];
      recvBuffers.push_back(rBuffer);

      request[requestCounter++] = 
	mpiComm.Irecv(rBuffer, nValues, MPI_INT, rankIt->first, 0);
    }

    // === Finish communication and delete all send buffers. ===

    MPI::Request::Waitall(requestCounter, request);
    for (int i = 0; i < static_cast<int>(sendBuffers.size()); i++)
      delete [] sendBuffers[i];


    // === Finally, check the results, i.e., the indices of element at the     === 
    // === boundaries, if they fit together. First check the interior bounds,  ===
    // === and after this the periodic ones.                                   ===

    int bufCounter = 0;
    for (RankToBoundMap::iterator rankIt = pdb.intBoundary[0].other.begin();
	 rankIt != pdb.intBoundary[0].other.end(); ++rankIt) {

      TEST_EXIT(rankIt->second.size() == 
		pdb.intBoundary[0].other[rankIt->first].size())
	("Boundaries does not fit together!\n");      

      for (unsigned int i = 0; i < rankIt->second.size(); i++) {
	int elIndex1 = recvBuffers[bufCounter][i];
	int elIndex2 = pdb.intBoundary[0].other[rankIt->first][i].neighObj.elIndex;

	TEST_EXIT(elIndex1 == elIndex2)("Wrong element index at interior boundary!\n");
      }

      delete [] recvBuffers[bufCounter++];
    }


    for (RankToBoundMap::iterator rankIt = pdb.intBoundary[0].periodic.begin();
	 rankIt != pdb.intBoundary[0].periodic.end(); ++rankIt) {
      if (rankIt->first == mpiRank)
	continue;

      for (unsigned int i = 0; i < rankIt->second.size(); i++) {
	int elIndex1 = recvBuffers[bufCounter][i];
	int elIndex2 = pdb.intBoundary[0].periodic[rankIt->first][i].neighObj.elIndex;

	TEST_EXIT(elIndex1 == elIndex2)
	  ("Wrong element index at periodic boundary el %d with rank %d: %d %d\n", 
	   pdb.intBoundary[0].periodic[rankIt->first][i].rankObj.elIndex,
	   rankIt->first, elIndex1, elIndex2);
      }

      delete [] recvBuffers[bufCounter++];
    }
    
    delete[] request;
  }


  void ParallelDebug::testPeriodicBoundary(MeshDistributor &pdb)
  {
    FUNCNAME("ParallelDebug::testPeriodicBoundary()");

    MSG("WARNING: TEST FOR PERIODIC BOUNDARIES IS DISABLED!\n");
//     for (unsigned int i = 0; i < pdb.feSpaces.size(); i++)
//       testPeriodicBoundary(pdb, pdb.feSpaces[i]);
  }


  void ParallelDebug::testPeriodicBoundary(MeshDistributor &pdb,
					   const FiniteElemSpace *feSpace)
  {
    FUNCNAME("ParallelDebug::testPeriodicBoundary()");

    // === 1. check: All periodic DOFs should have at least a correct number ===
    // === of periodic associations.                                         ===

    PeriodicMap &perMap = pdb.getPeriodicMap();
    for (map<int, std::set<BoundaryType> >::iterator it = 
	   perMap.periodicDofAssociations[feSpace].begin();
	 it != perMap.periodicDofAssociations[feSpace].end(); ++it) {
      WorldVector<double> c;
      pdb.macroMesh->getDofIndexCoords(it->first, pdb.feSpaces[0], c);
    }    


    // === 2. check: All periodic DOFs must be symmetric, i.e., if A is mapped ===
    // === to B, then B must be mapped to A.                                   ===

    MPI::Intracomm &mpiComm = pdb.levelData.getMpiComm(0);
    int mpiSize = mpiComm.Get_size();
    int mpiRank = mpiComm.Get_rank();

    StdMpi<PeriodicDofMap> stdMpi(mpiComm, true);

    if (mpiRank == 0) {
      for (int i = 1; i < mpiSize; i++)
	stdMpi.recv(i);
    } else {
      stdMpi.send(0, perMap.periodicDofMap[feSpace]);
    }

    stdMpi.startCommunication();

    int foundError = 0;

    // === The boundary DOFs are checked only on the zero rank. === 

    if (mpiRank == 0) {
      // Stores to each rank the periodic DOF mappings of this rank.
      map<int, PeriodicDofMap> rankToMaps;
      PeriodicDofMap dofMap = perMap.periodicDofMap[feSpace];
      rankToMaps[0] = dofMap;

      for (int i = 1; i < mpiSize; i++) {
	PeriodicDofMap &otherMap = stdMpi.getRecvData(i);
	rankToMaps[i] = otherMap;
	
	for (PeriodicDofMap::iterator it = otherMap.begin(); 
	     it != otherMap.end(); ++it) {
	  for (std::map<DegreeOfFreedom, DegreeOfFreedom>::iterator dofIt = it->second.begin();
	       dofIt != it->second.end(); ++dofIt) {
	    if (dofMap.count(it->first) == 1 &&
		dofMap[it->first].count(dofIt->first) == 1) {
	      TEST_EXIT_DBG(dofMap[it->first][dofIt->first] == dofIt->second)
		("Should not happen!\n");
	    } else {
	      dofMap[it->first][dofIt->first] = dofIt->second;
	    }
	  }
	}
      }


      // === Now we test if global DOF A is mapped to B, then B must be mapped ===
      // === to A for the same boundary type.                                  ===

      for (PeriodicDofMap::iterator it = dofMap.begin(); 
	   it != dofMap.end(); ++it) {
	for (std::map<DegreeOfFreedom, DegreeOfFreedom>::iterator dofIt = it->second.begin();
	     dofIt != it->second.end(); ++dofIt) {
	  if (it->second[dofIt->second] != dofIt->first) {
	    MSG("[DBG]  For boundary type %d: DOF %d -> %d, but %d -> %d!\n",
		it ->first, 
		dofIt->first, dofIt->second, 
		dofIt->second, it->second[dofIt->second]);

	    for (int i = 0; i < mpiSize; i++) {
	      if (rankToMaps[i][it->first].count(dofIt->first) == 1) {
		MSG("[DBG]    %d -> %d in rank %d\n", 
		    dofIt->first, rankToMaps[i][it->first][dofIt->first], i);
	      }

	      if (rankToMaps[i][it->first].count(dofIt->second) == 1) {
		MSG("[DBG]    %d -> %d in rank %d\n", 
		    dofIt->second, rankToMaps[i][it->first][dofIt->second], i);
	      }
	    }
	    
	    ERROR("Wrong periodic DOFs!\n");
	    foundError = 1;
	  }
	}
      }
    }

    mpi::globalAdd(foundError);
    TEST_EXIT(foundError == 0)("Error found on at least on rank!\n");


    // === 3. check: On all edge and face periodic DOFs, at least one        ===
    // === componant of coordinates of each periodic DOF pair must be equal  ===
    // === (at least as long we consider periodic boundaries only on         ===
    // === rectangulars and boxes.                                           ===

    RankToCoords sendCoords;
    map<int, vector<BoundaryType> > rankToDofType;

    for (RankToBoundMap::iterator it = pdb.intBoundary[0].periodic.begin();
	 it != pdb.intBoundary[0].periodic.end(); ++it) {
      if (it->first == mpiRank)
	continue;

      for (vector<AtomicBoundary>::iterator boundIt = it->second.begin();
	   boundIt != it->second.end(); ++boundIt) {	
	if (boundIt->rankObj.subObj == VERTEX)
	  continue;

	DofContainer dofs;
	boundIt->rankObj.el->getAllDofs(feSpace, boundIt->rankObj, dofs);
	
	for (unsigned int i = 0; i < dofs.size(); i++) {
	  WorldVector<double> c;
	  pdb.macroMesh->getDofIndexCoords(*(dofs[i]), feSpace, c);
	  sendCoords[it->first].push_back(c);
	  rankToDofType[it->first].push_back(boundIt->type);
	}
      }
    }

    // Each rank must receive exactly the same number of coordinates as it sends
    // to another rank.
    RankToCoords recvCoords;
    for (RankToCoords::iterator it = sendCoords.begin(); 
	 it != sendCoords.end(); ++it)
      recvCoords[it->first].resize(it->second.size());


    StdMpi<CoordsVec> stdMpiCoords(mpiComm, true);
    stdMpiCoords.send(sendCoords);
    stdMpiCoords.recv(recvCoords);   
    stdMpiCoords.startCommunication();


    for (RankToCoords::iterator it = sendCoords.begin(); 
	 it != sendCoords.end(); ++it) {
      for (unsigned int i = 0; i < it->second.size(); i++) {
	WorldVector<double> &c0 = it->second[i];
	WorldVector<double> &c1 = stdMpiCoords.getRecvData(it->first)[i];


	int nEqual = 0;
	for (int j = 0; j < pdb.macroMesh->getDim(); j++)
	  if (c0[j] == c1[j])
	    nEqual++;

	if (nEqual == 0) {
	  MSG("[DBG]  %d-ith periodic DOF in boundary between ranks %d <-> %d is not correct!\n",
	      i, mpiRank, it->first);
	  MSG("[DBG]  Coords on rank %d: %f %f %f\n", 
	      mpiRank, c0[0], c0[1], (pdb.macroMesh->getDim() == 3 ? c0[2] : 0.0));
	  MSG("[DBG]  Coords on rank %d: %f %f %f\n", 
	      it->first, c1[0], c1[1], (pdb.macroMesh->getDim() == 3 ? c1[2] : 0.0));

	  foundError = 1;
	}
      }
    }

    mpi::globalAdd(foundError);
    TEST_EXIT(foundError == 0)("Wrond periodic coordinates found on at least on rank!\n");
  }


  // Test if the coordinates of recv and send dofs in dofcomm are matched
  void ParallelDebug::testCommonDofs(MeshDistributor &pdb, Mesh* mesh, bool printCoords)
  {
    FUNCNAME("ParallelDebug::testCommonDofs()");

    Timer t;
    // Get FE space with basis functions of the highest degree
    const FiniteElemSpace *feSpace = 
      pdb.meshToFeSpaces[mesh][pdb.meshToFeSpaces[mesh].size() - 1];
    MultiLevelDofComm& dofComm = pdb.dofComms[mesh];

    int testCommonDofs = 1;
    Parameters::get("dbg->test common dofs", testCommonDofs);
    if (testCommonDofs == 0) {
      MSG("Skip test common dofs!\n");
      return;
    }

    int nLevels = pdb.levelData.getNumberOfLevels();
    for (int level = 0; level < nLevels; level++) {
      MPI::Intracomm &mpiComm = pdb.levelData.getMpiComm(level);
      if (mpiComm == MPI::COMM_SELF)
	continue;

      int mpiRank = mpiComm.Get_rank();
      int mpiSize = mpiComm.Get_size();
      std::set<int> &ranks = pdb.levelData.getLevelRanks(level);

      TEST_EXIT(mpiSize == static_cast<int>(ranks.size()))
	("Wrong mpi sizes:  Get_size() = %d   ranks.size() = %d\n",
	 mpiSize, ranks.size());

      /// Defines a mapping type from rank numbers to sets of DOFs.
//       typedef map<int, DofContainer> RankToDofContainer;
      
      // Maps to each neighbour rank an array of WorldVectors. This array contains the 
      // coordinates of all DOFs this rank shares on the interior boundary with the 
      // neighbour rank. A rank sends the coordinates to another rank, if it owns the
      // boundarys DOFs.
      RankToCoords sendCoords;
      
      // A rank receives all boundary DOFs that are at its interior boundaries but are
      // not owned by the rank. This map stores for each rank the coordinates of DOFs
      // this rank expectes to receive from.
      RankToCoords recvCoords;

      DOFVector<WorldVector<double> > coords(feSpace, "dofCorrds");
      mesh->getDofIndexCoords(coords);

      for (DofComm::Iterator it(dofComm[level].getSendDofs(), feSpace); 
	   !it.end(); it.nextRank())
	for (; !it.endDofIter(); it.nextDof())
	  sendCoords[it.getRank()].push_back(coords[it.getDofIndex()]);

      for (DofComm::Iterator it(dofComm[level].getRecvDofs(), feSpace); 
	   !it.end(); it.nextRank())
	for (; !it.endDofIter(); it.nextDof())
	  recvCoords[it.getRank()].push_back(coords[it.getDofIndex()]);
      
      map<int, int> sendSize;
      map<int, int> recvSize;
      map<int, int> recvSizeBuffer;
      MPI::Request *request = new MPI::Request[(mpiSize - 1) * 2];
      int requestCounter = 0;

      for (RankToCoords::iterator it = sendCoords.begin(); it != sendCoords.end(); ++it)
	sendSize[it->first] = it->second.size();

      for (RankToCoords::iterator it = recvCoords.begin(); it != recvCoords.end(); ++it)
	recvSize[it->first] = it->second.size();

      for (int rank = 0; rank < mpiSize; rank++) {
	if (rank == mpiRank)
	  continue;

	request[requestCounter++] = 
	  mpiComm.Isend(&(sendSize[rank]), 1, MPI_INT, rank, 0);
      }   

      for (int rank = 0; rank < mpiSize; rank++) {
	if (rank == mpiRank)
	  continue;

	request[requestCounter++] = 
	  mpiComm.Irecv(&(recvSizeBuffer[rank]), 1, MPI_INT, rank, 0);
      }
      
      MPI::Request::Waitall(requestCounter, request);

      int foundError = 0;
      for (std::set<int>::iterator it = ranks.begin(); it != ranks.end(); ++it) {
	if (*it == mpiRank)
	  continue;
	
	if (recvSize[*it] != recvSizeBuffer[*it]) {
	  ERROR("MPI rank %d expectes to receive %d DOFs from rank %d. But this rank sends %d DOFs!\n", 
		mpiRank, recvSize[*it], *it, recvSizeBuffer[*it]);	
	  foundError = 1;
	}
      }
      mpi::globalAdd(foundError);
      TEST_EXIT(foundError == 0)("Error found on at least on rank!\n");
      
      delete[] request;
      
      // === Now we know that the number of send and received DOFs fits together. ===
      // === So we can check if also the coordinates of the communicated DOFs are ===
      // === the same on both corresponding ranks.                                ===
      
      StdMpi<CoordsVec> stdMpi(mpiComm, true);
      stdMpi.send(sendCoords);
      stdMpi.recv(recvCoords);   
      stdMpi.startCommunication();
      
      int dimOfWorld = Global::getGeo(WORLD);
      
      // === Compare the received with the expected coordinates. ===
      
      for (RankToCoords::iterator it = stdMpi.getRecvData().begin(); 
	   it != stdMpi.getRecvData().end(); ++it) {
	for (unsigned int i = 0; i < it->second.size(); i++) {
	  WorldVector<double> tmp = (it->second)[i];
	  tmp -= recvCoords[it->first][i];
	  
	  if (norm(tmp) > 1e-8) {
	    // === Print error message if the coordinates are not the same. ===
	    if (printCoords) {
	      MSG("[DBG] i = %d\n", i);	  
	      stringstream oss;
	      oss.precision(5);
	      oss << "[DBG] Rank " << mpiRank << " from rank " << it->first
		  << " expect coords (";
	      for (int k = 0; k < dimOfWorld; k++) {
		oss << recvCoords[it->first][i][k];
		if (k + 1 < dimOfWorld)
		  oss << " / ";
	      }
	      oss << ")  received coords (";
	      for (int k = 0; k < dimOfWorld; k++) {
		oss << (it->second)[i][k];
		if (k + 1 < dimOfWorld)
		  oss << " / ";
	      }
	      oss << ")";
	      MSG("%s\n", oss.str().c_str());
	      
	      debug::printInfoByDof(feSpace, 
				    *(dofComm[level].getRecvDofs()[it->first][feSpace][i]));
	    }
	    ERROR("Wrong DOFs in rank %d!\n", mpiRank);
	    foundError = 1;
	  }	 
	}
      }
      mpi::globalAdd(foundError);
      TEST_EXIT(foundError == 0)("Error found on at least on rank!\n");
    }
      
    MSG("Test common dofs needed %.5f seconds\n", t.elapsed());
  }

  // Test dofmap dofs with same global index have same coords
  void ParallelDebug::testGlobalIndexByCoords(MeshDistributor &pdb, Mesh* mesh)
  {
    FUNCNAME("ParallelDebug::testGlobalIndexByCoords()");

    // Get FE space with basis functions of the highest degree
    MultiLevelDofComm& dofComm = pdb.dofComms[mesh];
    
    for (size_t i = 0; i < pdb.dofMaps.size(); i++) {
      vector<const FiniteElemSpace*> dofMapSpaces = pdb.dofMaps[i]->getFeSpaces(mesh);
      
      if(dofMapSpaces.empty())
	continue;
      
      const FiniteElemSpace *feSpace = NULL;
      ParallelDofMapping* dfmap = pdb.dofMaps[i];
      
      size_t j = pdb.meshToFeSpaces[mesh].size() - 1;
      for( ; j >= 0; j--)
	if(find(dofMapSpaces.begin(), dofMapSpaces.end(), pdb.meshToFeSpaces[mesh][j]) 
	  != dofMapSpaces.end())
	  break;
      feSpace = pdb.meshToFeSpaces[mesh][j];
      
      TEST_EXIT(dfmap || feSpace)("Something is wrong.\n");

      DOFVector<WorldVector<double> > coords(feSpace, "tmp");
      mesh->getDofIndexCoords(coords);

      int nLevels = pdb.levelData.getNumberOfLevels();
      for (int level = 0; level < nLevels; level++) {
	MPI::Intracomm &mpiComm = pdb.levelData.getMpiComm(level);
	if (mpiComm == MPI::COMM_SELF)
	  continue;

  //       int mpiRank = mpiComm.Get_rank();
  //       int mpiSize = mpiComm.Get_size();

	typedef map<int, WorldVector<double> > CoordsIndexMap;
	CoordsIndexMap coordsToIndex;
	
	DOFIterator<WorldVector<double> > it(&coords, USED_DOFS);
	for (it.reset(); !it.end(); ++it) {
	  int idx = (*dfmap)[feSpace][it.getDOFIndex()].global;
	  coordsToIndex[idx] = *it;
  // 	MSG("   CHECK FOR DOF %d/%d AT COORDS %f %f %f\n",
  // 	    it.getDOFIndex(),
  // 	    coordsToIndex[(*it)], 
  // 	    (*it)[0], (*it)[1], (pdb.mesh->getDim() == 3 ? (*it)[2] : 0.0));
	}
	
	StdMpi<CoordsIndexMap> stdMpi(mpiComm, true);
	for (DofComm::Iterator it(dofComm[level].getSendDofs(), feSpace); 
	    !it.end(); it.nextRank())
	  stdMpi.send(it.getRank(), coordsToIndex);
	for (DofComm::Iterator it(dofComm[level].getRecvDofs(), feSpace); 
	    !it.end(); it.nextRank())
	  stdMpi.recv(it.getRank());
	
	stdMpi.startCommunication();
	
	int foundError = 0;
	for (DofComm::Iterator it(dofComm[level].getRecvDofs(), feSpace); 
	    !it.end(); it.nextRank()) {
	  CoordsIndexMap& otherCoords = stdMpi.getRecvData(it.getRank());

	  for (; !it.endDofIter(); it.nextDof()) {
	    WorldVector<double> recvCoords = coords[it.getDofIndex()];
	    int idx = (*dfmap)[feSpace][it.getDofIndex()].global;
	    
	    TEST_EXIT_DBG(otherCoords.count(idx))("Global index not found in neighbour partition\n");

	    WorldVector<double> diff = otherCoords[idx] - recvCoords;
	    double dist = sqrt(diff*diff);
	    if (dist > DBL_TOL) {
	      stringstream oss;
	      oss.precision(5);
	      oss << "DOF at coords ";
	      for (int i = 0; i < Global::getGeo(WORLD); i++)
		oss << recvCoords[i] << " ";
	      oss << " do not fit together on rank " 
		  << mpiComm.Get_rank() << " (global index: " 
		  << it.getDofIndex() << ") and on rank "
		  << it.getRank() << " "
		  << "  - LEVEL " << level;
	      
	      MSG("[DBG] %s\n", oss.str().c_str());
	      
	      foundError = 1;	    
	    }
	  }
	}

	mpi::globalAdd(foundError);
	TEST_EXIT(foundError == 0)("Error found on at least on rank!\n");
      }
    }
  }


  void ParallelDebug::testAllElements(MeshDistributor &pdb)
  {
    FUNCNAME("ParallelDebug::testAllElements()");
    
    TEST_EXIT(pdb.macroMesh) ("No macro mesh in mesh distributor.\n");

    std::set<int> macroElements;
    int minElementIndex = numeric_limits<int>::max();
    int maxElementIndex = numeric_limits<int>::min();   

    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(pdb.macroMesh, 0, Mesh::CALL_EL_LEVEL);
    while (elInfo) {
      int elIndex = elInfo->getElement()->getIndex();
      minElementIndex = std::min(minElementIndex, elIndex);
      maxElementIndex = std::max(maxElementIndex, elIndex);
      macroElements.insert(elInfo->getElement()->getIndex());
      elInfo = stack.traverseNext(elInfo);      
    }
      

    MPI::Intracomm mpiComm = MPI::COMM_WORLD;

    int globalMinIndex, globalMaxIndex;
    mpiComm.Allreduce(&minElementIndex, &globalMinIndex, 1, MPI_INT, MPI_MIN);
    mpiComm.Allreduce(&maxElementIndex, &globalMaxIndex, 1, MPI_INT, MPI_MAX);

    TEST_EXIT(globalMinIndex == 0)("No macro element with index 0!\n");
    for (int i = 0; i <= globalMaxIndex; i++) {
      int sendId = macroElements.count(i);
      int recvId = 0;
      mpiComm.Allreduce(&sendId, &recvId, 1, MPI_INT, MPI_SUM);

      if (recvId != 1 && mpiComm.Get_rank() == 0) {
	if (recvId == 0) {
	  ERROR_EXIT("Element %d has no member partition!\n", i);
	}

	if (recvId > 1) {
	  ERROR_EXIT("Element %d is member of more than pne partition!\n", i); 
	}
      }
    }
  }

  // Test dofcomm whether send and recv dofs match
  void ParallelDebug::testDofContainerCommunication(MeshDistributor &pdb, Mesh* mesh)
  {
    FUNCNAME("ParallelDebug::testDofContainerCommunication()");    

    MSG("WARNING: ONLY LEVEL 0 TEST!\n");

    typedef map<int, map<const FiniteElemSpace*, DofContainer> >::iterator it_type;

    MPI::Intracomm &mpiComm = pdb.levelData.getMpiComm(0);
    MultiLevelDofComm& dofComm = pdb.dofComms[mesh];
    
    map<int, int> sendNumber;
    for (it_type it = dofComm[0].getSendDofs().begin(); 
	 it != dofComm[0].getSendDofs().end(); ++it)
      for (map<const FiniteElemSpace*, DofContainer>::iterator dcIt = it->second.begin(); 
	   dcIt != it->second.end(); ++dcIt)
	sendNumber[it->first] += dcIt->second.size();
       
    map<int, int> recvNumber;
    for (it_type it = dofComm[0].getRecvDofs().begin(); 
	 it != dofComm[0].getRecvDofs().end(); ++it)
      for (map<const FiniteElemSpace*, DofContainer>::iterator dcIt = it->second.begin(); 
	   dcIt != it->second.end(); ++dcIt)
	recvNumber[it->first] += dcIt->second.size();
    
    StdMpi<int> stdMpi(mpiComm);
    stdMpi.send(sendNumber);
    for (it_type it = dofComm[0].getRecvDofs().begin(); 
	 it != dofComm[0].getRecvDofs().end(); ++it)
      stdMpi.recv(it->first);

    stdMpi.startCommunication();

    int foundError = 0;
    for (map<int, int>::iterator it = stdMpi.getRecvData().begin();
	 it != stdMpi.getRecvData().end(); ++it) {
      if (it->second != recvNumber[it->first]) {
	ERROR("Rank expectes %d DOFs to receive from rank %d, but %d DOFs are received!\n", 
	      recvNumber[it->first], it->first, it->second);
	foundError = 1;
      }
    }
    
    std::set<DegreeOfFreedom> sendDofs;
    for (DofComm::Iterator it(dofComm[0].getSendDofs(), pdb.meshToFeSpaces[mesh][0]);
       !it.end(); it.nextRank())
      for (; !it.endDofIter(); it.nextDof())
        sendDofs.insert(it.getDofIndex());

    for (DofComm::Iterator it2(dofComm[0].getRecvDofs(), pdb.meshToFeSpaces[mesh][0]);
       !it2.end(); it2.nextRank())
      for (; !it2.endDofIter(); it2.nextDof())
        TEST_EXIT(!sendDofs.count(it2.getDofIndex()))
          ("Send and recv dof contaniners share same dof %d.\n", it2.getDofIndex());

    mpi::globalAdd(foundError);
    TEST_EXIT(foundError == 0)("Error found on at least one rank!\n");
  }

  // Test whether one domain has two dofs with same coordinates
  void ParallelDebug::testDoubleDofs(Mesh *mesh)
  {
    FUNCNAME("ParallelDebug::testDoubleDofs()");

    map<WorldVector<double>, DegreeOfFreedom> cMap;
    int foundError = 0;

    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL | Mesh::FILL_COORDS);
    while (elInfo) {
      for (int i = 0; i < mesh->getGeo(VERTEX); i++) {
	WorldVector<double> &c = elInfo->getCoord(i);
	if (cMap.count(c) == 0) {
	  cMap[c] = elInfo->getElement()->getDof(i, 0);
	} else {
	  if (cMap[c] != elInfo->getElement()->getDof(i, 0)) {
	    MSG("[DBG] Found two DOFs %d and %d with the same coords %f %f %f!\n",
		cMap[c], elInfo->getElement()->getDof(i, 0), 
		c[0], c[1], mesh->getDim() == 3 ? c[2] : 0.0);
	    foundError = 1;
	  }
	}
      }
      
      elInfo = stack.traverseNext(elInfo);
    }
    
    mpi::globalAdd(foundError);
    TEST_EXIT(foundError == 0)("Error found on at least one rank!\n");
  }


  void ParallelDebug::printMapLocalGlobal(MeshDistributor &pdb, int rank)
  {   
    FUNCNAME("ParallelDebug::printMapLocalGlobal()");

    ERROR_EXIT("Rewrite this function!\n");

#if 0
    if (rank == -1 || MPI::COMM_WORLD.Get_rank() == rank) {
      const FiniteElemSpace *feSpace = pdb.feSpaces[0];

      cout << "====== DOF MAP LOCAL -> GLOBAL ====== " << endl;
            
      DofMap &dofMap = pdb.dofMap[feSpace].getMap();
      for (DofMap::iterator it = dofMap.begin(); it != dofMap.end(); ++it) {
	cout << "DOF " << it->first << " " << it->second.global << "\n";
	WorldVector<double> coords;
	pdb.mesh->getDofIndexCoords(it->first, feSpace, coords);
	coords.print();

	for (DofComm::Iterator rit(pdb.dofComm.getSendDofs(), feSpace); 
	     !rit.end(); rit.nextRank())
	  for (; !rit.endDofIter(); rit.nextDof())
	    if (it->first == rit.getDofIndex())
	      cout << "SEND DOF TO " << rit.getRank() << endl;
	
	for (DofComm::Iterator rit(pdb.dofComm.getRecvDofs(), feSpace); 
	     !rit.end(); rit.nextRank())
	  for (; !rit.endDofIter(); rit.nextDof())
	    if (it->first == rit.getDofIndex())
	      cout << "RECV DOF FROM " << rit.getRank() << endl;

	cout << "------" << endl;
      }
    }
#endif
  }


  void ParallelDebug::printMapPeriodic(MeshDistributor &pdb, int rank)
  {
    FUNCNAME("ParallelDebug::printMapPeriodic()");

    ERROR_EXIT("Function must be rewritten, check svn for old code!\n");
  }

  
  void ParallelDebug::printRankDofs(MeshDistributor &pdb, 
				    int rank, 
				    DofContainer& rankDofs,
				    DofContainer& rankAllDofs)
  {
    if (rank == -1 || MPI::COMM_WORLD.Get_rank() == rank) {
      cout << "====== RANK DOF INFORMATION ====== " << endl;

      cout << "  RANK OWNED DOFS: " << endl;
      for (DofContainer::iterator dofit = rankDofs.begin();
	   dofit != rankDofs.end(); ++dofit) {
	cout << "    " << **dofit << endl;
	WorldVector<double> coords;
	pdb.macroMesh->getDofIndexCoords(*dofit, pdb.feSpaces[0], coords);
	coords.print();
      }

      cout << "  RANK ALL DOFS: " << endl;
      for (DofContainer::iterator dofit = rankAllDofs.begin();
	   dofit != rankAllDofs.end(); ++dofit) {
	cout << "    " << **dofit << endl;
	WorldVector<double> coords;
	pdb.macroMesh->getDofIndexCoords(*dofit, pdb.feSpaces[0], coords);
	coords.print();
      }      
    }
  }


  void ParallelDebug::printBoundaryInfo(InteriorBoundary &intBoundary,
					bool force)
  {
    FUNCNAME("ParallelDebug::printBoundaryInfo()");

    int tmp = 0;
    Parameters::get("parallel->debug->print boundary info", tmp);
    if (tmp <= 0 && force == false)
      return;

    MSG("Interior boundary info:\n");

    for (InteriorBoundary::iterator it(intBoundary.own); !it.end(); ++it) {
      MSG("Rank owned boundary with rank %d: \n", it.getRank());
      MSG("  ranks obj-ind: %d  sub-obj: %d   ith-obj: %d\n",
	  it->rankObj.elIndex, it->rankObj.subObj, it->rankObj.ithObj);
      MSG("  neigh obj-ind: %d  sub-obj: %d   ith-obj: %d\n",
	  it->neighObj.elIndex, it->neighObj.subObj, it->neighObj.ithObj);
    }

    for (InteriorBoundary::iterator it(intBoundary.other); !it.end(); ++it) {
      MSG("Other owned boundary with rank %d: \n", it.getRank());
      MSG("  ranks obj-ind: %d  sub-obj: %d   ith-obj: %d\n",
	  it->rankObj.elIndex, it->rankObj.subObj, it->rankObj.ithObj);
      MSG("  neigh obj-ind: %d  sub-obj: %d   ith-obj: %d\n",
	  it->neighObj.elIndex, it->neighObj.subObj, it->neighObj.ithObj);
    }

    for (InteriorBoundary::iterator it(intBoundary.periodic); !it.end(); ++it) {
      MSG("Periodic boundary (ID %d) with rank %d: \n", 
	  it->type, it.getRank());
      MSG("  ranks obj-ind: %d  sub-obj: %d   ith-obj: %d\n",
	  it->rankObj.elIndex, it->rankObj.subObj, it->rankObj.ithObj);
      MSG("  neigh obj-ind: %d  sub-obj: %d   ith-obj: %d\n",
	  it->neighObj.elIndex, it->neighObj.subObj, it->neighObj.ithObj);
    }    
  }
  
  void ParallelDebug::writeDebugFile(MeshToFeSpaces& meshToFeSpaces, 
				     ParallelDofMapping &dofMap,
				     string debugOutputDir)
  {
    int i = 0;
    MeshToFeSpaces::iterator it = meshToFeSpaces.begin();
    while(it != meshToFeSpaces.end()) {
      string prefix = debugOutputDir + "mpi-dbg-mesh" + std::to_string(i);
      writeDebugFile(it->second[it->second.size() - 1],
		     dofMap, prefix, "dat");
      it++;
      i++;
    }
    
  }

  void ParallelDebug::writeDebugFile(const FiniteElemSpace *feSpace,
				     ParallelDofMapping &dofMap,
				     string prefix, 
				     string postfix)
  {
    FUNCNAME("ParallelDebug::writeDebugFile()");

    Mesh *mesh = feSpace->getMesh();

    stringstream filename;
    filename << prefix << "-" << MPI::COMM_WORLD.Get_rank() << "." << postfix;

    DOFVector<WorldVector<double> > coords(feSpace, "tmp");
    mesh->getDofIndexCoords(coords);

    typedef map<int, vector<DegreeOfFreedom> > ElDofMap;
    ElDofMap elDofMap;
    TraverseStack stack;
    const BasisFunction *basisFcts = feSpace->getBasisFcts();
    vector<DegreeOfFreedom> localIndices(basisFcts->getNumber());
    ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_LEAF_EL);
    while (elInfo) {
      basisFcts->getLocalIndices(elInfo->getElement(), 
				 feSpace->getAdmin(), localIndices);
      elDofMap[elInfo->getElement()->getIndex()] = localIndices;
      elInfo = stack.traverseNext(elInfo);
    }

    // === Write informations about all DOFs. ===

    ofstream file;
    file.open(filename.str().c_str());
//     file << "# First line contains number of DOFs, than each line has the format\n";
//     file << "# DOF index     Local DOF index     Global DOF index     Is rank DOF     x-coord     y-coord     z-coord\n";
    file << dofMap[feSpace].size() << "\n";
    DOFIterator<WorldVector<double> > it(&coords, USED_DOFS);
    for (it.reset(); !it.end(); ++it) {
      if (dofMap[feSpace].isSet(it.getDOFIndex())) {
	file << it.getDOFIndex() << " " 
	     << dofMap[feSpace][it.getDOFIndex()].local << " "
	     << dofMap[feSpace][it.getDOFIndex()].global << " "
	     << dofMap[feSpace].isRankDof(it.getDOFIndex());
	for (int i = 0; i < mesh->getDim(); i++)
	  file << " " << (*it)[i];
	file << "\n";
      }
    }

    // === Write to all elements in ranks mesh the included dofs. ===

//     file << "\n\n";
//     file << "# First line containes number of elements in mesh, second line contain the number of DOFs per element.\n";
//     file << "# Than, each entry contains of two lines. The first is the element index, the second line is a list with the local DOF indices of this element.\n";
    file << elDofMap.size() << "\n";
    file << basisFcts->getNumber() << "\n";
    for (ElDofMap::iterator it = elDofMap.begin(); it != elDofMap.end(); ++it) {
      file << it->first << "\n";
      for (int i = 0; i < basisFcts->getNumber(); i++) 
	file << it->second[i] << " ";
      file << "\n";
    }

    file.close();
  }


  void ParallelDebug::writePartitioning(MeshDistributor &pdb, string filename)
  {
    FUNCNAME("ParallelDebug::writeParitioning()");

    map<int, double> vec;    
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(pdb.macroMesh, -1, Mesh::CALL_LEAF_EL);
    
    while (elInfo) {		  
      int index = elInfo->getMacroElement()->getIndex();
      vec[index] = pdb.partitionMap[index];
      elInfo = stack.traverseNext(elInfo);
    }

    io::ElementFileWriter::writeFile(vec, pdb.macroMesh, filename);
  }


  void ParallelDebug::writePartitioningFile(string filename, 
					    int counter,
					    const FiniteElemSpace *feSpace)
  {
    FUNCNAME("ParallelDebug::writePartitioningFile()");

    stringstream oss;
    oss << filename;
    if (counter >= 0)
      oss << "-" << counter;
    oss << ".vtu";

    DOFVector<double> tmp(feSpace, "tmp");
    tmp.set(MPI::COMM_WORLD.Get_rank());
    io::VtkWriter::writeFile(&tmp, oss.str());
  }

  
  bool ParallelDebug::followThisBound(int rankElIndex, int neighElIndex)
  {
    FUNCNAME("ParallelDebug::followThisBound()");

    int el0 = std::min(rankElIndex, neighElIndex);
    int el1 = std::max(rankElIndex, neighElIndex);

    vector<int> els;
    Parameters::get("parallel->debug->follow boundary", els);
    if (els.size() != 2)
      return false;

    return (el0 == els[0] && el1 == els[1]);
  }


  void ParallelDebug::followBoundary(MeshDistributor &pdb)
  {
    FUNCNAME("ParallelDebug::followBoundary()");

    int mpiRank = MPI::COMM_WORLD.Get_rank();

    for (InteriorBoundary::iterator it(pdb.intBoundary[0].own); !it.end(); ++it)
      if (followThisBound(it->rankObj.elIndex, it->neighObj.elIndex))
	debug::writeLocalElementDofs(mpiRank, 
				     it->rankObj.elIndex, 
				     pdb.feSpaces[0]);

    for (InteriorBoundary::iterator it(pdb.intBoundary[0].other); !it.end(); ++it)
      if (followThisBound(it->rankObj.elIndex, it->neighObj.elIndex))
	debug::writeLocalElementDofs(mpiRank,
				     it->rankObj.elIndex,
				     pdb.feSpaces[0]);
  }


  void ParallelDebug::followBoundary(Mesh *mesh, 
				     AtomicBoundary &bound, 
				     MeshStructure &code)
  {
    FUNCNAME("ParallelDebug::followBoundary()");

    if (mesh->getDim() != bound.rankObj.subObj)
      return;

    if (!followThisBound(bound.rankObj.elIndex, bound.neighObj.elIndex))
      return;

    MSG("Mesh structure code of bound %d/%d <-> %d/%d: %s\n",
	bound.rankObj.elIndex, mesh->getDim(),
	bound.neighObj.elIndex, mesh->getDim(), 
	code.toStr().c_str());
  }


#if 0
  //  compare PETSc matrices

    if (MPI::COMM_WORLD.Get_rank() == 0) {
      Mat matOld, matNew;
      MatCreate(PETSC_COMM_SELF, &matOld);
      MatCreate(PETSC_COMM_SELF, &matNew);

      PetscViewer fdOld, fdNew;
      PetscViewerBinaryOpen(PETSC_COMM_SELF, "coarse_interior_old.mat", FILE_MODE_READ, &fdOld);
      PetscViewerBinaryOpen(PETSC_COMM_SELF, "coarse_interior_new.mat", FILE_MODE_READ, &fdNew);
      MatLoad(matOld, fdOld);
      MatLoad(matNew, fdNew);

      int m, n;
      MatGetSize(matOld, &m, &n);
      MSG("MAT OLD SIZE: %d %d\n", m, n);

      MatGetSize(matNew, &m, &n);
      MSG("MAT NEW SIZE: %d %d\n", m, n);

      for (int i = 0; i < m; i++) {
	PetscInt nColsOld, nColsNew;
	const PetscInt *colsOld, *colsNew;
	const PetscScalar *valsOld, *valsNew;

	MatGetRow(matOld, i, &nColsOld, &colsOld, &valsOld);
	MatGetRow(matNew, i, &nColsNew, &colsNew, &valsNew);

	MSG("  row: %d with cols %d %d\n", i, nColsOld, nColsNew);

	if (nColsOld != nColsNew) {
	  MSG("WRONG COLS NUMBER in ROW %d: %d %d\n", i, nColsOld, nColsNew);
	} else {
	  for (int j = 0; j < nColsOld; j++) {
	    if (colsOld[j] != colsNew[j]) {
	      MSG("WRONG COLS IN ROW %d: %d %d \n", i, colsOld[j], colsNew[j]);
	    } else {
	      if (fabs(colsOld[j] - colsNew[j]) > 1e-8) {
		MSG("WRONG VALUES IN ROW: %d\n", i);
	      }
	    }
	  }
	}

	MatRestoreRow(matOld, i, &nColsOld, &colsOld, &valsOld);
	MatRestoreRow(matNew, i, &nColsNew, &colsNew, &valsNew);
      }

      MSG("FINISHED-TEST\n");
    }
#endif


  void ParallelDebug::writeCsvElementMap(const FiniteElemSpace *feSpace,
					 ParallelDofMapping &dofMap,
					 string prefix, 
					 string postfix)
  {
    FUNCNAME("ParallelDebug::writeCsvElementMap()");

    return;

    MSG("writing local Element map to CSV File \n");

    Mesh *mesh = feSpace->getMesh();

    if (mesh->getDim() != 3)
      return;

    stringstream filename;
    filename << prefix << "-" << MPI::COMM_WORLD.Get_rank() << "." << postfix;

    ofstream file;
    file.open(filename.str().c_str());

    DOFVector<WorldVector<double> > dofCoords(feSpace, "DOF coords");
    mesh->getDofIndexCoords(dofCoords);

    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(mesh,	
					 -1, 
					 Mesh::CALL_EVERY_EL_POSTORDER);

    while (elInfo) {
      /*
	MSG("Start traverse element %d in level %d\n",
	elInfo->getElement()->getIndex(),
	elInfo->getLevel());
      */

      // check if Element is Leaflevel
      // because higherOrderDofs(=NON VERTICES DOFS) are not saved in nonleaf Elements
      if (elInfo->getElement()->isLeaf()) {

	//traverse ALL DOFS
	ElementDofIterator elDofIter(feSpace, true);
	elDofIter.reset(elInfo->getElement());
	
	int locDofNumber = 0;  

	do {
   	
	  WorldVector<double> c = dofCoords[elDofIter.getDof()];
	
	  file << elInfo->getElement()->getIndex() << "\t";	//element number
	  file << elInfo->getLevel() << "\t";					//element Level
	  file << elDofIter.getDof() << "\t";					//dof number
	  file << elDofIter.getCurrentPos()+1 << "\t";		//dof type (+1 to pseudo convert to geoIndex, does not work for CENTER DOFS( geoIndex=0))
	  file << c[0] << "\t";								//dof coords
	  file << c[1] << "\t";	
	  file << c[2] << "\t";
	  file << locDofNumber << "\t";						//local Dof number
	  file << elInfo->getType() << "\n";					//element Type
		
	  locDofNumber++;
	} while (elDofIter.next());
      } else {
	
	// traverse only VERTEX DOFS
	for (int i = 0; i < mesh->getGeo(VERTEX); i++) {
	  DegreeOfFreedom dof = elInfo->getElement()->getDof(i, 0);	

	  WorldVector<double> c = dofCoords[dof];
		
	  file << elInfo->getElement()->getIndex() << "\t";	//element number
	  file << elInfo->getLevel() << "\t";					//element Level
	  file << dof << "\t";								//dof number
	  file << VERTEX << "\t";								//dof type
	  file << c[0] << "\t";								//dof coords
	  file << c[1] << "\t";	
	  file << c[2] << "\t";
	  file << i << "\t";									//local Dof number
	  file << elInfo->getType() << "\n";					//element Type
	}		
      }

      elInfo = stack.traverseNext(elInfo);
    }

  }
  
  void ParallelDebug::writeDofMap(ParallelDofMapping &dofMap, int component, string filename, string postfix)
  {
    filename += std::to_string(MPI::COMM_WORLD.Get_rank());
    filename = filename + "." + postfix;
     
    ofstream file;
    file.open(filename.c_str(), ios::out | ios::trunc);
	
    DofToMatIndex::MapType& data = dofMap.getMatData(component);
    
    DofToMatIndex::MapType::iterator it = data.begin();
    
    file << &dofMap << "\n";
    file << dofMap[component].nRankDofs << "\n";
    
    //
//     for (DofComm::Iterator it(dofComm->getRecvDofs(), feSpace);
// 	 !it.end(); it.nextRank()) {
//       int rank = it.getRank();
//       // if (meshLevel > 0)
//       //   rank = levelData->mapRank(rank, 0, meshLevel);
// 
//       int i = 0;
//       for (; !it.endDofIter(); it.nextDof())
// 	if (nonRankDofs.count(it.getDofIndex()))
// 	  dofMap[it.getDofIndex()].global = stdMpi.getRecvData(rank)[i++];
//     }
    //
    
    while (it != data.end()) {
      if (dofMap.isMatIndexFromGlobal())
	file << it->first << " " << it->second << "\n";
      else
        file << it->first << " " << it->second << "\n";
      it++;
    }
    file.close();
  }
  
  void ParallelDebug::writePeriodicElObjInfo(MeshDistributor &pdb, string debugOutputDir)
  {
    string filename = debugOutputDir + "periodic-info.dat";
    ofstream file;
    file.open(filename.c_str());
    
    ElementObjectDatabase& elObjDb = pdb.getElementObjectDb();
    
    file << "periodicVertices:\n";
    for (PerBoundMap<DegreeOfFreedom>::iterator it = elObjDb.periodicVertices.begin();
	 it != elObjDb.periodicVertices.end(); it++) {
      file << it->first.first << ", " << it->first.second
	   << " : " << it->second << endl;
    }
    file << "periodicEdges:\n";
    for (PerBoundMap<DofEdge>::iterator it = elObjDb.periodicEdges.begin(); 
	 it != elObjDb.periodicEdges.end(); it++) {
      file << it->first.first.first << " " << it->first.first.second << ", "
	   << it->first.second.first << " " << it->first.second.second << " : "
	   << it->second << endl;
    }
    file << "periodicFaces:\n";
    for (PerBoundMap<DofFace>::iterator it = elObjDb.periodicFaces.begin(); 
	 it != elObjDb.periodicFaces.end(); it++) {
      file << it->first.first.get<0>() << it->first.first.get<1>()
	   << it->first.first.get<2>() << ", "
	   << it->first.second.get<0>() << it->first.second.get<1>()
	   << it->first.second.get<2>() << " : " 
	   << it->second << endl;
    }
    file << "periodicDofAssoc:\n";
    {std::map<DegreeOfFreedom, std::set<BoundaryType> >::iterator it;
    for (it = elObjDb.periodicDofAssoc.begin(); 
	 it != elObjDb.periodicDofAssoc.end(); it++) {
      file << it->first << endl;
      for(std::set<BoundaryType>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
	file << *it2 << " ";
      file << endl;
    }}
    file << "periodicEdgeAssoc:\n";
    {std::map<DofEdge, std::set<DofEdge> >::iterator it;
    for (it = elObjDb.periodicEdgeAssoc.begin(); 
	it != elObjDb.periodicEdgeAssoc.end(); it++) {
      file << it->first.first << " " << it->first.second << endl;
      for(std::set<DofEdge>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
	file << "(" <<  (*it2).first << ", " << (*it2).second << ") ";
      file << endl;
    }}
    file.close();
  }
  
  void ParallelDebug::writeInterchangeVector(MeshDistributor &pdb, string debugOutputDir)
  {
    for (size_t i = 0; i < pdb.interchangeVectors.size(); i++) {
      string filename = debugOutputDir + "-int" + std::to_string(i) + ".vtu";
      io::VtkWriter::writeFile(pdb.interchangeVectors[i], filename);
    }
  }
  
} } // end namespaces
