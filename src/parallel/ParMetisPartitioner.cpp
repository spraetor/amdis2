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


#include <queue>
#include <mpi.h>

#include "parallel/ParMetisPartitioner.hpp"
#include "parallel/ParallelDofMapping.hpp"
#include "parallel/MpiHelper.hpp"
#include "Serializer.h"
#include "Mesh.h"
#include "Traverse.h"
#include "ElInfo.h"
#include "Element.h"
#include "FixVec.h"
#include "DOFVector.h"

using namespace std;

namespace AMDiS
{
  namespace Parallel
  {

    ParMetisMesh::ParMetisMesh(Mesh* mesh, MPI::Intracomm* comm,
                               std::map<int, bool>& elementInRank,
                               DofMap* mapLocalGlobal)
      : dim(mesh->getDim()),
        nElements(0),
        mpiComm(comm)
    {
      FUNCNAME("ParMetisMesh::ParMetisMesh()");

      int mpiSize = mpiComm->Get_size();
      int elementCounter = 0;
      int dow = Global::getGeo(WORLD);

      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL);
      while (elInfo)
      {
        if (elementInRank[elInfo->getElement()->getIndex()])
          elementCounter++;

        elInfo = stack.traverseNext(elInfo);
      }

      nElements = elementCounter;

      TEST_EXIT(nElements > 0)("No elements in ParMETIS mesh!\n");

      // allocate memory
      eptr = new int[nElements + 1];
      eind = new int[nElements * (mesh->getGeo(VERTEX))];
      elmdist = new int[mpiSize + 1];
      elem_p2a = new int[nElements];

      if (dim == dow)
        xyz = new float[nElements * dim];
      else
        xyz = NULL;

      eptr[0] = 0;

      int* ptr_eptr = eptr + 1;
      int* ptr_eind = eind;
      float* ptr_xyz = xyz;

      // gather element numbers and create elmdist
      mpiComm->Allgather(&nElements, 1, MPI_INT, elmdist + 1, 1, MPI_INT);

      elmdist[0] = 0;
      for (int i = 2; i < mpiSize + 1; i++)
        elmdist[i] += elmdist[i - 1];

      // traverse mesh and fill distributed ParMETIS data
      DimVec<double> bary(dim, 1.0 / mesh->getGeo(VERTEX));
      WorldVector<double> world;

      elementCounter = 0;
      int nodeCounter = 0;

      elInfo = stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL | Mesh::FILL_COORDS);
      while (elInfo)
      {
        Element* element = elInfo->getElement();
        int index = element->getIndex();

        // if element in partition
        if (elementInRank[index])
        {
          // remember index
          setParMetisIndex(index, elementCounter);
          setAMDiSIndex(elementCounter, index);

          // write eptr entry
          nodeCounter += mesh->getGeo(VERTEX);
          *ptr_eptr = nodeCounter;
          ptr_eptr++;

          // write eind entries (element nodes)
          for (int i = 0; i < dim + 1; i++)
          {
            if (mapLocalGlobal)
              *ptr_eind = (*mapLocalGlobal)[element->getDof(i, 0)].global;
            else
              *ptr_eind = element->getDof(i, 0);

            ptr_eind++;
          }

          // write xyz element coordinates
          if (ptr_xyz)
          {
            elInfo->coordToWorld(bary, world);
            for (int i = 0; i < dim; i++)
            {
              *ptr_xyz = static_cast<float>(world[i]);
              ptr_xyz++;
            }
          }

          elementCounter++;
        }

        elInfo = stack.traverseNext(elInfo);
      }
    }


    ParMetisMesh::~ParMetisMesh()
    {
      if (eptr)
        delete [] eptr;

      if (eind)
        delete [] eind;

      if (elmdist)
        delete [] elmdist;

      if (xyz)
        delete [] xyz;

      if (elem_p2a)
        delete [] elem_p2a;
    }


    ParMetisGraph::ParMetisGraph(ParMetisMesh* parMesh,
                                 MPI::Intracomm* comm,
                                 int ncommonnodes)
      : parMetisMesh(parMesh)
    {
      FUNCNAME("ParMetisGraph::ParMetisGraph()");

      TEST_EXIT(parMesh)("No ParMetisMesh defined!\n");
      TEST_EXIT(comm)("No MPI communicator defined!\n");

      int numflag = 0;

      if (ncommonnodes == -1)
        ncommonnodes = parMetisMesh->getDim();

      MPI_Comm tmpComm = MPI_Comm(*comm);

      ParMETIS_V3_Mesh2Dual(parMetisMesh->getElementDist(),
                            parMetisMesh->getElementPtr(),
                            parMetisMesh->getElementInd(),
                            &numflag,
                            &ncommonnodes,
                            &xadj,
                            &adjncy,
                            &tmpComm);
    }


    ParMetisGraph::~ParMetisGraph()
    {
      free(xadj);
      free(adjncy);
    }


    void ParMetisGraph::print()
    {
      FUNCNAME("ParMetisGraph::print()");

      stringstream oss;
      for (int i = 0; i <= MPI::COMM_WORLD.Get_size(); i++)
        oss << parMetisMesh->getElementDist()[i] << " ";

      MSG("Element dist = %s\n", oss.str().c_str());

      int mpiRank = MPI::COMM_WORLD.Get_rank();
      int nElements = parMetisMesh->getElementDist()[mpiRank + 1] -
                      parMetisMesh->getElementDist()[mpiRank];

      MSG("nElements = %d   in index range %d - %d\n",
          nElements,
          parMetisMesh->getElementDist()[mpiRank],
          parMetisMesh->getElementDist()[mpiRank + 1]);

      oss.str("");
      oss.clear();

      for (int i = 0; i <= nElements; i++)
        oss << xadj[i] << ", ";

      MSG("xadj = {%s}\n", oss.str().c_str());

      oss.str("");
      oss.clear();

      for (int i = 0; i <= xadj[nElements] - 1; i++)
        oss << adjncy[i] << ", ";

      MSG("adjncy = {%s}\n", oss.str().c_str());
    }


    ParMetisPartitioner::~ParMetisPartitioner()
    {
      if (parMetisMesh)
        delete parMetisMesh;
    }


    bool ParMetisPartitioner::partition(map<int, double>& elemWeights,
                                        PartitionMode mode)
    {
      FUNCNAME("ParMetisPartitioner::partition()");

      int mpiSize = mpiComm->Get_size();


      // === Create parmetis mesh ===

      if (parMetisMesh)
        delete parMetisMesh;

      TEST_EXIT_DBG(elementInRank.size() != 0)("Should not happen!\n");

      parMetisMesh = new ParMetisMesh(mesh, mpiComm, elementInRank, mapLocalGlobal);

      int nElements = parMetisMesh->getNumElements();


      // === Create weight array ===

      vector<int> wgts(nElements);
      vector<float> floatWgts(nElements);
      unsigned int floatWgtsPos = 0;
      float maxWgt = 0.0;

      TraverseStack stack;
      ElInfo* elInfo = stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL);
      while (elInfo)
      {
        int index = elInfo->getElement()->getIndex();

        if (elementInRank[index])
        {
          // get weight
          float wgt = static_cast<float>(elemWeights[index]);
          maxWgt = std::max(wgt, maxWgt);

          // write float weight
          TEST_EXIT_DBG(floatWgtsPos < floatWgts.size())("Should not happen!\n");
          floatWgts[floatWgtsPos++] = wgt;
        }
        elInfo = stack.traverseNext(elInfo);
      }

      TEST_EXIT_DBG(floatWgtsPos == floatWgts.size())("Should not happen!\n");

      float tmp;
      mpiComm->Allreduce(&maxWgt, &tmp, 1, MPI_FLOAT, MPI_MAX);
      maxWgt = tmp;


      // === Create dual graph ===

      ParMetisGraph parMetisGraph(parMetisMesh, mpiComm);


      // === Partitioning of dual graph ===

      int wgtflag = 2; // weights at vertices only!
      int numflag = 0; // c numbering style!
      int ncon = 1; // one weight at each vertex!
      int nparts = mpiSize; // number of partitions

      vector<double> tpwgts(mpiSize);
      double ubvec = 1.05;
      int options[4] = {0, 0, 15, PARMETIS_PSR_COUPLED}; // default options
      int edgecut = -1;
      vector<int> part(nElements);

      // set tpwgts
      for (int i = 0; i < mpiSize; i++)
        tpwgts[i] = 1.0 / static_cast<double>(nparts);

      //     float scale = 10000.0 / maxWgt;
      for (int i = 0; i < nElements; i++)
        wgts[i] = floatWgts[i];
      //      wgts[i] = static_cast<int>(floatWgts[i] * scale);


      // === Start ParMETIS. ===

      MPI_Comm tmpComm = MPI_Comm(*mpiComm);

      switch (mode)
      {
      case INITIAL:
        ParMETIS_V3_PartKway(parMetisMesh->getElementDist(),
                             parMetisGraph.getXAdj(),
                             parMetisGraph.getAdjncy(),
                             &(wgts[0]),
                             NULL,
                             &wgtflag,
                             &numflag,
                             &ncon,
                             &nparts,
                             &(tpwgts[0]),
                             &ubvec,
                             options,
                             &edgecut,
                             &(part[0]),
                             &tmpComm);
        break;
      case ADAPTIVE_REPART:
      {
        vector<int> vsize(nElements);
        for (int i = 0; i < nElements; i++)
          vsize[i] = static_cast<int>(floatWgts[i]);

        ParMETIS_V3_AdaptiveRepart(parMetisMesh->getElementDist(),
                                   parMetisGraph.getXAdj(),
                                   parMetisGraph.getAdjncy(),
                                   &(wgts[0]),
                                   NULL,
                                   &(vsize[0]),
                                   &wgtflag,
                                   &numflag,
                                   &ncon,
                                   &nparts,
                                   &(tpwgts[0]),
                                   &ubvec,
                                   &itr,
                                   options,
                                   &edgecut,
                                   &(part[0]),
                                   &tmpComm);
      }
      break;
      case REFINE_PART:
        ParMETIS_V3_RefineKway(parMetisMesh->getElementDist(),
                               parMetisGraph.getXAdj(),
                               parMetisGraph.getAdjncy(),
                               &(wgts[0]),
                               NULL,
                               &wgtflag,
                               &numflag,
                               &ncon,
                               &nparts,
                               &(tpwgts[0]),
                               &ubvec,
                               options,
                               &edgecut,
                               &(part[0]),
                               &tmpComm);

        break;
      default:
        ERROR_EXIT("unknown partitioning mode\n");
      }


      // === Distribute new partition data. ===

      return distributePartitioning(&(part[0]));
    }


    void ParMetisPartitioner::createPartitionMap(map<int, int>& pMap)
    {
      FUNCNAME("ParMetisPartitioner::createPartitionMap()");

      partitionMap.clear();

      // update ParMETIS mesh to new partitioning
      if (!parMetisMesh)
        parMetisMesh = new ParMetisMesh(mesh, mpiComm, elementInRank, mapLocalGlobal);

      int mpiRank = mpiComm->Get_rank();
      int mpiSize = mpiComm->Get_size();
      vector<int> nPartitionElements(mpiSize);
      int* elmdist = parMetisMesh->getElementDist();

      for (int i = 0; i < mpiSize; i++)
        nPartitionElements[i] = elmdist[i + 1] - elmdist[i];

      // === count number of elements ===
      int nElements = 0;
      int localElements = parMetisMesh->getNumElements();
      mpiComm->Allreduce(&localElements, &nElements, 1, MPI_INT, MPI_SUM);

      vector<int> partitionElements(nElements);

      // distribute partition elements
      mpiComm->Allgatherv(parMetisMesh->getAMDiSIndices(),
                          nPartitionElements[mpiRank],
                          MPI_INT,
                          &(partitionElements[0]),
                          &(nPartitionElements[0]),
                          elmdist,
                          MPI_INT);

      // fill partitionMap
      for (int i = 0; i < mpiSize; i++)
        for (int j = 0; j < nPartitionElements[i]; j++)
          partitionMap[partitionElements[elmdist[i] + j]] = i;

      pMap = partitionMap;
    }


    bool ParMetisPartitioner::distributePartitioning(int* part)
    {
      FUNCNAME("ParMetisPartitioner::distributePartitioning()");

      int mpiSize = mpiComm->Get_size();
      int mpiRank = mpiComm->Get_rank();
      int nElements = parMetisMesh->getNumElements();

      // nPartitionElements[i] is the number of elements for the i-th partition
      int* nPartitionElements = new int[mpiSize];
      for (int i = 0; i < mpiSize; i++)
        nPartitionElements[i] = 0;
      for (int i = 0; i < nElements; i++)
        nPartitionElements[part[i]]++;

      // collect number of partition elements from all ranks for this rank
      int* nRankElements = new int[mpiSize];
      mpiComm->Alltoall(nPartitionElements, 1, MPI_INT, nRankElements, 1, MPI_INT);


      // sum up partition elements over all ranks
      int* sumPartitionElements = new int[mpiSize];
      mpiComm->Allreduce(nPartitionElements, sumPartitionElements, mpiSize,
                         MPI_INT, MPI_SUM);

      // Test if there exists an empty partition
      bool emptyPartition = false;
      for (int i = 0; i < mpiSize; i++)
        emptyPartition |= (sumPartitionElements[i] == 0);

      if (emptyPartition)
        return false;

      // prepare distribution (fill partitionElements with AMDiS indices)
      int* bufferOffset = new int[mpiSize];
      bufferOffset[0] = 0;
      for (int i = 1; i < mpiSize; i++)
        bufferOffset[i] = bufferOffset[i - 1] + nPartitionElements[i - 1];

      int* partitionElements = new int[nElements];
      int** partitionPtr = new int* [mpiSize];

      for (int i = 0; i < mpiSize; i++)
        partitionPtr[i] = partitionElements + bufferOffset[i];

      sendElements.clear();
      for (int i = 0; i < nElements; i++)
      {
        int partition = part[i];
        int amdisIndex = parMetisMesh->getAMDiSIndex(i);

        if (partition != mpiRank)
          sendElements[partition].push_back(amdisIndex);

        *(partitionPtr[partition]) = amdisIndex;
        ++(partitionPtr[partition]);
      }

      // all to all: partition elements to rank elements
      int* rankElements = new int[sumPartitionElements[mpiRank]];
      int* recvBufferOffset = new int[mpiSize];
      recvBufferOffset[0] = 0;
      for (int i = 1; i < mpiSize; i++)
        recvBufferOffset[i] = recvBufferOffset[i - 1] + nRankElements[i - 1];

      mpiComm->Alltoallv(partitionElements,
                         nPartitionElements,
                         bufferOffset,
                         MPI_INT,
                         rankElements,
                         nRankElements,
                         recvBufferOffset,
                         MPI_INT);

      TEST_EXIT(elementInRank.size() != 0)("Should not happen!\n");
      for (map<int, bool>::iterator it = elementInRank.begin();
           it != elementInRank.end(); ++it)
        elementInRank[it->first] = false;

      // Create map which stores for each element index on macro level
      // if the element is in the partition of this rank.
      recvElements.clear();
      for (int i = 0; i < mpiSize; i++)
      {
        int* rankStart = rankElements + recvBufferOffset[i];
        int* rankEnd = rankStart + nRankElements[i];
        for (int* rankPtr = rankStart; rankPtr < rankEnd; ++rankPtr)
        {
          elementInRank[*rankPtr] = true;
          if (i != mpiRank)
            recvElements[i].push_back(*rankPtr);
        }
      }

      delete parMetisMesh;
      parMetisMesh = NULL;

      delete [] rankElements;
      delete [] nPartitionElements;
      delete [] nRankElements;
      delete [] sumPartitionElements;
      delete [] partitionElements;
      delete [] partitionPtr;
      delete [] bufferOffset;
      delete [] recvBufferOffset;

      return true;
    }

  }
}
