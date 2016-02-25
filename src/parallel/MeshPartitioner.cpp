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


#include "parallel/MeshPartitioner.hpp"
#include "io/ArhReader.hpp"
#include "io/Arh2Reader.hpp"
#include "Mesh.h"
#include "Traverse.h"
#include "Serializer.h"

using namespace std;

namespace AMDiS
{
  namespace Parallel
  {

    bool MeshPartitioner::createInitialPartitioning()
    {
      FUNCNAME("MeshPartitioner::createInitialPartitioning()");

      int mpiRank = mpiComm->Get_rank();
      int mpiSize = mpiComm->Get_size();
      int nLeaves = mesh->getNumberOfLeaves();
      int elPerRank = nLeaves / mpiSize;
      bool useInitialPartitioning = false;

      // === Check for reading ARH meta file to create initial partitioning. ===

      map<int, int> mapElInRank;
      map<int, int> arhElCodeSize;

      string partitioningFile = "";

      Parameters::get(initFileStr + "->initial partitioning file",
                      partitioningFile);
      if (partitioningFile != "")
      {
        MSG("Read initial partitioning file: %s\n", partitioningFile.c_str());

        ifstream file;
        file.open(partitioningFile.c_str());
        TEST_EXIT(file.is_open())("Should not happen!\n");

        int nElements = 0;
        file >> nElements;
        for (int i = 0; i < nElements; i++)
        {
          int rank = -1;
          file >> rank;
          mapElInRank[i] = rank;
        }
        file.close();

        useInitialPartitioning = true;
      }
      else
      {
        string arhMetaFile = "";
        Parameters::get(initFileStr + "->read meta arh",
                        arhMetaFile);
        bool partitioningMetaArhBased = (arhMetaFile != "");

        string arhFile = "";
        Parameters::get(initFileStr + "->read arh",
                        arhFile);
        bool partitioningArhBased = (arhFile != "");

        TEST_EXIT(!(partitioningMetaArhBased && partitioningArhBased))
        ("decide for meta arh or arh files for initialise partition!\n");

        if (partitioningMetaArhBased)
        {
          MSG("Read Meta-Arh partitioning file: %s\n", arhMetaFile.c_str());
          int nProc = io::Arh2Reader::readMetaData(arhMetaFile, mapElInRank, arhElCodeSize);
          if (nProc != mpiSize)
            useInitialPartitioning = false;
          else
            useInitialPartitioning = true;
        }
        if (partitioningArhBased)
        {
          MSG("Read Arh partitioning files: %s\n", arhFile.c_str());
          int nProc = io::Arh2Reader::readMetaFromArh(arhFile, mapElInRank, arhElCodeSize);
          if (nProc != mpiSize)
            useInitialPartitioning = false;
          else
            useInitialPartitioning = true;
        }
      }


      // === Create initial partitioning of the AMDiS mesh. ===

      elementInRank.clear();

      // Is used in box partitioning mode only.
      map<DofEdge, std::set<int>> vertexElements;

      TraverseStack stack;
      ElInfo* elInfo =
        stack.traverseFirst(mesh, 0,
                            Mesh::CALL_EL_LEVEL | Mesh::FILL_NEIGH | Mesh::FILL_BOUND);
      while (elInfo)
      {
        Element* el = elInfo->getElement();
        int elIndex = el->getIndex();

        // === Store for all macro elements the interior neighbours (thus, no ===
        // === periodic neighbours) in the map elNeighbours.                  ===

        for (int i = 0; i < mesh->getGeo(NEIGH); i++)
          if (elInfo->getNeighbour(i) && elInfo->getBoundary(i) == INTERIOR)
            elNeighbours[elIndex].push_back(elInfo->getNeighbour(i)->getIndex());

        // === Create initial partitioning. ===

        if (!boxPartitioning)
        {
          // In standard mode assign to each macro element an arbitrary but unique
          // rank number.

          int elInRank = 0;
          if (mapElInRank.empty())
            elInRank = std::min(elIndex / elPerRank, mpiSize - 1);
          else
            elInRank = mapElInRank[elIndex];

          elementInRank[elIndex] = (elInRank == mpiRank);
          partitionMap[elIndex] = elInRank;
        }
        else
        {
          // In box partitioning mode, we do the assignment of boxes to ranks later.

          vertexElements[el->getEdge(0)].insert(elIndex);
        }

        elInfo = stack.traverseNext(elInfo);
      }


      // === Do initial partitioning in box partitioning mode. ===

      if (boxPartitioning)
      {
        TEST_EXIT(mesh->getDim() == 3)("Box partitioning only implemented for 3D!\n");


        // === Create boxes: all elements sharing the same 0-edge are assumed ===
        // === to be in the same box.                                         ===

        int boxCounter = 0;
        for (map<DofEdge, std::set<int>>::iterator it = vertexElements.begin();
             it != vertexElements.end(); ++it)
        {
          TEST_EXIT_DBG(it->second.size() == 6)("Should not happen!\n");

          boxSplitting[boxCounter] = it->second;

          for (std::set<int>::iterator elIt = it->second.begin();
               elIt != it->second.end(); ++elIt)
            elInBox[*elIt] = boxCounter;

          boxCounter++;
        }


        // === Calculate box neighbourhood relation. ===

        for (map<int, int>::iterator it = elInBox.begin(); it != elInBox.end(); ++it)
        {
          int elBoxNo = it->second;

          for (vector<int>::iterator neighIt = elNeighbours[it->first].begin();
               neighIt != elNeighbours[it->first].end(); ++neighIt)
          {
            int neighBoxNo = elInBox[*neighIt];

            if (elBoxNo != neighBoxNo)
              boxNeighbours[elBoxNo].insert(neighBoxNo);
          }
        }

        MSG("Box partitioning with %d boxes enabled!\n", boxCounter);


        /// === Create initial partitioning of the boxes. ====

        int boxPerRank = boxCounter / mpiSize;

        for (map<int, std::set<int>>::iterator it = boxSplitting.begin();
             it != boxSplitting.end(); ++it)
        {
          int boxInRank = std::min(it->first / boxPerRank, mpiSize - 1);

          // check if data compatible to boxpartition mode.
          // That is, every element of a box belongs to the same processor
          if(!mapElInRank.empty())
          {
            std::set<int>::iterator elIt = it->second.begin();
            boxInRank=mapElInRank[*elIt];
            for ( ; elIt != it->second.end(); ++elIt)
            {
              TEST_EXIT(mapElInRank[*elIt]== boxInRank)
              ("initial partion read from arh-file is not compatible to boxpartitioning. \n ");
            }
          }

          for (std::set<int>::iterator elIt = it->second.begin();
               elIt != it->second.end(); ++elIt)
          {
            elementInRank[*elIt] = (boxInRank == mpiRank);
            partitionMap[*elIt] = boxInRank;
          }
        }
      }

      return useInitialPartitioning;
    }


    void MeshPartitioner::serialize(std::ostream& out)
    {
      SerUtil::serialize(out, elementInRank);
      SerUtil::serialize(out, boxPartitioning);
      SerUtil::serialize(out, boxSplitting);
      SerUtil::serialize(out, boxNeighbours);
      SerUtil::serialize(out, elInBox);
    }


    void MeshPartitioner::deserialize(std::istream& in)
    {
      SerUtil::deserialize(in, elementInRank);
      SerUtil::deserialize(in, boxPartitioning);
      SerUtil::deserialize(in, boxSplitting);
      SerUtil::deserialize(in, boxNeighbours);
      SerUtil::deserialize(in, elInBox);
    }

  }
}
