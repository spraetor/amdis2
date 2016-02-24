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

#include "Arh3Writer.h"
#include "Mesh.h"
#include "MeshStructure.h"
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/StdMpi.hpp"
#endif

namespace AMDiS
{
  namespace io
  {

    namespace Arh3Writer
    {

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
      void writeMetaData(Mesh* mesh, std::string metaFilename)
      {
        FUNCNAME("Arh3Writer::writeMetaData()");

        using namespace std;
        using namespace AMDiS::Parallel;

        int mpiSize = MPI::COMM_WORLD.Get_size();
        vector<std::set<pair<int, int>>> overallData;
        std::set<pair<int, int>> data;

        // Calculate local data

        MeshStructure elementStructure;
        int macroElIndex = -1;

        TraverseStack stack;
        ElInfo* elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
        while (elInfo)
        {
          if (elInfo->getLevel() == 0)
          {
            if (macroElIndex != -1)
            {
              elementStructure.commit();

              data.insert(make_pair(macroElIndex, elementStructure.getNumElements()));
            }
            elementStructure.clear();

            macroElIndex = elInfo->getElement()->getIndex();
          }
          elementStructure.insertElement(elInfo->getElement()->isLeaf());
          elInfo = stack.traverseNext(elInfo);
        }

        TEST_EXIT_DBG(macroElIndex != -1)("Should not happen!\n");
        elementStructure.commit();
        data.insert(make_pair(macroElIndex, elementStructure.getNumElements()));

        // Collect data from other processors

        StdMpi<std::set<pair<int, int>>> stdMpi(MPI::COMM_WORLD);

        if(MPI::COMM_WORLD.Get_rank() == 0)
        {
          for(int rank = 1; rank < mpiSize; rank++)
            stdMpi.recv(rank);
        }
        else
        {
          stdMpi.send(0, data);
        }

        stdMpi.startCommunication();

        if(MPI::COMM_WORLD.Get_rank() == 0)
        {
          overallData.push_back(data);

          for(int rank = 1; rank < mpiSize; rank++)
          {
            std::set<pair<int, int>>& recvData = stdMpi.getRecvData(rank);
            overallData.push_back(recvData);
          }

          // Write to meta file

          ofstream file;
          file.open(metaFilename.c_str());
          file << "METAARH\n";
          file << "p" << "\n";	//edit 150420: just write something
          file << mpiSize << "\n";
          for (int i = 0; i < mpiSize; i++)
          {
            file << i << " " << overallData[i].size() << "\n";
            for (std::set<pair<int, int>>::iterator it = overallData[i].begin(); it != overallData[i].end(); ++it)
              file << it->first << " " << it->second << "\n";
          }
          file.close();
        }
      }
#endif

    }
  }
}
