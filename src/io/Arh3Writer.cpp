#include "io/Arh3Writer.hpp"

#include "Mesh.hpp"
#include "MeshStructure.hpp"
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/StdMpi.hpp"
#endif

namespace AMDiS
{
  namespace io
  {

    namespace Arh3Writer
    {

      void writeFile(SystemVector* sysVec,
                    std::string filename,
                    bool writeParallel,
                    Cpsformat cps,
                    std::string dataformat)
      {
        std::vector<DOFVector<double>*> vecs;
        for (int i = 0; i < sysVec->getSize(); i++)
          vecs.push_back(sysVec->getDOFVector(i));
        detail::write(filename, NULL, vecs, writeParallel, cps, dataformat);
      }

      /// write the meshstructure and the dof values of DOFVectors in vec0
      /// the behavior is equal to writeFile(SystemVector* sysVec, string filename).
      void writeFile(DOFVector<double>* vec0,
                    std::string filename,
                    bool writeParallel,
                    Cpsformat cps,
                    std::string dataformat)
      {
        std::vector<DOFVector<double>*> vecs;
        vecs.push_back(vec0);
        detail::write(filename, NULL, vecs, writeParallel, cps, dataformat);
      }

      /// write the meshstructure and the dof values of DOFVectors in vecs
      /// the behavior is equal to writeFile(SystemVector* sysVec, string filename).
      void writeFile(std::vector<DOFVector<double>*> vecs,
                    std::string filename,
                    bool writeParallel,
                    Cpsformat cps,
                    std::string dataformat)
      {
        detail::write(filename, NULL, vecs, writeParallel, cps, dataformat);
      }

      /// write the meshstructure of the mesh to arh file.
      void writeFile(Mesh* mesh,
                    std::string filename,
                    bool writeParallel,
                    Cpsformat cps,
                    std::string dataformat)
      {
        std::vector<DOFVector<double>*> vecs;
        detail::write(filename, mesh, vecs, writeParallel, cps, dataformat);
      }

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
} // end namespace AMDiS
