#include "io/detail/Arh2Reader.hpp"

#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#ifdef HAVE_COMPRESSION
#include <boost/iostreams/filter/zlib.hpp>
#endif

#include "Mesh.hpp"
#include "MeshStructure.hpp"
#include "Traverse.hpp"
#include "DOFVector.hpp"
#include "SystemVector.hpp"
#include "Debug.hpp"
#include "io/Arh2Reader.hpp"
#include "MacroElement.hpp"

namespace AMDiS
{
  namespace io
  {

    using namespace std;

    namespace Arh2Reader
    {
      namespace detail
      {

        uint8_t firstRead(ifstream& file)
        {
          FUNCNAME("Arh2Reader::detail::firstRead");

          string typeId(4, ' ');
          uint8_t major = 0, minor = 0;

          file.read(const_cast<char*>(typeId.data()), 4);
          TEST_EXIT(typeId == "arh2")
          ("Cannot read Arh2 format: this file is not \"arh2\" format.\n");

          file.read(reinterpret_cast<char*>(&major), 1);
          TEST_EXIT(major == AMDiS::io::Arh2Reader::MAJOR)
          ("Cannot read Arh2 format: Arh2Reader major version is %i, the file major version is %i. \n",
           AMDiS::io::Arh2Reader::MAJOR, major);

          file.read(reinterpret_cast<char*>(&minor), 1);
          TEST_EXIT(minor <= AMDiS::io::Arh2Reader::MINOR)
          ("Cannot read Arh2 format: Arh2Reader minor version is %i is smaller than the file minor version %i.\n",
           AMDiS::io::Arh2Reader::MINOR, minor);
          return minor;
        }

        void setDofValues(int macroElIndex, Mesh* mesh,
                          vector<vector<double>>& values,
                          vector<DOFVector<double>*>& vecs,
                          vector<vector<int>>& feSpaces)
        {
          if(vecs.size() == 0)
            return;

          std::set<int> unique;
          std::pair<std::set<int>::iterator,bool> ret;
          std::vector<std::set<int>> visited;
          std::vector<int> dofIndex;
          int globalDof = 0;

          for(size_t i = 0; i < feSpaces.size(); i++)
          {
            visited.push_back(unique);
            dofIndex.push_back(0);
          }

          TraverseStack stack;
          ElInfo* elInfo = stack.traverseFirstOneMacro(mesh, macroElIndex, -1,
                           Mesh::CALL_EVERY_EL_PREORDER);
          while (elInfo)
          {
            Element* el = elInfo->getElement();
            if (el->isLeaf())
            {
              int valuePos = 0;
              for(size_t i = 0; i < feSpaces.size(); i++)
              {
                int num = -1;
                for(size_t j = 0; j < feSpaces[i].size(); j++)
                {
                  if((vecs.size() > (unsigned)feSpaces[i][j]) && vecs[feSpaces[i][j]])
                  {
                    num = feSpaces[i][j];
                    break;
                  }
                }
                if(num == -1)
                {
                  valuePos += feSpaces[i].size();
                  continue;
                }
                DOFAdmin* admin = vecs[num]->getFeSpace()->getAdmin();
                TEST_EXIT(admin)("the DOFAdmin of DOFVector is null, this should not happen.");

                int n0, nd, nd0;

                if ((nd = admin->getNumberOfDofs(VERTEX)))
                {
                  int vertices = mesh->getGeo(VERTEX);
                  nd0 = admin->getNumberOfPreDofs(VERTEX);
                  n0 = mesh->getNode(VERTEX);
                  for (int n = 0; n < vertices; n++)
                    for(int d = 0; d < nd; d++)
                    {
                      globalDof = elInfo->getElement()->getDof(n0 + n, nd0 + d);
                      ret = visited[i].insert(globalDof);
                      if(ret.second)
                      {
                        for(size_t j = 0 ; j < feSpaces[i].size(); j++)
                        {
                          if((vecs.size() > (unsigned)feSpaces[i][j]) && vecs[feSpaces[i][j]])
                          {
                            (*vecs[feSpaces[i][j]])[globalDof] = values[valuePos + j][dofIndex[i]];
                          }
                        }
                        dofIndex[i]++;
                      }
                    }
                }
                if (mesh->getDim() > 1)
                {
                  if ((nd = admin->getNumberOfDofs(EDGE)))
                  {
                    int edges = mesh->getGeo(EDGE);
                    nd0 = admin->getNumberOfPreDofs(EDGE);
                    n0 = mesh->getNode(EDGE);
                    for (int n = 0; n < edges; n++)
                      for(int d = 0; d < nd; d++)
                      {
                        globalDof = elInfo->getElement()->getDof(n0 + n, nd0 + d);
                        ret = visited[i].insert(globalDof);
                        if(ret.second)
                        {
                          for(size_t j = 0 ; j < feSpaces[i].size(); j++)
                          {
                            if((vecs.size() > (unsigned)feSpaces[i][j]) && vecs[feSpaces[i][j]])
                              (*vecs[feSpaces[i][j]])[globalDof] = values[valuePos + j][dofIndex[i]];
                          }
                          dofIndex[i]++;
                        }
                      }
                  }
                }
                if (mesh->getDim() == 3)
                {
                  if ((nd = admin->getNumberOfDofs(FACE)))
                  {
                    int faces = mesh->getGeo(FACE);
                    nd0 = admin->getNumberOfPreDofs(FACE);
                    n0 = mesh->getNode(FACE);
                    for (int n = 0; n < faces; n++)
                      for(int d = 0; d < nd; d++)
                      {
                        globalDof = elInfo->getElement()->getDof(n0 + n, nd0 + d);
                        ret = visited[i].insert(globalDof);
                        if(ret.second)
                        {
                          for(size_t j = 0 ; j < feSpaces[i].size(); j++)
                          {
                            if((vecs.size() > (unsigned)feSpaces[i][j]) && vecs[feSpaces[i][j]])
                              (*vecs[feSpaces[i][j]])[globalDof] = values[valuePos + j][dofIndex[i]];
                          }
                          dofIndex[i]++;
                        }
                      }
                  }
                }
                if ((nd = admin->getNumberOfDofs(CENTER)))
                {
                  nd0 = admin->getNumberOfPreDofs(CENTER);
                  n0 = mesh->getNode(CENTER);
                  for(int d = 0; d < nd; d++)
                  {
                    globalDof = elInfo->getElement()->getDof(n0, nd0 + d);
                    ret = visited[i].insert(globalDof);
                    if(ret.second)
                    {
                      for(size_t j = 0 ; j < feSpaces[i].size(); j++)
                      {
                        if((vecs.size() > (unsigned)feSpaces[i][j]) && vecs[feSpaces[i][j]])
                          (*vecs[feSpaces[i][j]])[globalDof] = values[valuePos + j][dofIndex[i]];
                      }
                      dofIndex[i]++;
                    }
                  }
                }
                valuePos += feSpaces[i].size();
              }
            }
            elInfo = stack.traverseNext(elInfo);
          }
        }

        void read(string filename,
                  Mesh* mesh,
                  vector<DOFVector<double>*> vecs,
                  bool byName)
        {
          FUNCNAME("Arh2Reader::detail::read()");

          // Get set of all macro elements in mesh.
          std::set<int> macroInMesh;
          for (std::deque<MacroElement*>::iterator it = mesh->getMacroElements().begin();
               it != mesh->getMacroElements().end(); ++it)
            macroInMesh.insert((*it)->getIndex());

          RefinementManager* refManager = NULL;
          switch (mesh->getDim())
          {
          case 1:
            refManager = new RefinementManager1d();
            break;
          case 2:
            refManager = new RefinementManager2d();
            break;
          case 3:
            refManager = new RefinementManager3d();
            break;
          default:
            ERROR_EXIT("Should not happen!\n");
          }
          ifstream file;
          file.open(filename.c_str(), ios::in | ios::binary);
          TEST_EXIT(file.is_open())
          ("Cannot open file %s\n", filename.c_str());

          file.seekg (0, file.end);
          int fileSize = file.tellg();
          file.seekg (0, file.beg);

          string cps = "null";
          uint32_t headerLen = 0;
          uint32_t nMacroElements = 0;
          uint32_t nValueVectors = 0;
          uint32_t nFeSpaces = 0;
          uint32_t dim = 0, dow = 0;

          // Read fixed header
          uint8_t minor = firstRead(file);

          file.read(reinterpret_cast<char*>(&headerLen), 4);
          file.read(reinterpret_cast<char*>(&dow), 4);
          file.read(reinterpret_cast<char*>(&dim), 4);
          file.read(reinterpret_cast<char*>(&nFeSpaces), 4);
          file.read(reinterpret_cast<char*>(&nValueVectors), 4);
          file.read(reinterpret_cast<char*>(&nMacroElements), 4);

          if(minor >= 1)
            file.read(const_cast<char*>(cps.data()), 4);

#ifdef HAVE_COMPRESSION
          TEST_EXIT(cps == "null" || cps == "zlib")
          ("Cannot read Arh2 file. Currently only support zlib compression.\n");
#else
          TEST_EXIT(cps == "null")
          ("HAVE_COMPRESSION OFF. Cannot read compressed Arh2 file.\n");
#endif

          TEST_EXIT(dow == (unsigned)mesh->getGeo(WORLD))
          ("Dow is not match!\n");
          TEST_EXIT(dim == (unsigned)mesh->getDim())
          ("Dim is not match!\n");
          TEST_EXIT(nValueVectors >= vecs.size())
          ("File %s has %d vector(s), which is less than the number of DOFVectors %i in vecs!\n",
           filename.c_str(), nValueVectors, vecs.size());

          vector<int> vecsNameLen(0);
          vector<string> vecsName(0);
          vector<int> vecsFeSpaceNum(0);
          vector<int> macroElIndex(nMacroElements);
          vector<int> macroElSize(nMacroElements);

          vector<vector<int>> feSpaceDOFs(0);
          vector<int> perDOFs(4, 0);
          vector<vector<int>> sortedFeSpaces(nFeSpaces);

          // Read macro table
          for(unsigned int i = 0; i < nMacroElements; i++)
          {
            file.read(reinterpret_cast<char*>(&macroElIndex[i]), 4);
            file.read(reinterpret_cast<char*>(&macroElSize[i]), 4);
          }
          // Read feSpace table
          for(unsigned int i = 0; i < nFeSpaces; i++)
          {
            for(int j = 0; j < 4; j++)
            {
              file.read(reinterpret_cast<char*>(&perDOFs[j]), 4);
            }
            feSpaceDOFs.push_back(perDOFs);
          }
          // Read value table
          for(unsigned int i = 0; i < nValueVectors; i++)
          {
            string tmpString("");
            uint32_t tmpInt = 0;
            file.read(reinterpret_cast<char*>(&tmpInt), 4);
            vecsNameLen.push_back(tmpInt);
            tmpString.resize(tmpInt, ' ');
            file.read(const_cast<char*>(tmpString.data()), tmpInt); //
            vecsName.push_back(tmpString);
            file.read(reinterpret_cast<char*>(&tmpInt), 4);
            sortedFeSpaces[tmpInt].push_back(i);
            vecsFeSpaceNum.push_back(tmpInt);
          }
          // Adjust and check vecs
          if(byName)
          {
            if(!vecs.empty())
            {
              vector<DOFVector<double>*> tmpVecs = vecs;
              vecs.clear();
              vecs.resize(nValueVectors, NULL);

              for(size_t k = 0; k < tmpVecs.size(); k++)
              {
                if(!tmpVecs[k])
                  break;

                unsigned int i;
                TEST_EXIT(tmpVecs[k]->getFeSpace()->getAdmin())
                ("Vecs number %i has no DOFAdmin. Should not happen.\n", k);
                DimVec<int>* nDOF = tmpVecs[k]->getFeSpace()->getBasisFcts()->getNumberOfDofs();

                for(i = 0; i < nValueVectors; i++)
                {
                  if(tmpVecs[k]->getName() != vecsName[i])
                  {
                    continue;
                  }
                  int j;
                  for(j = 0; j < nDOF->getSize(); j++)
                  {
                    if((*nDOF)[j] != feSpaceDOFs[vecsFeSpaceNum[i]][j])
                    {
                      break;
                    }
                  }
                  if(j == nDOF->getSize())
                  {
                    vecs[i] = tmpVecs[k];
                    break;
                  }
                }
                TEST_EXIT(i < nValueVectors)
                ("NO DOFVector with the same name and feSpace is found in the file.\n");
              }
            }
          }
          else
          {
            for(size_t i = 0; i < vecs.size(); i++)
            {
              if(vecs[i])
              {
                TEST_EXIT(vecs[i]->getFeSpace()->getAdmin())
                ("Vecs number %i has no DOFAdmin. Should not happen.\n", i);

                DimVec<int>* nDOF = vecs[i]->getFeSpace()->getBasisFcts()->getNumberOfDofs();
                for(int j = 0; j < nDOF->getSize(); j++)
                {
                  TEST_EXIT((*nDOF)[j] == feSpaceDOFs[vecsFeSpaceNum[i]][j])
                  ("The fespace of vec number %i is not equal to the correspond fespace.\n", i+1);
                }
              }
            }
          }
          // Read data: meshstructure and dof values
          for (unsigned int i = 0; i < nMacroElements; i++)
          {
            stringstream dataStream(ios::out | ios::in | ios::binary);
            int size = 0;
            if(minor >= 1)
              size = macroElSize[i];
            else
              size = (i != nMacroElements - 1) ? macroElSize[i + 1] - macroElSize[i] : fileSize - macroElSize[i];
            char* buffer = new char[size];
            file.read(buffer, size);
            dataStream.write(buffer, size);
            delete[] buffer;
#ifdef HAVE_COMPRESSION
            if(cps == "zlib")
            {
              stringstream tmp(ios::out | ios::in);
              boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
              in.push(boost::iostreams::zlib_decompressor());
              in.push(dataStream);
              boost::iostreams::copy(in, tmp);
              dataStream.str(tmp.str());
            }
#endif
            uint32_t nStructureCodes = 0;
            uint32_t codeSize = 0;

            dataStream.read(reinterpret_cast<char*>(&nStructureCodes), 4);
            dataStream.read(reinterpret_cast<char*>(&codeSize), 4);

            vector<uint64_t> structureCode(nStructureCodes);
            dataStream.read(reinterpret_cast<char*>(&(structureCode[0])), 8 * nStructureCodes);

            MeshStructure elementStructure;
            elementStructure.init(structureCode, codeSize);
            if (macroInMesh.count(macroElIndex[i]))
              elementStructure.fitMeshToStructure(mesh, refManager, false, macroElIndex[i]);

            if (nValueVectors > 0)
            {
              vector<vector<double>> values(nValueVectors);
              int valuePos = 0;
              for(unsigned int j = 0; j < nFeSpaces; j++)
              {
                uint32_t nValuesPerVector = 0;
                dataStream.read(reinterpret_cast<char*>(&nValuesPerVector), 4);

                for(size_t k = 0; k < sortedFeSpaces[j].size(); k++)
                {
                  values[valuePos + k].resize(nValuesPerVector);
                  dataStream.read(reinterpret_cast<char*>(&(values[valuePos + k][0])), 8 * nValuesPerVector);
                }
                valuePos += sortedFeSpaces[j].size();
              }
              if(!vecs.empty())
              {
                if(macroInMesh.count(macroElIndex[i]))
                {
                  setDofValues(macroElIndex[i], mesh, values, vecs, sortedFeSpaces);
                }
              }
            }
          }
          file.close();
          delete refManager;
        }

        void readFile(string filename, Mesh* mesh,
                      vector<DOFVector<double>*> vecs,
                      bool writeParallel,
                      int nProcs,
                      bool byName)
        {
          FUNCNAME("Arh2Reader::detail::readFile()");
          //The last element in vecs must NOT be NULL
          std::set<string> nameSet;
          pair<std::set<string>::iterator,bool> ret;

          while(!vecs.empty())
          {
            if(vecs.back() == NULL)
              vecs.pop_back();
            else
              break;
          }
          // This is the precondition
          for(size_t i = 0; i < vecs.size(); i++)
          {
            if(vecs[i])
            {
              if(!mesh)
                mesh = vecs[i]->getFeSpace()->getMesh();
              else
                TEST_EXIT(mesh == vecs[i]->getFeSpace()->getMesh())
                ("The mesh of the DOFVectors should be the same for Reader because in one file there is only one mesh.\n");
              ret = nameSet.insert(vecs[i]->getName());
              TEST_EXIT(ret.second)("DOFVectors in vecs cannot have idential name. Please check.\n");
            }
          }
          if(!mesh)
          {
            WARNING("You haven't specified the target, no mesh or DOFVectors is given.\n");
            return;
          }
          if (writeParallel)
          {
            int sPos = filename.find(".arh");
            TEST_EXIT(sPos >= 0)("Failed to find file postfix!\n");
            string name = filename.substr(0, sPos);

            if (nProcs == -1)
            {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
              string procFilename = name + "-p" + std::to_string(MPI::COMM_WORLD.Get_rank()) + "-.arh";
              read(procFilename, mesh, vecs, byName);
              MSG("ARH file read from: %s\n", procFilename.c_str());
#else
              ERROR_EXIT("Reading parallel ARH files in sequential computations requires to specify the number of nodes on which the ARH file was created!\n");
#endif
            }
            else
            {

              // check if there are no more or less files as nProcs
              int n = 0;
              for (; n < nProcs + 1; n++)
              {
                string fn = name + "-p" + std::to_string(n) + "-.arh";
                if(!boost::filesystem::exists(fn)) break;
              }
              TEST_EXIT(n == nProcs)
              ("Number of arh files doesn't match number of processors \n");

              //
              // Some test should be checked. This is because of the variation of the
              // number of macro elements per processor:
              //
              //    There should be at least 10 macro Elements per processor, therefore:
              //        nMacroElements * 2^gr >= nProcs * 10
              //          =>  gr = log_2(nProcs * 10 / nMacroElements)

              int allMacros = mesh->getNumberOfMacros();
              int allMacrosFromProcFiles = 0;

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
              int sendValue = static_cast<int>(mesh->getNumberOfMacros());
              MPI::COMM_WORLD.Allreduce(&sendValue, &allMacros, 1, MPI_INT, MPI_SUM);

              if(MPI::COMM_WORLD.Get_rank() == 0)
              {
#endif
                for(int i = 0; i < nProcs; i++)
                {
                  allMacrosFromProcFiles += readNumOfMacrosFromSgArh(filename, i);
                }

                TEST_EXIT(allMacros == allMacrosFromProcFiles)
                ("Number of macro elements in parallel ARH files doesn't match to the current runtime.\n");
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
              }
#endif
              for (int i = 0; i < nProcs; i++)
              {
                string procFilename = name + "-p" + std::to_string(i) + "-.arh";
                read(procFilename, mesh, vecs, byName);
                MSG("ARH file read from: %s\n", procFilename.c_str());
              }
            }
          }
          else
          {
            read(filename, mesh, vecs, byName);
          }
          MSG("ARH file read from: %s\n", filename.c_str());
        }

        void readMetaFromSgArh(std::string filename, int nProc,
                               std::vector<std::set<std::pair<int, int>>>& data)
        {
          ifstream file;
          file.open(filename.c_str(), ios::in | ios::binary);
          TEST_EXIT(file.is_open())
          ("Cannot open file %s\n", filename.c_str());

          file.seekg (0, file.end);
          int fileSize = file.tellg();
          file.seekg (0, file.beg);

          uint32_t dow = 0, dim = 0;
          uint32_t headerLen = 0;
          uint32_t nMacroElements = 0;
          uint32_t nValueVectors = 0;
          uint32_t nFeSpaces = 0;
          string cps = "null";

          uint8_t minor = firstRead(file);

          file.read(reinterpret_cast<char*>(&headerLen), 4);
          file.read(reinterpret_cast<char*>(&dow), 4);
          file.read(reinterpret_cast<char*>(&dim), 4);
          file.read(reinterpret_cast<char*>(&nFeSpaces), 4);
          file.read(reinterpret_cast<char*>(&nValueVectors), 4);
          file.read(reinterpret_cast<char*>(&nMacroElements), 4);

          if(minor > 0)
            file.read(const_cast<char*>(cps.data()), 4);

#ifdef HAVE_COMPRESSION
          TEST_EXIT(cps == "null" || cps == "zlib")
          ("Cannot read Arh2 file. Currently only support zlib compression.\n");
#else
          TEST_EXIT(cps == "null")
          ("HAVE_COMPRESSION OFF. Cannot read compressed Arh2 file.\n");
#endif

          vector<int> macroElIndex(nMacroElements);
          vector<int> macroElSize(nMacroElements);
          vector<vector<int>> sortedFeSpaces(nFeSpaces);

          // Read macro table
          for(uint32_t i = 0; i < nMacroElements; i++)
          {
            file.read(reinterpret_cast<char*>(&macroElIndex[i]), 4);
            file.read(reinterpret_cast<char*>(&macroElSize[i]), 4);
          }

          // Read feSpace table
          file.seekg(nFeSpaces * 4 * 4, ios_base::cur);

          // Read value table
          for(uint32_t i = 0; i < nValueVectors; i++)
          {
            string tmpString("");
            uint32_t tmpInt = 0;
            file.read(reinterpret_cast<char*>(&tmpInt), 4);
            tmpString.resize(tmpInt, ' ');
            file.read(const_cast<char*>(tmpString.data()), tmpInt); //
            file.read(reinterpret_cast<char*>(&tmpInt), 4);
            sortedFeSpaces[tmpInt].push_back(i);
          }

          for (uint32_t i = 0; i < nMacroElements; i++)
          {
            stringstream dataStream(ios::out | ios::in | ios::binary);
            int size = 0;

            if(minor > 0)
              size = macroElSize[i];
            else
              size = (i != nMacroElements - 1) ? macroElSize[i + 1] - macroElSize[i] : fileSize - macroElSize[i];

            char* buffer = new char[size];
            file.read(buffer, size);
            dataStream.write(buffer, size);
            delete[] buffer;

#ifdef HAVE_COMPRESSION
            if(cps == "zlib")
            {
              stringstream tmp(ios::out | ios::in);
              boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
              in.push(boost::iostreams::zlib_decompressor());
              in.push(dataStream);
              boost::iostreams::copy(in, tmp);
              dataStream.str(tmp.str());
            }
#endif
            uint32_t nStructureCodes = 0;
            uint32_t codeSize = 0;

            dataStream.read(reinterpret_cast<char*>(&nStructureCodes), 4);
            dataStream.read(reinterpret_cast<char*>(&codeSize), 4);

            data[nProc].insert(make_pair(macroElIndex[i], codeSize));

            // We only need codeSize for each macro element data block, so skip all of the unnecessary.
            dataStream.seekg(8 * nStructureCodes, ios_base::cur);

            if (nValueVectors > 0)
            {

              for(uint32_t j = 0; j < nFeSpaces; j++)
              {
                uint32_t nValuesPerVector = 0;
                dataStream.read(reinterpret_cast<char*>(&nValuesPerVector), 4);
                dataStream.seekg(sortedFeSpaces[j].size() * 8 * nValuesPerVector, ios_base::cur);
              }
            }
          }
        } // end readMetaFromSgArh


        int readNumOfMacrosFromSgArh(std::string filename, int nProc)
        {
          FUNCNAME("Arh2Reader::readHeaderSize");

          int sPos = filename.find(".arh");
          TEST_EXIT(sPos >= 0)("Failed to find file postfix!\n");

          if (nProc >= 0)
          {
            string name = filename.substr(0, sPos);
            filename = name + "-p" + std::to_string(nProc) + "-.arh";
          }

          ifstream file;

          file.open(filename.c_str(), ios::in | ios::binary);
          TEST_EXIT(file.is_open())
          ("Cannot open file %s\n", filename.c_str());

          int nMacroElements = 0;
          detail::firstRead(file);
          file.seekg(20, ios_base::cur);
          file.read(reinterpret_cast<char*>(&nMacroElements), 4);
          file.close();

          return nMacroElements;
        }


      } // end namespace detail
    } // end namespace Arh2Reader
  }
} // end namespace io, AMDiS