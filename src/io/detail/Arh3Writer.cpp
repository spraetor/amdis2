#include "io/detail/Arh3Writer.hpp"

#include <fstream>
#include <stdint.h>
#include <iostream>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#ifdef HAVE_COMPRESSION
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#endif

#include "Mesh.hpp"
#include "MeshStructure.hpp"
#include "Traverse.hpp"
#include "DOFVector.hpp"
#include "io/Arh3Reader.hpp"
#include "MacroElement.hpp"

namespace AMDiS
{
  namespace io
  {

    using namespace std;

    namespace Arh3Writer
    {
      namespace detail
      {
        void write(string filename,
                   DOFVector<double>* vec0,
                   DOFVector<double>* vec1,
                   DOFVector<double>* vec2,
                   bool writeParallel,
                   Cpsformat cps,
                   string dataformat,
                   string filenameType)
        {
          vector<DOFVector<double>*> vecs(0);
          if (vec0 != NULL)
            vecs.push_back(vec0);
          if (vec1 != NULL)
            vecs.push_back(vec1);
          if (vec2 != NULL)
            vecs.push_back(vec2);

          write(filename, NULL, vecs, writeParallel, cps, dataformat, filenameType);
        }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
        void writeParallelFile(string filename, Mesh* mesh, string filenameType)
        {
          TEST_EXIT(filenameType == "cont")("Unsupported filename type.\n");

          ofstream file;
          file.open(filename.c_str(), ios::out | ios::binary | ios::trunc);

          string typeId = "parh", macroFilename = "";
          Parameters::get(mesh->getName() + "->macro file name", macroFilename);
          uint8_t major = 1;
          uint8_t minor = 0;
          uint32_t nFiles = MPI::COMM_WORLD.Get_size();
          uint32_t macroFile_nl = macroFilename.length();
          map<int, int> partitionMap =
            Parallel::MeshDistributor::globalMeshDistributor->getPartitionMap();
          uint32_t nMacros =partitionMap.size();

          file.write(typeId.c_str(), 4);
          file.write(reinterpret_cast<char*>(&major), 1);
          file.write(reinterpret_cast<char*>(&minor), 1);
          file.write(reinterpret_cast<char*>(&nFiles), 4);
          file.write(filenameType.c_str(), 4);
          file.write(reinterpret_cast<char*>(&nMacros), 4);
          file.write(reinterpret_cast<char*>(&macroFile_nl), 4);
          file.write(macroFilename.c_str(), macroFile_nl);

          map<int, int>::const_iterator it = partitionMap.begin();
          uint32_t rank = 0;
          for (; it != partitionMap.end(); it++)
          {
            rank = it->second;
            file.write(reinterpret_cast<char*>(&rank), 4);
          }
        }
#endif

        void write(string filename,
                   Mesh* mesh,
                   vector<DOFVector<double>*> vecs,
                   bool writeParallel,
                   Cpsformat cps,
                   string dataformat,
                   string /*filenameType*/) // NOTE: param filenametype only in parallel-mode
        {
          FUNCNAME("Arh3Writer::detail::write()");

          if (!mesh && vecs.empty())
          {
            WARNING("There is nothing to be writen.\n");
            return;
          }

          map<string,Valformat>::const_iterator it = dataformatMap.find(dataformat);
          TEST_EXIT(it != dataformatMap.end())("Wrong data format.\n");

          std::set<string> nameSet;
          pair<std::set<string>::iterator,bool> ret;

          for(size_t i = 0; i < vecs.size(); i++)
          {
            TEST_EXIT(vecs[i] != NULL)("Vecs[%i] is NULL. Please check.\n", i);
            ret = nameSet.insert(vecs[i]->getName());
            TEST_EXIT(ret.second)("DOFVectors in vecs cannot have idential name. Please check.\n");
          }

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
          if (writeParallel)
          {
            int sPos = filename.find(".arh");
            TEST_EXIT(sPos >= 0)("Failed to find file postfix!\n");
            string name = filename.substr(0, sPos);
            std::vector<int> macroIdx;

            Mesh* mesh_ = mesh ? mesh : vecs[0]->getFeSpace()->getMesh();

            if (MPI::COMM_WORLD.Get_rank() == 0)
            {
              writeParallelFile(name + ".parh", mesh_, filenameType);
            }

            TEST_EXIT(filenameType == "cont")("Only filename type \"cont\".\n");
            filename = name + "-p" + std::to_string(MPI::COMM_WORLD.Get_rank()) + "-.arh";
          }
#endif
          //if mesh exists, the meshes in vecs should be the same.
          if(mesh)
          {
            for(size_t i = 0; i < vecs.size(); i++)
            {
              TEST_EXIT(mesh == vecs[i]->getFeSpace()->getMesh())
              ("The mesh of DOFVector %i in vecs is not equal to the second parameter.\n", i);
            }
            writeAux(filename, mesh, vecs, writeParallel, cps, dataformat);
          }
          //multiple meshes are allowed here.
          else
          {
            vector<bool> visited(vecs.size(), false);
            vector<DOFVector<double>*> splitedVecs(0);
            bool moreMesh = false;
            Mesh* tmpMesh = vecs[0]->getFeSpace()->getMesh();
            for(size_t i = 1; i < vecs.size(); i++)
            {
              if(vecs[i]->getFeSpace()->getMesh() != tmpMesh)
              {
                moreMesh = true;
                break;
              }
            }
            for(size_t i = 0; i < vecs.size(); i++)
            {
              if(!visited[i])
              {
                splitedVecs.clear();
                splitedVecs.push_back(vecs[i]);
                visited[i] = true;
                tmpMesh = vecs[i]->getFeSpace()->getMesh();
                for(size_t j = i + 1; j < vecs.size(); j++)
                {
                  if(vecs[j]->getFeSpace()->getMesh() == tmpMesh)
                  {
                    splitedVecs.push_back(vecs[j]);
                    visited[j] = true;
                  }
                }
                string newfilename = filename;
                if(moreMesh)
                {
                  if(filename.length() > 4 && filename.substr(filename.length()- 4, filename.length()) == ".arh")
                    newfilename = filename.substr(0, filename.length() - 4) +
                                  "." + tmpMesh->getName() +
                                  filename.substr(filename.length()-4 , filename.length());
                  else
                    newfilename = filename + "." + tmpMesh->getName() + ".arh";

                }
                writeAux(newfilename, splitedVecs[0]->getFeSpace()->getMesh(), splitedVecs, writeParallel, cps, dataformat);
              }
            }
          }
        }

        int writeHeader(ofstream& file,
                        Mesh* mesh,
                        vector<DOFVector<double>*> vecs,
                        map<const FiniteElemSpace*, vector<int>>& feSpaces,
                        Cpsformat cps,
                        string dataformat)
        {
          FUNCNAME("Arh3Writer::detail::writeHeader()");
          // 	int nbits = boost::lexical_cast<int>(dataformat.substr(2, 2));
          TEST_EXIT(file.is_open())("the file is not open. should not happen.\n");

          map<const FiniteElemSpace*, string> AFEDfileName;

          uint32_t valueNamesLen = 0, fileNamesLen = 0;
          for (size_t i = 0; i < vecs.size(); i++)
            valueNamesLen += vecs[i]->getName().length();

          map<const FiniteElemSpace*, vector<int>>::iterator feSpaceIt;
          for (feSpaceIt = feSpaces.begin(); feSpaceIt != feSpaces.end(); feSpaceIt++)
            AFEDfileName.insert(make_pair(feSpaceIt->first, string()));

          uint32_t nValueVectors = vecs.size();
          uint32_t nFeSpaces = feSpaces.size();
          uint32_t nMacroElements = mesh->getNumberOfMacros();

          uint32_t nMacro = 0;
          TraverseStack st;
          ElInfo* el = st.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL);
          while (el)
          {
            nMacro++;
            el = st.traverseNext(el);
          }

          uint32_t dow = mesh->getGeo(WORLD);
          uint32_t dim = mesh->getDim();
          uint32_t headerLen = 34 +                         //fixed part of header
                               nMacroElements * 12 + 12 +   //macroElemnts table
                               fileNamesLen +		  //feSpaces table
                               nFeSpaces * 20 + 	          //feSpaces table
                               valueNamesLen + 	          //value vector table
                               nValueVectors * 12;           //also value vector table
          string typeId = "sarh";
#ifndef HAVE_COMPRESSION
          cps = NONE;
#endif
          uint8_t* major = const_cast<uint8_t*>(&(AMDiS::io::Arh3Reader::MAJOR));
          uint8_t* minor = const_cast<uint8_t*>(&(AMDiS::io::Arh3Reader::MINOR));
          int cpsflag = static_cast<int>(cps);
          uint32_t minus1 = -1;

          //fixed header
          file.write(typeId.c_str(), 4);
          file.write(reinterpret_cast<char*>(major), 1);
          file.write(reinterpret_cast<char*>(minor), 1);
          file.write(reinterpret_cast<char*>(&headerLen), 4);
          file.write(reinterpret_cast<char*>(&dow), 4);
          file.write(reinterpret_cast<char*>(&dim), 4);
          file.write(reinterpret_cast<char*>(&nFeSpaces), 4);
          file.write(reinterpret_cast<char*>(&nValueVectors), 4);
          file.write(reinterpret_cast<char*>(&nMacroElements), 4);
          file.write(reinterpret_cast<char*>(&cpsflag), 4);
          //macro table
          deque<MacroElement*>::const_iterator macroIter = mesh->firstMacroElement();

          while(macroIter != mesh->endOfMacroElements())
          {
            uint32_t macroIndex = (*macroIter)->getIndex();
            file.write(reinterpret_cast<char*>(&macroIndex), 4);
            file.write(reinterpret_cast<char*>(&minus1), 4);
            file.write(reinterpret_cast<char*>(&minus1), 4);
            macroIter++;
          }
          file.write(reinterpret_cast<char*>(&minus1), 4);
          file.write(reinterpret_cast<char*>(&minus1), 4);
          file.write(reinterpret_cast<char*>(&minus1), 4);

          vector<int> feSpaceNumOfVecs(vecs.size());
          uint32_t posDOFs = 0, nameStrLen = 0;
          string nameStr("");
          size_t i = 0;

          //feSpace table
          for(feSpaceIt = feSpaces.begin(); feSpaceIt != feSpaces.end(); feSpaceIt++, i++)
          {
            nameStr = AFEDfileName[feSpaceIt->first];
            nameStrLen = nameStr.length();
            file.write(reinterpret_cast<char*>(&nameStrLen), 4);
            file.write(nameStr.c_str(), nameStrLen);
            DimVec<int>* nDOF = feSpaceIt->first->getBasisFcts()->getNumberOfDofs();
            //
            for(int j = 1; j < nDOF->getSize(); j++)
            {
              posDOFs = (*nDOF)[j];
              file.write(reinterpret_cast<char*>(&posDOFs), 4);
            }
            for(int j = nDOF->getSize(); j < 4 ; j++)
            {
              posDOFs = 0;
              file.write(reinterpret_cast<char*>(&posDOFs), 4);
            }
            posDOFs = (*nDOF)[0];
            file.write(reinterpret_cast<char*>(&posDOFs), 4);
            //
            for(size_t j = 0; j < feSpaceIt->second.size(); j++)
              feSpaceNumOfVecs[feSpaceIt->second[j]] = i;
          }

          //vector table
          for(i = 0; i < vecs.size(); i++)
          {
            nameStr = vecs[i]->getName();
            nameStrLen = nameStr.length();
            file.write(reinterpret_cast<char*>(&nameStrLen), 4);
            file.write(nameStr.c_str(), nameStrLen);
            file.write(reinterpret_cast<char*>(&feSpaceNumOfVecs[i]), 4);
            file.write(dataformat.c_str(), 4);
          }
          return headerLen;
        }

        void writeAux(string filename, Mesh* mesh,
                      vector<DOFVector<double>*> vecs,
                      bool /*writeParallel*/, // NOTE: param writeParallel only in parallel mode
                      Cpsformat cps,
                      string dataformat)
        {
          FUNCNAME("Arh3Writer::detail::writeAux()");

          //initialization
          ofstream file;
          file.open(filename.c_str(), ios::out | ios::binary | ios::trunc);

          map<const FiniteElemSpace*, vector<int>> sortedFeSpaces;
          map<const FiniteElemSpace*, vector<int>>::iterator feSpaceIt;
          vector<pair<int, int>> macroSize; // (uncompressed size, compressed size)

          DegreeOfFreedom globalDof;
          size_t i = 0, j = 0;

          for(i = 0; i < vecs.size(); i++)
          {
            sortedFeSpaces[vecs[i]->getFeSpace()].push_back(i);
          }
          vector<std::set<DegreeOfFreedom>> visited(sortedFeSpaces.size());
          pair<std::set<DegreeOfFreedom>::iterator,bool> ret;
          //file header information
          int headerLen = writeHeader(file, mesh, vecs, sortedFeSpaces, cps, dataformat);

          //macro elements information
          MeshStructure elementStructure;
          vector<vector<double>> values(vecs.size());
          int32_t macroElIndex = -1;

          TraverseStack stack;
          ElInfo* elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
          while (elInfo)
          {
            if (elInfo->getLevel() == 0)
            {
              if (macroElIndex != -1)
              {
                elementStructure.commit();
                macroSize.push_back(writeMacroElement(file, elementStructure, values, sortedFeSpaces, cps, dataformat));
              }
              elementStructure.clear();

              macroElIndex = elInfo->getElement()->getIndex();

              for (j = 0; j < vecs.size(); j++)
                values[j].clear();

              for (j = 0; j < sortedFeSpaces.size(); j++)
                visited[j].clear();
            }
            elementStructure.insertElement(elInfo->getElement()->isLeaf());

            if (elInfo->getElement()->isLeaf())
            {
              int valuePos = 0;
              for(feSpaceIt = sortedFeSpaces.begin(), i = 0; feSpaceIt != sortedFeSpaces.end(); feSpaceIt++, i++)
              {
                DOFAdmin* admin = feSpaceIt->first->getAdmin();
                TEST_EXIT(admin)("the DOFAdmin of DOFVector is null, this should not happen.\n");

                int n0, nd, nd0;

                if ((nd = admin->getNumberOfDofs(VERTEX)))
                {
                  int vertices = mesh->getGeo(VERTEX);
                  nd0 = admin->getNumberOfPreDofs(VERTEX);
                  n0 = mesh->getNode(VERTEX); //
                  for (int n = 0; n < vertices; n++)
                    for(int d = 0; d < nd; d++)
                    {
                      globalDof = elInfo->getElement()->getDof(n0 + n, nd0 + d);
                      ret = visited[i].insert(globalDof);
                      if(ret.second)
                      {
                        for(j = 0 ; j < feSpaceIt->second.size(); j++)
                        {
                          values[valuePos + j].push_back((*vecs[feSpaceIt->second[j]])[globalDof]);
                        }
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
                          for(j = 0 ; j < feSpaceIt->second.size(); j++)
                          {
                            values[valuePos + j].push_back((*vecs[feSpaceIt->second[j]])[globalDof]);
                          }
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
                          for(j = 0 ; j < feSpaceIt->second.size(); j++)
                          {
                            values[valuePos + j].push_back((*vecs[feSpaceIt->second[j]])[globalDof]);
                          }
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
                      for(j = 0 ; j < feSpaceIt->second.size(); j++)
                      {
                        values[valuePos + j].push_back((*vecs[feSpaceIt->second[j]])[globalDof]);
                      }
                    }
                  }
                }
                valuePos += feSpaceIt->second.size();
              }//feSpace loop
            }//isLeaf

            elInfo = stack.traverseNext(elInfo);
          }
          // And write the last macro element to file.
          TEST_EXIT_DBG(macroElIndex != -1)("Should not happen!\n");
          elementStructure.commit();
          macroSize.push_back(writeMacroElement(file, elementStructure, values, sortedFeSpaces, cps, dataformat));
          TEST_EXIT(macroSize.size() == (unsigned)mesh->getNumberOfMacros())("Should not happen.\n");
          //reset the macro positions in file
          setMacrosPos(file, headerLen, macroSize);
          file.close();
          MSG("ARH file written to: %s\n", filename.c_str());
        }

        void setMacrosPos(ofstream& file, int headerLen,
                          vector<pair<int, int>>& macroSize)
        {
          file.seekp(34);
          long pos = 0;
          int startPos = headerLen;
          for(size_t i = 0; i < macroSize.size(); i++)
          {
            pos = file.tellp();
            file.seekp(pos + 4);
            file.write(reinterpret_cast<char*>(&startPos), 4);
            file.write(reinterpret_cast<char*>(&macroSize[i].first), 4);
            startPos += macroSize[i].second;
          }
          pos = file.tellp();
          file.seekp(pos + 4);
          file.write(reinterpret_cast<char*>(&startPos), 4);
        }

        pair<int, int> writeMacroElement(ofstream& file,
                                         MeshStructure& code,
                                         vector<vector<double>>& values,
                                         map<const FiniteElemSpace*, vector<int>>& feSpaces,
                                         Cpsformat /*cps*/, // NOTE: param cps only with HAVE_COMPRESSION
                                         string dataformat)
        {
          stringstream dataStream(ios::out | ios::in | ios::binary);

          uint32_t nStructureCodes = code.getCode().size();
          dataStream.write(reinterpret_cast<char*>(&nStructureCodes), 4);

          uint32_t codeSize = code.getNumElements();
          dataStream.write(reinterpret_cast<char*>(&codeSize), 4);

          dataStream.write(reinterpret_cast<char*>(&(const_cast<vector<uint64_t>&>(code.getCode())[0])),
                           8 * nStructureCodes);

          if (values.size() > 0)
          {
            int moreSize = 0, valuePos = 0;
            map<const FiniteElemSpace*, vector<int>>::iterator it;
            for(it = feSpaces.begin(); it != feSpaces.end(); it++)
            {
              uint32_t nValuesPerVector = values[valuePos].size();
              dataStream.write(reinterpret_cast<char*>(&nValuesPerVector), 4);
              moreSize += 4;

              for (size_t i = 0; i < it->second.size(); i++)
                moreSize += writeValues(dataStream, dataformat, values[valuePos + i]);

              valuePos += it->second.size();
            }
          }
          stringstream tmp(ios::out | ios::in | ios::binary);
          boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
#ifdef HAVE_COMPRESSION
          switch(cps)
          {
          case ZLIB:
            in.push(boost::iostreams::zlib_compressor());
            break;
          case BZIP2:
            in.push(boost::iostreams::bzip2_compressor());
            break;
          case NONE:
            break;
          default:
            MSG("NOT correct compression flag.\n");
          }
#endif
          in.push(dataStream);
          boost::iostreams::copy(in, tmp);
          file << tmp.rdbuf();
          return make_pair(dataStream.str().length(), tmp.str().length());
        }

        template<typename T>
        int writeValues(stringstream& file, vector<double>& values)
        {
          size_t size = values.size();
          T* data = new T[size];
          for (size_t i = 0; i < size; i++)
            data[i] = static_cast<T>(values[i]);

          file.write(reinterpret_cast<char*>(&data[0]), sizeof(T) * size);
          delete [] data;
          return sizeof(T) * size;
        }

        int writeValues(stringstream& file,
                        string dataformat,
                        vector<double>& values)
        {
          FUNCNAME("Arh3Writer::detail::writeValues()");
          int newsize = 0;
          std::map<string,Valformat>::const_iterator it = dataformatMap.find(dataformat);

          TEST_EXIT(it != dataformatMap.end())("Wrong data format.\n");

          switch(it->second)
          {
          case SI08:
            newsize = writeValues<int8_t>(file, values);
            break;
          case SI16:
            newsize = writeValues<int16_t>(file, values);
            break;
          case SI32:
            newsize = writeValues<int32_t>(file, values);
            break;
          case SI64:
            newsize = writeValues<int64_t>(file, values);
            break;
          case UI08:
            newsize = writeValues<uint8_t>(file, values);
            break;
          case UI16:
            newsize = writeValues<uint16_t>(file, values);
            break;
          case UI32:
            newsize = writeValues<uint32_t>(file, values);
            break;
          case UI64:
            newsize = writeValues<uint64_t>(file, values);
            break;
          case SF32:
            newsize = writeValues<float>(file, values);
            break;
          case SF64:
            newsize = writeValues<double>(file, values);
            break;
          default:
            ERROR_EXIT("Wrong data format.\n");
          }
          return newsize;
        }
      }//end namespace detail
    } // end namespace Arh3Writer
  }
} // end namespace io, AMDiS
