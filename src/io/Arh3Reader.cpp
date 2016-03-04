#include "io/Arh3Reader.hpp"

#include <fstream>
#include <cstdint>

#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/zlib.hpp>

#include "io/detail/Arh3Reader.hpp"

namespace AMDiS
{
  namespace io
  {

    namespace Arh3Reader
    {
      using namespace std;

      void readByName(string filename,
                      DOFVector<double>& vec,
                      bool writeParallel,
                      int nProcs)
      {
        readByName(filename, &vec, writeParallel, nProcs);
      }

      void readByName(string filename,
                      DOFVector<double>* vec,
                      bool writeParallel,
                      int nProcs)
      {
        vector<DOFVector<double>*> vecs(0);
        vecs.push_back(vec);
        detail::readFile(filename, NULL, vecs, writeParallel, nProcs, true);
      }

      void readByName(string filename,
                      vector<DOFVector<double>*> vecs,
                      bool writeParallel,
                      int nProcs)
      {
        detail::readFile(filename, NULL, vecs, writeParallel, nProcs, true);
      }

      void readByName(string filename,
                      SystemVector* sysVec,
                      bool writeParallel,
                      int nProcs)
      {
        vector<DOFVector<double>*> vecs(0);
        for(int i = 0 ; i < sysVec->getSize(); i++)
        {
          vecs.push_back(sysVec->getDOFVector(i));
        }
        detail::readFile(filename, NULL, vecs, writeParallel, nProcs, true);
      }

      void readFile(string filename,
                    SystemVector* sysVec,
                    bool writeParallel,
                    int nProcs)
      {
        vector<DOFVector<double>*> vecs(0);
        for(int i = 0 ; i < sysVec->getSize(); i++)
        {
          vecs.push_back(sysVec->getDOFVector(i));
        }
        detail::readFile(filename, NULL, vecs, writeParallel, nProcs);
      }

      void readFile(string filename,
                    DOFVector<double>* vec0,
                    DOFVector<double>* vec1,
                    DOFVector<double>* vec2,
                    bool writeParallel,
                    int nProcs)
      {
        vector<DOFVector<double>*> vecs(0);

        if(vec0 || vec1 || vec2)
          vecs.push_back(vec0);
        if(vec1 || vec2)
          vecs.push_back(vec1);
        if(vec2)
          vecs.push_back(vec2);

        detail::readFile(filename, NULL, vecs, writeParallel, nProcs);
      }

      void readFile(string filename,
                    vector<DOFVector<double>*> vecs,
                    bool writeParallel,
                    int nProcs)
      {
        detail::readFile(filename, NULL, vecs, writeParallel, nProcs);
      }

      void readFile(string filename,
                    Mesh* mesh,
                    bool writeParallel,
                    int nProcs)
      {
        vector<DOFVector<double>*> vecs(0);
        detail::readFile(filename, mesh, vecs, writeParallel, nProcs);
      }

      int readNumOfValueVectors(string filename, bool writeParallel)
      {
        FUNCNAME("Arh3Reader::readNumOfValueVectors");

        ifstream file;

        if(writeParallel)
        {
          int sPos = filename.find(".arh");
          TEST_EXIT(sPos >= 0)("Failed to find file postfix!\n");
          string name = filename.substr(0, sPos);
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
          filename = name + "-p" + std::to_string(MPI::COMM_WORLD.Get_rank()) + "-.arh";
#endif
        }

        file.open(filename.c_str(), ios::in | ios::binary);
        TEST_EXIT(file.is_open())
        ("Cannot open file %s\n", filename.c_str());

        uint32_t nValueVectors = 0;

        detail::firstRead(file, "sarh", MAJOR, MINOR);
        long pos = file.tellg();
        file.seekg(pos + 16);
        file.read(reinterpret_cast<char*>(&nValueVectors), 4);

        file.close();
        return nValueVectors;
      }

      int readHeaderSize(string filename, bool writeParallel)
      {
        FUNCNAME("Arh3Reader::readHeaderSize");

        ifstream file;

        if (writeParallel)
        {
          int sPos = filename.find(".arh");
          TEST_EXIT(sPos >= 0)("Failed to find file postfix!\n");
          string name = filename.substr(0, sPos);
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
          filename = name + "-p" + std::to_string(MPI::COMM_WORLD.Get_rank()) + "-.arh";
#endif
        }

        file.open(filename.c_str(), ios::in | ios::binary);
        TEST_EXIT(file.is_open())
        ("Cannot open file %s\n", filename.c_str());

        uint32_t headerLen = 0;
        detail::firstRead(file, "sarh", MAJOR, MINOR);
        file.read(reinterpret_cast<char*>(&headerLen), 4);

        file.close();
        return headerLen;
      }

      bool isReadable(string filename, bool writeParallel)
      {
        FUNCNAME("Arh3Reader::isReadable");

        ifstream file;

        if(writeParallel)
        {
          int sPos = filename.find(".arh");
          TEST_EXIT(sPos >= 0)("Failed to find file postfix!\n");
          string name = filename.substr(0, sPos);
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
          filename = name + "-p" + std::to_string(MPI::COMM_WORLD.Get_rank()) + "-.arh";
#endif
        }

        file.open(filename.c_str(), ios::in | ios::binary);
        TEST_EXIT(file.is_open())
        ("Cannot open file %s\n", filename.c_str());

        string typeId(4, ' ');
        uint8_t major = 0, minor = 0;

        file.read(const_cast<char*>(typeId.data()), 4);
        file.read(reinterpret_cast<char*>(&major), 1);
        file.read(reinterpret_cast<char*>(&minor), 1);

        file.close();
        return (typeId == "sarh" && major == MAJOR && minor <= MINOR) ? true : false;
      }

      void readMeta(string filename,
                    DOFVector<double>* vec0,
                    DOFVector<double>* vec1,
                    DOFVector<double>* vec2)
      {
        vector<DOFVector<double>*> vecs;
        if (vec0)
          vecs.push_back(vec0);
        if (vec1)
          vecs.push_back(vec1);
        if (vec2)
          vecs.push_back(vec2);

        readMeta(filename, vecs);
      }

      void readMeta(string filename, vector<DOFVector<double>*> vecs)
      {
        FUNCNAME("Arh3Reader::readMeta()");

        Mesh* mesh = NULL;
        for (size_t i = 0; i < vecs.size(); i++)
        {
          if(vecs[i])
          {
            if(!mesh)
              mesh = vecs[i]->getFeSpace()->getMesh();
            else
              TEST_EXIT(mesh == vecs[i]->getFeSpace()->getMesh())
              ("The mesh of the DOFVectors should be the same for Reader because in one file there is only one mesh.\n");
          }
        }
        if(!mesh)
        {
          WARNING("You haven't specified the target.\n");
          return;
        }
        // === Read the meta arh file. ===

        string arhPrefix = "";
        map<int, int> elInRank;
        map<int, int> elCodeSize;
        int nProc = readMetaData(filename, elInRank, elCodeSize, arhPrefix);


        // === Check which arh files must be read by current rank. ===

        // Set of all file indices which should be read to restore all macro elements.
        std::set<int> readArhFiles;

        TraverseStack stack;
        ElInfo* elInfo = stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL);
        while (elInfo)
        {
          int macroElIndex = elInfo->getElement()->getIndex();
          TEST_EXIT(elInRank.count(macroElIndex))("Should not happen!\n");
          readArhFiles.insert(elInRank[macroElIndex]);
          elInfo = stack.traverseNext(elInfo);
        }


        // === Read the individual arh files. ===

        boost::filesystem::path p(filename.c_str());
        boost::filesystem::path directory = p.parent_path();

        for (std::set<int>::iterator it = readArhFiles.begin();
             it != readArhFiles.end(); ++it)
        {
          string arhFilename = directory.string() + "/" + arhPrefix;
          if (nProc == 1)
            arhFilename += ".arh";
          else
            arhFilename += "-p" + std::to_string(*it) + "-.arh";

          MSG("ARH2 file read from: %s\n", arhFilename.c_str());
          detail::read(arhFilename, mesh, vecs);
        }
      }

      int readMetaData(string filename,
                       map<int, int>& elInRank,
                       map<int, int>& elCodeSize,
                       string& arhPrefix)
      {
        FUNCNAME("Arh3Reader::readMetaData()");

        ifstream file;
        file.open(filename.c_str());
        TEST_EXIT(file.is_open())
        ("Cannot open arh meta file \"%s\"\n", filename.c_str());

        string readStr = "";
        file >> readStr;
        arhPrefix = "";
        file >> arhPrefix;
        int nProc;
        file >> nProc;

        MSG("found nProc %d in file %s \n",nProc,filename.c_str());

        // Maps to each macro element index the arh file index it is stored in.
        for (int i = 0; i < nProc; i++)
        {
          int tmp;
          file >> tmp;
          TEST_EXIT(tmp == i)("Should not happen!\n");
          int nMacroEl;
          file >> nMacroEl;
          for (int j = 0; j < nMacroEl; j++)
          {
            int elIndex, codeSize;
            file >> elIndex;
            file >> codeSize;
            elInRank[elIndex] = i;
            elCodeSize[elIndex] = codeSize;
          }
        }

        file.close();

        return nProc;
      }

      int readMetaData(string filename)
      {
        FUNCNAME("Arh3Reader::readMetaData()");

        ifstream file;
        file.open(filename.c_str());
        TEST_EXIT(file.is_open())
        ("Cannot open arh meta file \"%s\"\n", filename.c_str());

        string readStr = "";
        file >> readStr;
        file >> readStr;
        int nProc;
        file >> nProc;
        file.close();
        return nProc;
      }

      int readMetaFromArh(std::string filename,
                          std::map<int, int>& elInRank,
                          std::map<int, int>& elCodeSize)
      {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
        int nProcs = MPI::COMM_WORLD.Get_size();
#else
        int nProcs = -1;
#endif

        // check for consistency
        // That is, if there are no more or less files as nProcs
        int i = 0;
        for (; i < nProcs + 1; i++)
        {
          string fn = filename + "-p" + std::to_string(i) + "-.arh";
          if(!boost::filesystem::exists(fn)) break;
        }
        TEST_EXIT(i == nProcs)
        ("Number of arh files doesn't match number of processors \n");

        // data format: (rank ; (elIndex; elCodeSize) )
        vector<std::set<pair<int, int>>> data(nProcs);

        //collect data
        for (int i = 0; i < nProcs; i++)
        {
          string fn = filename + "-p" + std::to_string(i) + "-.arh";
          detail::readMetaFromSgArh(fn, i, data);
        }

        //make elInRank and elCodeSize-Map
        for (int i = 0; i < nProcs; i++)
        {
          for (std::set<pair<int, int>>::iterator it = data[i].begin();
               it != data[i].end(); ++it)
          {
            elInRank[it->first]=i;
            elCodeSize[it->first]=it->second;
          }
        }
        return nProcs;
      }

    } // end namespace Arh3Reader
  }
} // end namespace io, AMDiS
