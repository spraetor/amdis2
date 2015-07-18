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


#include <fstream>
#include <stdint.h>
#include <boost/filesystem.hpp>

#include "ArhReader.h"
#include "Mesh.h"
#include "MeshStructure.h"
#include "Traverse.h"
#include "DOFVector.h"
#include "Debug.h"
#include "detail/ArhReader.h"

namespace AMDiS { namespace io {
  
  namespace ArhReader
  {
    using namespace std;
  
    void readFile(string filename,
			DOFVector<double>* vec0,
			DOFVector<double>* vec1,
			DOFVector<double>* vec2,
			bool writeParallel,
			int nProcs)
    {
      vector<DOFVector<double>*> vecs(0);
      if (vec0)
	vecs.push_back(vec0);
      if (vec1)
	vecs.push_back(vec1);
      if (vec2)
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
      FUNCNAME("ArhReader::readMeta()");

      Mesh* mesh = NULL;
      for(size_t i = 0; i < vecs.size(); i++)
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
      ElInfo *elInfo = stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL);
      while (elInfo) {
	int macroElIndex = elInfo->getElement()->getIndex();
	TEST_EXIT(elInRank.count(macroElIndex))("Should not happen!\n");
	readArhFiles.insert(elInRank[macroElIndex]);
	elInfo = stack.traverseNext(elInfo);
      } 


      // === Read the individual arh files. ===

      boost::filesystem::path p(filename.c_str());
      boost::filesystem::path directory = p.parent_path();

      for (std::set<int>::iterator it = readArhFiles.begin();
	  it != readArhFiles.end(); ++it) {
	string arhFilename = directory.string() + "/" + arhPrefix;
	if (nProc == 1)
	  arhFilename += ".arh";
	else
	  arhFilename += "-p" + boost::lexical_cast<string>(*it) + "-.arh";

	MSG("ARH file read from: %s\n", arhFilename.c_str());
	detail::read(arhFilename, mesh, vecs);
      }
    }

    int readMetaData(string filename,
				map<int, int> &elInRank,
				map<int, int> &elCodeSize,
				string &arhPrefix)
    {
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

      // Maps to each macro element index the arh file index it is stored in.
      for (int i = 0; i < nProc; i++) {
	int tmp;
	file >> tmp;
	TEST_EXIT(tmp == i)("Should not happen!\n");
	int nMacroEl;
	file >> nMacroEl;
	for (int j = 0; j < nMacroEl; j++) {
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
      FUNCNAME("ArhReader::readMetaData()");

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
    
    int readNumOfValueVectors(string filename)
    {
      ifstream file;
      file.open(filename.c_str(), ios::in | ios::binary);
      
      string typeId = "";
      uint32_t nMacroElements = 0;
      uint32_t nValueVectors = 0;
      
      file.read(const_cast<char*>(typeId.data()), 4);
      file.read(reinterpret_cast<char*>(&nMacroElements), 4);
      file.read(reinterpret_cast<char*>(&nValueVectors), 4);
      
      file.close();
      
      return nValueVectors;
    }

    //the following three functions are identical to read/readfile except that they read from
    //a block of memory instead of from a file
    void readFromMemoryBlock(vector<char> &data, Mesh *mesh,
					DOFVector<double>* vec0,
					DOFVector<double>* vec1,
					DOFVector<double>* vec2,
					bool writeParallel,
					int nProcs)
    {
      uint32_t nValueVectors;
      memcpy(reinterpret_cast<char*>(&nValueVectors), &data[8], 4);
      vector<DOFVector<double>*> vecs(0);
      if (nValueVectors > 0)
	vecs.push_back(vec0);
      if (nValueVectors > 1)
	vecs.push_back(vec1);
      if (nValueVectors > 2)
	vecs.push_back(vec2);
      for (uint32_t i = 3; i < nValueVectors; i++)
	vecs.push_back(NULL);

      readFromMemoryBlock(data, mesh, vecs, writeParallel, nProcs);
    }
    
    void readFromMemoryBlock(vector<char> &data, Mesh *mesh,
					vector<DOFVector<double>*> vecs,
					bool writeParallel,
					int nProcs)
    {
      FUNCNAME("ArhReader::readFromMemoryBlock()");

      if (writeParallel) {
	ERROR_EXIT("Reading a hierarchy from several memory blocks is not implemented yet!\n");
      } else {
	detail::readBlock(data, mesh, vecs);
      }
    }
    
    int getNumValueVectors(string filename)
    {
      FUNCNAME("ArhReader::getNumValueVectors()")
      WARNING("this function is obsolete.\n");
      return readNumOfValueVectors(filename);
    }
    
    void readFile(string filename,
			    Mesh *mesh,
			    vector<DOFVector<double>*> vecs)
    {
      FUNCNAME("ArhReader::readFile()");
      WARNING("this function is obsolete.\n");
      detail::read(filename, mesh, vecs);
    }
    
    void read(string filename, Mesh *mesh,
			vector<DOFVector<double>*> vecs,
			bool writeParallel,
			int nProcs)
    {
      FUNCNAME("ArhReader::read()");
      WARNING("this function is obsolete.\n");
      readFile(filename, vecs, writeParallel, nProcs);
    }

    void read(string filename, Mesh *mesh,
			DOFVector<double>* vec0,
			DOFVector<double>* vec1,
			DOFVector<double>* vec2,
			bool writeParallel,
			int nProcs)
    {
      FUNCNAME("ArhReader::read()")
      WARNING("this function is obsolete.\n");
      readFile(filename, vec0, vec1, vec2, writeParallel, nProcs);
    }
    
    void setDofValues(int macroElIndex, Mesh *mesh,
				vector<double>& values, DOFVector<double>* vec)
    {
      FUNCNAME("ArhReader::seDofValues()");
      WARNING("this function is obsolete.\n");
      detail::setDofValues(macroElIndex, mesh, values, vec);
    }

    void readBlock(vector<char> &data,
			      Mesh *mesh,
			      vector<DOFVector<double>*> vecs)
    {
      FUNCNAME("ArhReader::readBlock()");
      WARNING("this function is obsolete.\n");
      detail::readBlock(data, mesh, vecs);
    }
    
    void readMeta(string filename,
			 Mesh *mesh,
			 vector<DOFVector<double>*> vecs)
    {
       FUNCNAME("ArhReader::readMeta()");
       WARNING("this function is obsolete.\n");
       readMeta(filename, vecs);
    }

    void readMeta(string filename,
		         Mesh* mesh,
			 DOFVector<double>* vec0,
			 DOFVector<double>* vec1,
			 DOFVector<double>* vec2)
    {
       FUNCNAME("ArhReader::readMeta()");
       WARNING("this function is obsolete.\n");
       readMeta(filename, vec0, vec1, vec2);
    }
  } // end namespace ArhReader
} } // end namespace io, AMDiS
