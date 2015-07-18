#include <fstream>
#include <stdint.h>
#include <iostream>

#include "Arh2Writer.h"
#include "Mesh.h"
#include "MeshStructure.h"
#include "Traverse.h"
#include "DOFVector.h"
#include "../Arh2Writer.h"

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#ifdef HAVE_COMPRESSION
#include <boost/iostreams/filter/zlib.hpp>
#endif

namespace AMDiS { namespace io {
  
  using namespace std;

  namespace Arh2Writer
  {
    namespace detail
    {
      void write(string filename, 
  			DOFVector<double>* vec0, 
  			DOFVector<double>* vec1,
  			DOFVector<double>* vec2,
			bool writeParallel)
      {
        vector<DOFVector<double>*> vecs(0);
        if (vec0 != NULL)
          vecs.push_back(vec0);
        if (vec1 != NULL)
          vecs.push_back(vec1);
        if (vec2 != NULL)
          vecs.push_back(vec2);

        write(filename, NULL, vecs, writeParallel);
      }
      
      void write(std::string filename, 
		 Mesh* mesh, 
		 std::vector<DOFVector<double>*> vecs,
		 bool writeParallel)
      {
	FUNCNAME("Arh2Writer::detail::write()");
	
	std::set<string> nameSet;
	pair<std::set<string>::iterator,bool> ret;
	
	for(size_t i = 0; i < vecs.size(); i++)
	{
	  TEST_EXIT(vecs[i] != NULL)("Vecs[%i] is NULL. Please check.\n", i);
	  ret = nameSet.insert(vecs[i]->getName());
	  TEST_EXIT(ret.second)("DOFVectors in vecs cannot have idential name. Please check.\n");
	} 
	//if mesh exists, the meshes in vecs should be the same.
	if(mesh) 
	{
	  for(size_t i = 0; i < vecs.size(); i++)
	  {
	    TEST_EXIT(mesh == vecs[i]->getFeSpace()->getMesh())
	    ("The mesh of DOFVector %i in vecs is not equal to the second parameter.\n", i);
	  }
	  writeAux(filename, mesh, vecs, writeParallel); 
	}
	//multiple meshes are allowed here.
	else
	{
	  if(vecs.empty())
	  {
	    WARNING("There is nothing to be writen.\n");
	    return;
	  }
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
	      writeAux(newfilename, splitedVecs[0]->getFeSpace()->getMesh(), splitedVecs, writeParallel); 
	    }
	  }
        }
      }
  
      int writeHeader(ofstream& file,
                               Mesh *mesh,
                               vector<DOFVector<double>*> vecs,
                               map<const FiniteElemSpace*, vector<int> >& feSpaces)
      {
        FUNCNAME("Arh2Writer::detail::writeHeader()");
        TEST_EXIT(file.is_open())("the file is not open. should not happen.\n");
   
        uint32_t namesLen = 0;
        for(size_t i = 0; i < vecs.size(); i++)
          namesLen += vecs[i]->getName().length();
    
        uint32_t nValueVectors = vecs.size();
        uint32_t nFeSpaces = feSpaces.size();
	uint32_t nMacroElements = mesh->getNumberOfMacros();
	
	uint32_t nMacro = 0;
	TraverseStack st;
	ElInfo *el = st.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL);
	while (el) {
	  nMacro++;
	  el = st.traverseNext(el);
	}
	
        uint32_t dow = mesh->getGeo(WORLD);
        uint32_t dim = mesh->getDim();
        uint32_t headerLen = 34 +                    //fixed part of header 
			     nMacroElements * 8 +    //macroElemnts table
			     nFeSpaces * 16 + 	     //feSpaces table
			     namesLen + 	     //value vector table
			     nValueVectors * 8;      //also value vector table
        string typeId = "arh2", cps("null");
#ifdef HAVE_COMPRESSION
	cps = "zlib";
#endif
        uint8_t *major = const_cast<uint8_t*>(&(AMDiS::io::Arh2Writer::MAJOR)); 
        uint8_t *minor = const_cast<uint8_t*>(&(AMDiS::io::Arh2Writer::MINOR));
        
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
	file.write(cps.c_str(), 4);
	//macro table
	deque<MacroElement*>::const_iterator macroIter = mesh->firstMacroElement();
	while(macroIter != mesh->endOfMacroElements())
	{
	  uint32_t macroIndex = (*macroIter)->getIndex(), macroPos = 0;
	  file.write(reinterpret_cast<char*>(&macroIndex), 4);
	  file.write(reinterpret_cast<char*>(&macroPos), 4);
	  macroIter++;
	}
	
	map<const FiniteElemSpace*, vector<int> >::iterator feSpaceIt;
	vector<int> feSpaceNumOfVecs(vecs.size());
	uint32_t posDOFs = 0, vecNameLen = 0;
	string vecName("");
	size_t i = 0;
	
	//feSpace table
	for(feSpaceIt = feSpaces.begin(); feSpaceIt != feSpaces.end(); feSpaceIt++, i++)
	{
	  DimVec<int>* nDOF = feSpaceIt->first->getBasisFcts()->getNumberOfDofs();
	  for(int j = 0; j < nDOF->getSize(); j++)
          {
            posDOFs = (*nDOF)[j];
            file.write(reinterpret_cast<char*>(&posDOFs), 4);
          }
          for(size_t j = nDOF->getSize(); j < 4 ; j++)
          {
            posDOFs = 0;
            file.write(reinterpret_cast<char*>(&posDOFs), 4);
          }
          for(size_t j = 0; j < feSpaceIt->second.size(); j++)
	  {
	    feSpaceNumOfVecs[feSpaceIt->second[j]] = i;
	  }
	}
	
	//vector table
        for(i = 0; i < vecs.size(); i++)
        {
          vecName = vecs[i]->getName();
          vecNameLen = vecs[i]->getName().length();
          file.write(reinterpret_cast<char*>(&vecNameLen), 4);
          file.write(vecName.c_str(), vecNameLen);
	  file.write(reinterpret_cast<char*>(&feSpaceNumOfVecs[i]), 4);
        }
        return headerLen;
      }
 
      void writeAux(string filename, Mesh *mesh,
                        vector<DOFVector<double>*> vecs,
                        bool writeParallel)
      {
        FUNCNAME("Arh2Writer::detail::writeAux()");
	
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	if (writeParallel) {	  
	  int sPos = filename.find(".arh");
	  TEST_EXIT(sPos >= 0)("Failed to find file postfix!\n");
	  string name = filename.substr(0, sPos);      
	  filename = name + "-p" + std::to_string(MPI::COMM_WORLD.Get_rank()) + "-.arh";
	}
#endif
        //initialization
        ofstream file;
        file.open(filename.c_str(), ios::out | ios::binary | ios::trunc);

        map<const FiniteElemSpace*, vector<int> > sortedFeSpaces;
	map<const FiniteElemSpace*, vector<int> >::iterator feSpaceIt;
        vector<int> macroBlockSize;
        
        DegreeOfFreedom globalDof;  
        size_t i = 0, j = 0;
    
        for(i = 0; i < vecs.size(); i++)
	{
          sortedFeSpaces[vecs[i]->getFeSpace()].push_back(i);
	}
	vector<std::set<DegreeOfFreedom> > visited(sortedFeSpaces.size());
	pair<std::set<DegreeOfFreedom>::iterator,bool> ret;
        //file header information 
        writeHeader(file, mesh, vecs, sortedFeSpaces);
    
        //macro elements information
        MeshStructure elementStructure;
        vector<vector<double> > values(vecs.size());
        int32_t macroElIndex = -1;
    
        TraverseStack stack;
        ElInfo *elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
        while (elInfo) {
          if (elInfo->getLevel() == 0) {
            if (macroElIndex != -1) {
              elementStructure.commit();
	      macroBlockSize.push_back(writeMacroElement(file, elementStructure, values, sortedFeSpaces));
            }
	    elementStructure.clear();
 	
	    macroElIndex = elInfo->getElement()->getIndex();
	    for (j = 0; j < vecs.size(); j++) {
	      values[j].clear();
	    }
            for (j = 0; j < sortedFeSpaces.size(); j++)
            {
              visited[j].clear();       
            }
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

              if ((nd = admin->getNumberOfDofs(VERTEX)))  {                                                      
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
              if (mesh->getDim() > 1) {
                if ((nd = admin->getNumberOfDofs(EDGE)))  {                                                      
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
              if (mesh->getDim() == 3) {
                if ((nd = admin->getNumberOfDofs(FACE)))  {                                                      
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
              if ((nd = admin->getNumberOfDofs(CENTER)))  {                                                      
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
	macroBlockSize.push_back(writeMacroElement(file, elementStructure, values, sortedFeSpaces));
	TEST_EXIT(macroBlockSize.size() == (unsigned)mesh->getNumberOfMacros())("Should not happen.\n");
	//reset the macro positions in file
	setMacrosPos(file, macroBlockSize);
        file.close();
        MSG("ARH file written to: %s\n", filename.c_str());
      }
      
      void setMacrosPos(ofstream& file,
	                    vector<int>& macroBlockSize)
      {
	file.seekp(34);
	for(size_t i = 0; i < macroBlockSize.size(); i++)
	{
	  long pos = file.tellp();
	  file.seekp(pos + 4);
	  file.write(reinterpret_cast<char*>(&macroBlockSize[i]), 4);
	}
      }
  
      int writeMacroElement(ofstream &file, 
				    MeshStructure &code,
				    vector<vector<double> >& values,
                                    map<const FiniteElemSpace*, vector<int> >& feSpaces)
      {
	stringstream dataStream(ios::out | ios::in | ios::binary);
    
        uint32_t nStructureCodes = code.getCode().size();
        dataStream.write(reinterpret_cast<char*>(&nStructureCodes), 4);
    
        uint32_t codeSize = code.getNumElements();
        dataStream.write(reinterpret_cast<char*>(&codeSize), 4);

        dataStream.write(reinterpret_cast<char*>(&(const_cast<vector<uint64_t>&>(code.getCode())[0])), 
	       8 * nStructureCodes);
	
        if (values.size() > 0) {
	  int moreSize = 0, valuePos = 0;
          map<const FiniteElemSpace*, vector<int> >::iterator it; 
          for(it = feSpaces.begin(); it != feSpaces.end(); it++)
          {
            uint32_t nValuesPerVector = values[valuePos].size();
            dataStream.write(reinterpret_cast<char*>(&nValuesPerVector), 4);
	    moreSize += 4;
	    
            for (size_t i = 0; i < it->second.size(); i++)
	    {
              dataStream.write(reinterpret_cast<char*>(&(values[valuePos + i][0])), 8 * nValuesPerVector);
	      moreSize += 8 * nValuesPerVector;
	    }
            valuePos += it->second.size();
          }
        }
	stringstream tmp(ios::out | ios::in | ios::binary);
        boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
#ifdef HAVE_COMPRESSION
	in.push(boost::iostreams::zlib_compressor());
#endif
	in.push(dataStream);
	boost::iostreams::copy(in, tmp);
	file << tmp.rdbuf();
        return tmp.str().length();
      }
      
    }//end namespace detail
  } // end namespace Arh2Writer
} } // end namespace io, AMDiS