#include <fstream>
#include <stdint.h>
#include <boost/filesystem.hpp>

#include "ArhReader.h"
#include "Mesh.h"
#include "MeshStructure.h"
#include "Traverse.h"
#include "DOFVector.h"
#include "MacroElement.h"

namespace AMDiS { namespace io {

  using namespace std;
  
  namespace ArhReader
  {
    namespace detail
    {
      
      void setDofValues(int macroElIndex, Mesh *mesh,
				vector<double>& values, DOFVector<double>* vec)
      {
	bool macroElement = false;
	int valuePos = 0;
	TraverseStack stack;
	ElInfo *elInfo = stack.traverseFirstOneMacro(mesh, macroElIndex, -1,
						    Mesh::CALL_EVERY_EL_PREORDER);
	while (elInfo) {
	  if (!macroElement) {
	    Element *mEl = elInfo->getMacroElement()->getElement();
	    for (int i = 0; i < mesh->getGeo(VERTEX); i++)
	      (*vec)[mEl->getDof(i, vec->getFeSpace()->getAdmin()->getNumberOfPreDofs(VERTEX))] = values[valuePos++];
	    macroElement = true;
	  }

	  Element *el = elInfo->getElement();
	  if (!el->isLeaf())
	    (*vec)[el->getChild(0)->getDof(mesh->getDim(), 
					   vec->getFeSpace()->getAdmin()->getNumberOfPreDofs(VERTEX))] = values[valuePos++];

	  elInfo = stack.traverseNext(elInfo);
	}
      }
      
      void read(string filename,
			    Mesh *mesh,
			    vector<DOFVector<double>*> vecs)
      {
	FUNCNAME("ArhReader::detail::read()");

	// === Get set of all macro elements in mesh. ===
	std::set<int> macroInMesh;
	for (std::deque<MacroElement*>::iterator it = mesh->getMacroElements().begin();
	    it != mesh->getMacroElements().end(); ++it)
	  macroInMesh.insert((*it)->getIndex());


	RefinementManager *refManager = NULL;
	switch (mesh->getDim()) {
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

	string typeId(4, ' ');
	uint32_t nMacroElements = 0;
	uint32_t nValueVectors = 0;
	uint32_t nAllValues = 0;

	file.read(const_cast<char*>(typeId.data()), 4);
	TEST_EXIT(typeId == "arhP" || typeId == "arhS")
	  ("Cannot read Arh format: this file is not \"arhP\" or \"arhS\" format.\n");
	
	file.read(reinterpret_cast<char*>(&nMacroElements), 4);
	file.read(reinterpret_cast<char*>(&nValueVectors), 4);
	file.read(reinterpret_cast<char*>(&nAllValues), 4);

	TEST_EXIT(nValueVectors == vecs.size())
	  ("File %s has %d vectors, but %d DOFVectors are provided!\n",
	  filename.c_str(), nValueVectors, vecs.size());

	for (unsigned int i = 0; i < nMacroElements; i++) {
	  uint32_t elIndex = 0;
	  uint32_t nStructureCodes = 0;
	  uint32_t codeSize = 0;

	  file.read(reinterpret_cast<char*>(&elIndex), 4);
	  file.read(reinterpret_cast<char*>(&nStructureCodes), 4);
	  file.read(reinterpret_cast<char*>(&codeSize), 4);

	  vector<uint64_t> structureCode(nStructureCodes);
	  file.read(reinterpret_cast<char*>(&(structureCode[0])), 8 * nStructureCodes);

	  MeshStructure elementStructure;
	  elementStructure.init(structureCode, codeSize);
	  if (macroInMesh.count(elIndex))
	    elementStructure.fitMeshToStructure(mesh, refManager, false, elIndex);

	  if (nValueVectors > 0) {
	    uint32_t nValuesPerVector = 0;
	    file.read(reinterpret_cast<char*>(&nValuesPerVector), 4);

	    for (unsigned int j = 0; j < nValueVectors; j++) {
	      vector<double> values(nValuesPerVector);
	      file.read(reinterpret_cast<char*>(&(values[0])), 8 * nValuesPerVector);
	      if (vecs[j] != NULL) {
		if (macroInMesh.count(elIndex))
		  setDofValues(elIndex, mesh, values, vecs[j]);
	      }
	    }
	  }
	}

	file.close();
	delete refManager;
	
	MSG("ARH file read from: %s\n", filename.c_str());
      }
      
      void readBlock(vector<char> &data,
			      Mesh *mesh,
			      vector<DOFVector<double>*> vecs)
      {
	FUNCNAME("ArhReader::detail::readBlock()");

	// === Get set of all macro elements in mesh. ===
	std::set<int> macroInMesh;
	for (std::deque<MacroElement*>::iterator it = mesh->getMacroElements().begin();
					    it != mesh->getMacroElements().end(); ++it)
	  macroInMesh.insert((*it)->getIndex());


	RefinementManager *refManager = NULL;
	switch (mesh->getDim()) {
	case 2:
	  refManager = new RefinementManager2d();
	  break;
	case 3:
	  refManager = new RefinementManager3d();
	  break;
	default:
	  ERROR_EXIT("Should not happen!\n");
	}

	string typeId = "";
	uint32_t nMacroElements = 0;
	uint32_t nValueVectors = 0;
	uint32_t nAllValues = 0;
	uint32_t mem_index = 0;

	memcpy(const_cast<char*>(typeId.data()), &data[mem_index], 4);
	mem_index += 4;
	memcpy(reinterpret_cast<char*>(&nMacroElements), &data[mem_index], 4);
	mem_index += 4;
	memcpy(reinterpret_cast<char*>(&nValueVectors), &data[mem_index], 4);
	mem_index += 4;
	memcpy(reinterpret_cast<char*>(&nAllValues), &data[mem_index], 4);
	mem_index += 4;

	TEST_EXIT(nValueVectors == vecs.size())
	  ("data has %d vectors, but %d DOFVectors are provided!\n", nValueVectors, vecs.size());

	for (unsigned int i = 0; i < nMacroElements; i++) {
	  uint32_t elIndex = 0;
	  uint32_t nStructureCodes = 0;
	  uint32_t codeSize = 0;

	  memcpy(reinterpret_cast<char*>(&elIndex), &data[mem_index], 4);
	  mem_index += 4;
	  memcpy(reinterpret_cast<char*>(&nStructureCodes), &data[mem_index], 4);
	  mem_index += 4;
	  memcpy(reinterpret_cast<char*>(&codeSize), &data[mem_index], 4);
	  mem_index += 4;

	  vector<uint64_t> structureCode(nStructureCodes);
	  memcpy(reinterpret_cast<char*>(&structureCode[0]), &data[mem_index], 8 * nStructureCodes);
	  mem_index += 8 * nStructureCodes;

	  MeshStructure elementStructure;
	  elementStructure.init(structureCode, codeSize);
	  if (macroInMesh.count(elIndex) == 1)
	    elementStructure.fitMeshToStructure(mesh, refManager, false, elIndex);

	  if(nValueVectors > 0){
	    uint32_t nValuesPerVector = 0;
	    memcpy(reinterpret_cast<char*>(&nValuesPerVector), &data[mem_index], 4);
	    mem_index += 4;

	    for (unsigned int j = 0; j < nValueVectors; j++) {
	      vector<double> values(nValuesPerVector);
	      memcpy(reinterpret_cast<char*>(&values[0]), &data[mem_index], 8 * nValuesPerVector);
	      mem_index += 8 * nValuesPerVector;
	      if (vecs[j] != NULL) {
		if (macroInMesh.count(elIndex) == 1)
		  setDofValues(elIndex, mesh, values, vecs[j]);
	      }
	    }
	  }
	}

	delete refManager;
      }
      
      void readFile(string filename, Mesh *mesh,
			vector<DOFVector<double>*> vecs,
			bool writeParallel,
			int nProcs)
      {
	FUNCNAME("ArhReader::detail::readFile()");

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
	  WARNING("You haven't specified the target, no mesh or DOFVectors is given.\n");
	  return;
	}
	
	if (writeParallel) {
	  int sPos = filename.find(".arh");
	  TEST_EXIT(sPos >= 0)("Failed to find file postfix!\n");
	  string name = filename.substr(0, sPos);

	  if (nProcs == -1) {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	    string procFilename = name + "-p" + std::to_string(MPI::COMM_WORLD.Get_rank()) + "-.arh";
	    read(procFilename, mesh, vecs);
#else
	    ERROR_EXIT("Reading parallel ARH files in sequential computations requires to specify the number of nodes on which the ARH file was created!\n");
#endif
	  } else {
	    for (int i = 0; i < nProcs; i++) {
	      string procFilename = name + "-p" + std::to_string(i) + "-.arh";
	      read(procFilename, mesh, vecs);
	    }
	  }
	} else {
	  read(filename, mesh, vecs);
	}
      }
    }//end namespace detail
  } // end namespace ArhReader
} } // end namespace io, AMDiS
