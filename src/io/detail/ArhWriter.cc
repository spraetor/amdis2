#include <fstream>
#include <stdint.h>

#include "ArhWriter.h"
#include "Mesh.h"
#include "MeshStructure.h"
#include "Traverse.h"
#include "DOFVector.h"

namespace AMDiS { namespace io {
  
  using namespace std;

  namespace ArhWriter
  {
    namespace detail
    {
      void write(string filename, Mesh *mesh,
			    DOFVector<double>* vec0, 
			    DOFVector<double>* vec1,
			    DOFVector<double>* vec2)
      {
	vector<DOFVector<double>*> vecs(0);
	if (vec0 != NULL)
	  vecs.push_back(vec0);
	if (vec1 != NULL)
	  vecs.push_back(vec1);
	if (vec2 != NULL)
	  vecs.push_back(vec2);

	write(filename, mesh, vecs);
      }
      
      void write(string filename, Mesh *mesh, 
			  vector<DOFVector<double>*> vecs,
			  bool writeParallel)
      {
	FUNCNAME("ArhWriter::detail::write()");

    #ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	if (writeParallel) {
	  using boost::lexical_cast;
	  
	  int sPos = filename.find(".arh");
	  TEST_EXIT(sPos >= 0)("Failed to find file postfix!\n");
	  string name = filename.substr(0, sPos);      
	  filename = name + "-p" + lexical_cast<string>(MPI::COMM_WORLD.Get_rank()) + "-.arh";
	}
    #endif

	ofstream file;
	file.open(filename.c_str(), ios::out | ios::binary | ios::trunc);

	string typeId = "arhP";
	file.write(typeId.c_str(), 4);

	uint32_t nMacroElements = 0;
	TraverseStack stack;
	ElInfo *elInfo = stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL);
	while (elInfo) {
	  nMacroElements++;
	  elInfo = stack.traverseNext(elInfo);
	}

	file.write(reinterpret_cast<char*>(&nMacroElements), 4);
	
	uint32_t nValueVectors = vecs.size();
	file.write(reinterpret_cast<char*>(&nValueVectors), 4);

	uint32_t nAllValues = 
	  (vecs.size() > 0 ? vecs[0]->getFeSpace()->getAdmin()->getUsedDofs() : 0);
	file.write(reinterpret_cast<char*>(&nAllValues), 4);

	MeshStructure elementStructure;
	vector<vector<double> > values(vecs.size());
	int32_t macroElIndex = -1;
	
	elInfo = stack.traverseFirst(mesh, -1, Mesh::CALL_EVERY_EL_PREORDER);
	while (elInfo) {
	  if (elInfo->getLevel() == 0) {
	    if (macroElIndex != -1) {
	      elementStructure.commit();
	      writeMacroElement(file, elementStructure, values, macroElIndex);
	    }

	    elementStructure.clear();
	    macroElIndex = elInfo->getElement()->getIndex();
	    for (size_t i = 0; i < vecs.size(); i++) {
	      values[i].clear();

	      for (int j = 0; j < mesh->getGeo(VERTEX); j++)
		values[i].push_back((*vecs[i])[elInfo->getElement()->getDof(
		  j, vecs[i]->getFeSpace()->getAdmin()->getNumberOfPreDofs(VERTEX)
		)]);
	    }
	  }

	  elementStructure.insertElement(elInfo->getElement()->isLeaf());
	  
	  if (!elInfo->getElement()->isLeaf()) {
	    for (size_t i = 0; i < vecs.size(); i++)
	      values[i].push_back((*vecs[i])[elInfo->getElement()->getChild(0)->getDof(
		mesh->getDim(), vecs[i]->getFeSpace()->getAdmin()->getNumberOfPreDofs(VERTEX)
	      )]);
	  }

	  elInfo = stack.traverseNext(elInfo);
	}

	// And write the last macro element to file.
	TEST_EXIT_DBG(macroElIndex != -1)("Should not happen!\n");
	elementStructure.commit();
	writeMacroElement(file, elementStructure, values, macroElIndex);
	
	file.close();

	MSG("ARH file written to: %s\n", filename.c_str());
      }
      
      void writeMacroElement(ofstream &file, 
				      MeshStructure &code,
				      vector<vector<double> >& values,
				      uint32_t elIndex)
      {
	file.write(reinterpret_cast<char*>(&elIndex), 4); 
	
	uint32_t nStructureCodes = code.getCode().size();
	file.write(reinterpret_cast<char*>(&nStructureCodes), 4);
	
	uint32_t codeSize = code.getNumElements();
	file.write(reinterpret_cast<char*>(&codeSize), 4);

	file.write(reinterpret_cast<char*>(&(const_cast<vector<uint64_t>&>(code.getCode())[0])), 
		  8 * nStructureCodes);

	if (values.size() > 0) {
	  uint32_t nValuesPerVector = values[0].size();
	  file.write(reinterpret_cast<char*>(&nValuesPerVector), 4);

	  for (size_t i = 0; i < values.size(); i++)
	    file.write(reinterpret_cast<char*>(&(values[i][0])), 8 * nValuesPerVector);
	}
      }
      
    }//end namespace detail
  } // end namespace ArhWriter
} } // end namespace io, AMDiS
