#include "Arh3Reader.h"
#include "Mesh.h"
#include "MeshStructure.h"
#include "Traverse.h"
#include "DOFVector.h"
#include "SystemVector.h"
#include "Debug.h"
#include "../Arh3Reader.h"
#include "Arh3Writer.h"

#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#ifdef HAVE_COMPRESSION
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#endif

namespace AMDiS { namespace io {

  using namespace std;

  namespace Arh3Reader
  {
    namespace detail
    {
      
      void firstRead(ifstream& file, string id, uint8_t major, uint8_t minor)
      {
	FUNCNAME("Arh3Reader::detail::firstRead");
	
	uint8_t major_ = 0;
	uint8_t minor_ = 0;
	string typeId(4, ' ');
    
	file.read(const_cast<char*>(typeId.data()), 4);
	TEST_EXIT(typeId == id)
	  ("Cannot read Arh format: this file is not \"%s\" format.\n", id.c_str());

	file.read(reinterpret_cast<char*>(&major_), 1);
	TEST_EXIT(major == major_)
	  ("Cannot read Arh format: Arh2Reader major version is %i, the file major version is %i. \n",
	  major, major_);

	file.read(reinterpret_cast<char*>(&minor_), 1);
	TEST_EXIT(minor <= minor_)
	  ("Cannot read Arh format: ArhReader minor version is %i is smaller than the file minor version %i.\n", 
	  minor, minor_);
      }
      
      void setDofValues(int macroElIndex, Mesh *mesh,
			     vector<vector<double> >& values, 
                             vector<DOFVector<double>*>& vecs,
                             vector<vector<int> >& feSpaces)
      {
	if(vecs.size() == 0)
	  return;

	std::set<int> unique;
	std::pair<std::set<int>::iterator,bool> ret;
	std::vector<std::set<int> > visited;
	std::vector<int> dofIndex;
	int globalDof = 0; 

	for(size_t i = 0; i < feSpaces.size(); i++)
	{
	  visited.push_back(unique);
	  dofIndex.push_back(0);
	}

	TraverseStack stack;
	ElInfo *elInfo = stack.traverseFirstOneMacro(mesh, macroElIndex, -1, 
						    Mesh::CALL_EVERY_EL_PREORDER);
	while (elInfo) {
	  Element *el = elInfo->getElement();
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

	      if ((nd = admin->getNumberOfDofs(VERTEX)))  {
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
	      if ((nd = admin->getNumberOfDofs(CENTER)))  {
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
                         Mesh *mesh,
                         vector<DOFVector<double>*> vecs,
                         bool byName)
       {
	FUNCNAME("Arh3Reader::detail::read()");
	
	using namespace ::AMDiS::io::Arh3Writer;

	// Get set of all macro elements in mesh.
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
      
	Cpsformat cps = NONE;
	uint32_t headerLen = 0;
	uint32_t nMacroElements = 0;
	uint32_t nValueVectors = 0;
	uint32_t nFeSpaces = 0;
	uint32_t dim = 0, dow = 0;
	int cpsflag = 0;

	// Read fixed header
	firstRead(file, "sarh", MAJOR, MINOR);
      
	file.read(reinterpret_cast<char*>(&headerLen), 4);
	file.read(reinterpret_cast<char*>(&dow), 4);
	file.read(reinterpret_cast<char*>(&dim), 4);
	file.read(reinterpret_cast<char*>(&nFeSpaces), 4);
	file.read(reinterpret_cast<char*>(&nValueVectors), 4);
	file.read(reinterpret_cast<char*>(&nMacroElements), 4);
	file.read(reinterpret_cast<char*>(&cpsflag), 4);
	cps = static_cast<Cpsformat>(cpsflag);
	
#ifdef HAVE_COMPRESSION
	TEST_EXIT(cps == ZLIB  || 
		  cps == BZIP2 ||
		  cps == NONE)
	  ("Cannot read Arh2 file. Currently only support zlib and bzip2 compression.\n");
#else
	TEST_EXIT(cps == NONE)
	  ("HAVE_COMPRESSION OFF. Cannot read compressed Arh2 file.\n");	
#endif
	  
	TEST_EXIT(dow == (unsigned)mesh->getGeo(WORLD))
	  ("Dow is not match!\n");
	TEST_EXIT(dim == (unsigned)mesh->getDim())
	  ("Dim is not match!\n");
	TEST_EXIT(nValueVectors >= vecs.size())
	  ("File %s has %d vector(s), which is less than the number of DOFVectors %i in vecs!\n",
	   filename.c_str(), nValueVectors, vecs.size());

	vector<int> vecsNameLen;
	vector<string> vecsName;
	
	vector<int> filesNameLen;
	vector<string> filesName;
	
	vector<int> vecsFeSpaceNum;
	vector<string> dataformat;
	vector<int> macroElIndex(nMacroElements);
	vector<pair<int,int> > macroElSize(nMacroElements); // pos, uncompressed size
	int fileSize = 0;
	
	vector<vector<int> > feSpaceDOFs;
	vector<string> AFEDfileName(nFeSpaces);
	vector<int> perDOFs(4, 0);
	vector<vector<int> > sortedFeSpaces(nFeSpaces);
	
	// Read macro table
	for(size_t i = 0; i < nMacroElements; i++)
	{
	  file.read(reinterpret_cast<char*>(&macroElIndex[i]), 4);
	  file.read(reinterpret_cast<char*>(&macroElSize[i].first), 4);
	  file.read(reinterpret_cast<char*>(&macroElSize[i].second), 4);
	}
	file.seekg(4, ios_base::cur);
	file.read(reinterpret_cast<char*>(&fileSize), 4);
	file.seekg(4, ios_base::cur);
	// Read feSpace table
	for(size_t i = 0; i < nFeSpaces; i++)
	{
	  uint32_t tmpInt = 0;
	  file.read(reinterpret_cast<char*>(&tmpInt), 4);
	  file.read(const_cast<char*>(AFEDfileName[i].data()), tmpInt);
	  for(int j = 0; j < 4; j++)
	  {
	    file.read(reinterpret_cast<char*>(&perDOFs[j]), 4);
	  }
	  feSpaceDOFs.push_back(perDOFs);
	}
	// Read value table
	for(size_t i = 0; i < nValueVectors; i++)
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
	  tmpString.resize(4, ' ');
	  file.read(const_cast<char*>(tmpString.data()), 4);
	  dataformat.push_back(tmpString);
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
	      
	      size_t i;
	      TEST_EXIT(tmpVecs[k]->getFeSpace()->getAdmin())
		("Vecs number %i has no DOFAdmin. Should not happen.\n", k);
	      DimVec<int>* nDOF = tmpVecs[k]->getFeSpace()->getBasisFcts()->getNumberOfDofs();
	      
	      for(i = 0; i < nValueVectors; i++) {
		
		if(tmpVecs[k]->getName() != vecsName[i])  
		  continue;
		
		bool matchdof = true;
		if ((*nDOF)[0] != feSpaceDOFs[vecsFeSpaceNum[i]][3])
		  matchdof = false;
		else {
		  for(int j = 1; j < nDOF->getSize(); j++) {
		    if((*nDOF)[j] != feSpaceDOFs[vecsFeSpaceNum[i]][j-1]){
		      matchdof = false;
		      break;
		    }
		  }
		}
		if(matchdof) {
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
	  for(size_t i = 0; i < vecs.size(); i++) {
	    if(vecs[i]) {
	      
	      TEST_EXIT(vecs[i]->getFeSpace()->getAdmin())
	      ("Vecs number %i has no DOFAdmin. Should not happen.\n", i);
    
	      DimVec<int>* nDOF = vecs[i]->getFeSpace()->getBasisFcts()->getNumberOfDofs();
	      
	      TEST_EXIT((*nDOF)[0] == feSpaceDOFs[vecsFeSpaceNum[i]][3])
		  ("The fespace of vec number %i is not equal to the correspond fespace.\n", i+1);
	      for(int j = 1; j < nDOF->getSize(); j++)
		TEST_EXIT((*nDOF)[j] == feSpaceDOFs[vecsFeSpaceNum[i]][j-1])
		  ("The fespace of vec number %i is not equal to the correspond fespace.\n", i+1);
	    }
	  }
	}
	// Read data: meshstructure and dof values
	for (size_t i = 0; i < nMacroElements; i++) {
	  stringstream dataStream(ios::out | ios::in | ios::binary);
	  int size = (i == nMacroElements - 1) ?  fileSize - macroElSize[i].first
	    : (macroElSize[i+1].first - macroElSize[i].first);
	    
	  char* buffer = new char[size];
	  file.read(buffer, size);
	  dataStream.write(buffer, size);
	  delete[] buffer;
#ifdef HAVE_COMPRESSION
	  stringstream tmp(ios::out | ios::in);
	  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
	  switch(cps)
	  {
	    case ZLIB:
	      in.push(boost::iostreams::zlib_decompressor());
	      break;
	    case BZIP2:
	      in.push(boost::iostreams::bzip2_decompressor());
	      break;
	    case NONE:
	      break;
	    default:
	      MSG("NOT correct compression flag.\n");
	  }
	  in.push(dataStream);
	  boost::iostreams::copy(in, tmp);
	  dataStream.str(tmp.str());
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

	  if (nValueVectors > 0) {
	    vector<vector<double> > values(nValueVectors);
	    int valuePos = 0;
	    for(unsigned int j = 0; j < nFeSpaces; j++)
	    {
	      uint32_t nValuesPerVector = 0;
	      dataStream.read(reinterpret_cast<char*>(&nValuesPerVector), 4);
	      
	      for(size_t k = 0; k < sortedFeSpaces[j].size(); k++)
	      {
		values[valuePos + k].resize(nValuesPerVector); 
		readValues(dataStream, dataformat[valuePos + k], values[valuePos + k]);
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
	delete refManager;
      }
      
      void readValues(stringstream& file, 
		      string dataformat, 
		      vector<double>& values)
      {
	using namespace ::AMDiS::io::Arh3Writer;
	
	std::map<string,Valformat>::const_iterator it = dataformatMap.find(dataformat);
	TEST_EXIT(it != dataformatMap.end())("Wrong data format.\n");
	
	switch(it->second) {
	  case SI08:
	    readValues<int8_t>(file, values);
	    break;
	  case SI16:
	    readValues<int16_t>(file, values);
	    break;
	  case SI32:
	    readValues<int32_t>(file, values);
	    break;	    
	  case SI64:
	    readValues<int64_t>(file, values);
	    break;	
	  case UI08:
	    readValues<uint8_t>(file, values);
	    break;	  
	  case UI16:
	    readValues<uint16_t>(file, values);
	    break;
	  case UI32:
	    readValues<uint32_t>(file, values);
	    break;	  
	  case UI64:
	    readValues<uint64_t>(file, values);
	    break;
	  case SF32:
	    readValues<float>(file, values);
	    break;	
	  case SF64:
	    readValues<double>(file, values);
	    break;
	  default:
	    ERROR_EXIT("Wrong data format.\n");
	}
      }
      
      template<typename T>
      void readValues(stringstream& file, vector<double>& values)
      {
	size_t size = values.size();
	T* data = new T[size];
	file.read(reinterpret_cast<char*>(&data[0]), sizeof(T) * size);
	
	for (size_t i = 0; i < size; i++)
	  values[i] = static_cast<double>(data[i]);
	delete data;
      }
    
      void readFile(string filename, Mesh *mesh,
			  vector<DOFVector<double>*> vecs,
			  bool writeParallel,
			  int nProcs,
			  bool byName)
      {
	FUNCNAME("Arh3Reader::detail::readFile()");
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
	if (writeParallel) {
	  int sPos = filename.find(".arh");
	  TEST_EXIT(sPos >= 0)("Failed to find file postfix!\n");
	  string name = filename.substr(0, sPos);

	  if (nProcs == -1) {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	    string procFilename = name + "-p" + std::to_string(MPI::COMM_WORLD.Get_rank()) + "-.arh";
	    read(procFilename, mesh, vecs, byName);
	    MSG("ARH file read from: %s\n", procFilename.c_str());
#else
	    ERROR_EXIT("Reading parallel ARH files in sequential computations requires to specify the number of nodes on which the ARH file was created!\n");
#endif
	  } else {
	    string parhfn = name + ".parh", filenameType = "";
	    int nProcs_ = 0, nMacros_ = 0, nMacros = 0;
	    vector<int> partition;
	    
	    bool parh = boost::filesystem::exists(parhfn);
	    if (parh)
	      readParallelFile(parhfn, filenameType, partition, nProcs_, nMacros_);
	    else {
	      for (; nProcs_ < nProcs + 1; nProcs_++) {
		string fn = name + "-p" + std::to_string(nProcs_) + "-.arh";
		if(!boost::filesystem::exists(fn)) break;
	      }
	    }
	    TEST_EXIT(nProcs_ == nProcs)
	      ("Number of arh files doesn't match number of processors: %d vs %d\n", nProcs_, nProcs);
	      
	    if (!parh) {
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS	      
	      if(MPI::COMM_WORLD.Get_rank() == 0)
#endif
	        for(int i = 0; i < nProcs; i++)
		  nMacros_ += readNumOfMacrosFromSgArh(filename, i);
	    }
	      
	    nMacros = mesh->getNumberOfMacros();
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	    Parallel::mpi::globalAdd(nMacros);
	    if(MPI::COMM_WORLD.Get_rank() == 0)
#endif	      
	    {
	      TEST_EXIT(nMacros == nMacros_)
		  ("Number of macro elements in parallel ARH files doesn't match to the current runtime. %d vs %d\n",
		  nMacros, nMacros_); 
	    }

	    if (!parh) {
	      for (int i = 0; i < nProcs; i++) {
		string procFilename = name + "-p" + std::to_string(i) + "-.arh";
		read(procFilename, mesh, vecs, byName);
	      }
	    } else {
	      std::set<int> needFiles;
	      deque<MacroElement*>::iterator it = mesh->firstMacroElement();
	      for (;it != mesh->endOfMacroElements(); it++)
		needFiles.insert(partition[(*it)->getIndex()]);
	      
	      std::set<int>::iterator it2 = needFiles.begin();
	      for (;it2 != needFiles.end(); it2++) {
		string procFilename = name + "-p" + std::to_string(*it2) + "-.arh";
		read(procFilename, mesh, vecs, byName);
	      }
	    }
	  }
	} else {
	  read(filename, mesh, vecs, byName);
	}
	MSG("ARH file read from: %s\n", filename.c_str());
      }
      
      void readMetaFromSgArh(std::string filename, int nProc, 
			   std::vector< std::set<std::pair<int, int> > >& data)
      {
	using namespace ::AMDiS::io::Arh3Writer;
	
	ifstream file;
	file.open(filename.c_str(), ios::in | ios::binary);
	TEST_EXIT(file.is_open())
	  ("Cannot open file %s\n", filename.c_str());
	  
	uint32_t dow = 0, dim = 0;
	uint32_t headerLen = 0;
	uint32_t nMacroElements = 0;
	uint32_t nValueVectors = 0;
	uint32_t nFeSpaces = 0;
	int cpsflag = 0, fileSize = 0;
	Cpsformat cps = NONE;    
	
	firstRead(file, "sarh", MAJOR, MINOR);

	file.read(reinterpret_cast<char*>(&headerLen), 4);
	file.read(reinterpret_cast<char*>(&dow), 4);
	file.read(reinterpret_cast<char*>(&dim), 4);
	file.read(reinterpret_cast<char*>(&nFeSpaces), 4);
	file.read(reinterpret_cast<char*>(&nValueVectors), 4);
	file.read(reinterpret_cast<char*>(&nMacroElements), 4);
	file.read(reinterpret_cast<char*>(&cpsflag), 4);
	cps = static_cast<Cpsformat>(cpsflag);
	
#ifdef HAVE_COMPRESSION
	TEST_EXIT(cps == ZLIB  || 
		  cps == BZIP2 ||
		  cps == NONE)
	  ("Cannot read Arh2 file. Currently only support zlib and bzip2 compression.\n");
#else
	TEST_EXIT(cps == NONE)
	  ("HAVE_COMPRESSION OFF. Cannot read compressed Arh2 file.\n");	
#endif
      
	vector<int> macroElIndex(nMacroElements);
	vector<pair<int, int> > macroElSize(nMacroElements);
	vector<vector<int> > sortedFeSpaces(nFeSpaces);
	vector<string> dataformat;
	  
	// Read macro table
	for(uint32_t i = 0; i < nMacroElements; i++) {
	  file.read(reinterpret_cast<char*>(&macroElIndex[i]), 4);
	  file.read(reinterpret_cast<char*>(&macroElSize[i].first), 4);
	  file.read(reinterpret_cast<char*>(&macroElSize[i].second), 4);
	}
	file.seekg(4, ios_base::cur);
	file.read(reinterpret_cast<char*>(&fileSize), 4);
	file.seekg(4, ios_base::cur);
	
	// Read feSpace table
	for(size_t i = 0; i < nFeSpaces; i++)
	{
	  uint32_t tmpInt = 0;
	  file.read(reinterpret_cast<char*>(&tmpInt), 4);
	  file.seekg(tmpInt, ios_base::cur);
	  file.seekg(16, ios_base::cur);
	}
	
	// Read value table
	for(uint32_t i = 0; i < nValueVectors; i++) {
	  string tmpString("");
	  uint32_t tmpInt = 0;
	  file.read(reinterpret_cast<char*>(&tmpInt), 4);
	  tmpString.resize(tmpInt, ' ');
	  file.read(const_cast<char*>(tmpString.data()), tmpInt); //
	  file.read(reinterpret_cast<char*>(&tmpInt), 4);
	  sortedFeSpaces[tmpInt].push_back(i);
	  tmpString.resize(4, ' ');
	  file.read(const_cast<char*>(tmpString.data()), 4);
	  dataformat.push_back(tmpString);
	} 
    
	for (uint32_t i = 0; i < nMacroElements; i++) {
	  stringstream dataStream(ios::out | ios::in | ios::binary);
	  int size = (i == nMacroElements - 1) ?  fileSize - macroElSize[i].first
	    : (macroElSize[i+1].first - macroElSize[i].first);
	    
	  char* buffer = new char[size];
	  file.read(buffer, size);
	  dataStream.write(buffer, size);
	  delete[] buffer;
	  
#ifdef HAVE_COMPRESSION      
	  stringstream tmp(ios::out | ios::in);
	  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
	  switch(cps)
	  {
	    case ZLIB:
	      in.push(boost::iostreams::zlib_decompressor());
	      break;
	    case BZIP2:
	      in.push(boost::iostreams::bzip2_decompressor());
	      break;
	    case NONE:
	      break;
	    default:
	      MSG("NOT correct compression flag.\n");
	  }
	  in.push(dataStream);
	  boost::iostreams::copy(in, tmp);
	  dataStream.str(tmp.str());
#endif	    
	  uint32_t nStructureCodes = 0;
	  uint32_t codeSize = 0;

	  dataStream.read(reinterpret_cast<char*>(&nStructureCodes), 4);
	  dataStream.read(reinterpret_cast<char*>(&codeSize), 4);

	  data[nProc].insert(make_pair(macroElIndex[i], codeSize));
	  
	  // We only need codeSize for each macro element data block, so skip all of the unnecessary.
          dataStream.seekg(8 * nStructureCodes, ios_base::cur); 

	  if (nValueVectors > 0) {
	    
	    int valuePos = 0;
	    for(uint32_t j = 0; j < nFeSpaces; j++) {
	      uint32_t nValuesPerVector = 0;
	      dataStream.read(reinterpret_cast<char*>(&nValuesPerVector), 4);
	      
	      for(size_t k = 0; k < sortedFeSpaces[j].size(); k++) {
		
		std::map<string,Valformat>::const_iterator it = dataformatMap.find(dataformat[valuePos + k]);
		TEST_EXIT(it != dataformatMap.end())("Wrong data format.\n");
		
		int unitsize = 0;
		switch(it->second) {
		  case SI08: 
		  case UI08: unitsize = 1;break;
		  case SI16:
		  case UI16: unitsize = 2;break;
		  case SI32:
		  case UI32:
		  case SF32: unitsize = 4;break;
		  case SI64:
		  case UI64:
		  case SF64: unitsize = 8;break;
		  default:ERROR_EXIT("Wrong data format.\n");
		}
		dataStream.seekg(unitsize * nValuesPerVector, ios_base::cur);
	      }
	      valuePos += sortedFeSpaces[j].size();
	    }
	  }
	}
      } // end readMetaFromSgArh
    
  
      int readNumOfMacrosFromSgArh(std::string filename, int nProc)
      {
	FUNCNAME("Arh3Reader::readHeaderSize");
	int sPos = filename.find(".arh");
	TEST_EXIT(sPos >= 0)("Failed to find file postfix!\n");
	
	if (nProc >= 0) {
	  string name = filename.substr(0, sPos);
	  filename = name + "-p" + std::to_string(nProc) + "-.arh";
	}
	
	ifstream file;
	  
	file.open(filename.c_str(), ios::in | ios::binary);
	TEST_EXIT(file.is_open())
	  ("Cannot open file %s\n", filename.c_str());
	  
	int nMacroElements = 0;
	firstRead(file, "sarh", MAJOR, MINOR);
	file.seekg(20, ios_base::cur); 
	file.read(reinterpret_cast<char*>(&nMacroElements), 4);
	file.close();
	
	return nMacroElements;
      }
        
      void readParallelFile(string filename, string& filenameType, vector<int>& partition, int& nFiles, int& nMacros)
      {
	ifstream file;
	file.open(filename.c_str(), ios::in | ios::binary);
	TEST_EXIT(file.is_open())
	  ("Cannot open file %s\n", filename.c_str());
	  
	firstRead(file, "parh", 1, 0);
	
	uint32_t macroFile_nl = 0;
	string macroFilename;
	
	filenameType.resize(4, ' ');
	
	file.read(reinterpret_cast<char*>(&nFiles), 4);
	file.read(const_cast<char*>(filenameType.data()), 4);
	file.read(reinterpret_cast<char*>(&nMacros), 4);
	file.read(reinterpret_cast<char*>(&macroFile_nl), 4);
	macroFilename.resize(macroFile_nl, ' ');
	file.read(const_cast<char*>(macroFilename.data()), macroFile_nl);
	
	uint32_t rank = 0;
	for (int i = 0; i < nMacros; i++) {
	  file.read(reinterpret_cast<char*>(&rank), 4);
	  partition.push_back(rank);
	}
      }

    } // end namespace detail
  } // end namespace Arh3Reader
} } // end namespace io, AMDiS


