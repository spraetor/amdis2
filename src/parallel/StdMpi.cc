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


#include "StdMpi.h"

using namespace std;

namespace AMDiS { namespace Parallel {

  MPI_Datatype StdMpiHelper<int>::mpiDataType = MPI_INT;
  MPI_Datatype StdMpiHelper<vector<int> >::mpiDataType = MPI_INT;
  MPI_Datatype StdMpiHelper<std::set<int> >::mpiDataType = MPI_INT;
  MPI_Datatype StdMpiHelper<vector<std::set<int> > >::mpiDataType = MPI_INT;
  MPI_Datatype StdMpiHelper<std::set<std::pair<int, int> > >::mpiDataType = MPI_INT;
  MPI_Datatype StdMpiHelper<vector<double> >::mpiDataType = MPI_DOUBLE;
  MPI_Datatype StdMpiHelper<vector<vector<double> > >::mpiDataType = MPI_DOUBLE;
  MPI_Datatype StdMpiHelper<vector<MeshStructure> >::mpiDataType = MPI_UNSIGNED_LONG_LONG;
  MPI_Datatype StdMpiHelper<vector<AtomicBoundary> >::mpiDataType = MPI_INT;
  MPI_Datatype StdMpiHelper<map<BoundaryType, map<DegreeOfFreedom, DegreeOfFreedom> > >::mpiDataType = MPI_INT;
  MPI_Datatype StdMpiHelper<vector<std::pair<int, int> > >::mpiDataType = MPI_INT;
  MPI_Datatype StdMpiHelper<vector<WorldVector<double> > >::mpiDataType = MPI_DOUBLE;
  MPI_Datatype StdMpiHelper<vector<WorldVector<int> > >::mpiDataType = MPI_INT;
  MPI_Datatype StdMpiHelper<vector<Vector<int> > >::mpiDataType = MPI_INT;
  MPI_Datatype StdMpiHelper<vector<WorldMatrix<double> > >::mpiDataType = MPI_DOUBLE;
  MPI_Datatype StdMpiHelper<vector<WorldVector<WorldVector<double> > > >::mpiDataType = MPI_DOUBLE;
  MPI_Datatype StdMpiHelper<vector<WorldMatrix<int> > >::mpiDataType = MPI_INT;
  MPI_Datatype StdMpiHelper<vector<Matrix<int> > >::mpiDataType = MPI_INT;
  MPI_Datatype StdMpiHelper<map<WorldVector<double>, int> >::mpiDataType = MPI_DOUBLE;
  MPI_Datatype StdMpiHelper<map<int, WorldVector<double> > >::mpiDataType = MPI_DOUBLE;
  MPI_Datatype StdMpiHelper<vector<vector<WorldVector<double> > > >::mpiDataType = MPI_DOUBLE;

  // T = int

  int StdMpiHelper<int>::getBufferSize(int &data)
  {
    return 1;
  }

  void StdMpiHelper<int>::createBuffer(int &data, int *buf)
  {
    buf[0] = data;
  }

  void StdMpiHelper<int>::makeFromBuffer(int &data, int *buf, int bufSize)
  {
    data = buf[0];
  }



  // T = vector<int>

  int StdMpiHelper<vector<int> >::getBufferSize(vector<int> &data)
  {
    return data.size();
  }

  int StdMpiHelper<vector<int> >::getBufferSize(vector<const int*> &data)
  {
    return data.size();
  }

  void StdMpiHelper<vector<int> >::createBuffer(vector<int> &data, int *buf)
  {
    for (unsigned int i = 0; i < data.size(); i++)
      buf[i] = data[i];
  }

  void StdMpiHelper<vector<int> >::makeFromBuffer(vector<int> &data, int *buf, int bufSize)
  {
    data.resize(bufSize);

    for (int i = 0; i < bufSize; i++)
      data[i] = buf[i];
  }



  // T = std::set<int>

  int StdMpiHelper<std::set<int> >::getBufferSize(std::set<int> &data)
  {
    return data.size();
  }

  void StdMpiHelper<std::set<int> >::createBuffer(std::set<int> &data, int *buf)
  {
    int i = 0;
    for (std::set<int>::iterator it = data.begin(); it != data.end(); ++it)
      buf[i++] = *it;
  }

  void StdMpiHelper<std::set<int> >::makeFromBuffer(std::set<int> &data, int *buf, int bufSize)
  {
    data.clear();

    for (int i = 0; i < bufSize; i++)
      data.insert(buf[i]);
  }


  // T = vector<std::set<int> >

  int StdMpiHelper<vector<std::set<int> > >::getBufferSize(vector<std::set<int> > &data)
  {
    int dataSize = 1;
    for (unsigned int i = 0; i < data.size(); i++)
      dataSize += data[i].size() + 1;

    return dataSize;
  }

  void StdMpiHelper<vector<std::set<int> > >::createBuffer(vector<std::set<int> > &data, 
							   int *buf)
  {
    int counter = 1;
    buf[0] = data.size();

    for (unsigned int i = 0; i < data.size(); i++) {
      buf[counter++] = data[i].size();
      for (std::set<int>::iterator it = data[i].begin(); it != data[i].end(); ++it)
	buf[counter++] = *it;
    }
  }

  void StdMpiHelper<vector<std::set<int> > >::makeFromBuffer(vector<std::set<int> > &data, 
							     int *buf, 
							     int bufSize)
  {
    FUNCNAME("StdMpiHelper<vector<set<int> > >::makeFromBuffer()");

    int counter = 1;
    data.resize(buf[0]);

    for (int i = 0; i < buf[0]; i++) {
      data[i].clear();
      int setSize = buf[counter++];
      TEST_EXIT_DBG(setSize <= bufSize - counter)
	("Should not happen: %d %d %d\n", setSize, bufSize, counter);

      for (int j = 0; j < setSize; j++)
	data[i].insert(buf[counter++]);
    }

    TEST_EXIT(counter == bufSize)("Counter is %d, but buffer size is %d!\n",
				  counter, bufSize);
  }

  // T = std::set<std::pair<int, int> > >
  
  int StdMpiHelper<std::set<std::pair<int, int> > >::getBufferSize(std::set<std::pair<int, int> > &data)
  {
    return data.size() * 2;
  }
  
  void StdMpiHelper<std::set<std::pair<int, int> > >::createBuffer(std::set<std::pair<int, int> > &data,
								   int *buf)
  {
    int i = 0;
      
    for(std::set<std::pair<int, int> >::iterator it = data.begin(); it != data.end(); ++it) {
      buf[i++] = it->first;
      buf[i++] = it->second;
    }
  }
  
  void StdMpiHelper<std::set<std::pair<int, int> > >::makeFromBuffer(std::set<std::pair<int, int> > &data,
								     int *buf,
								     int bufSize)
  {
    data.clear();
	
    for (int i = 0; i < bufSize; i += 2) 
      data.insert(std::make_pair(buf[i], buf[i + 1]));
  }
								     
  

  // T = vector<double>

  int StdMpiHelper<vector<double> >::getBufferSize(vector<double> &data)
  {
    return data.size();
  }

  void StdMpiHelper<vector<double> >::createBuffer(vector<double> &data, 
						   double *buf)
  {
    for (unsigned int i = 0; i < data.size(); i++)
      buf[i] = data[i];
  }

  void StdMpiHelper<vector<double> >::makeFromBuffer(vector<double> &data, double *buf, int bufSize)
  {
    data.resize(bufSize);

    for (int i = 0; i < bufSize; i++)
      data[i] = buf[i];
  }



  // T = vector<vector<double> >

  int StdMpiHelper<vector<vector<double> > >::getBufferSize(vector<vector<double> > &data)
  {
    int size = 1;

    for (unsigned int i = 0; i < data.size(); i++)
      size += data[i].size() + 1;
    
    return size;
  }

  void StdMpiHelper<vector<vector<double> > >::createBuffer(vector<vector<double> > &data, 
							    double *buf)
  {
    buf[0] = data.size();
    int counter = 1;

    for (unsigned int i = 0; i < data.size(); i++) {
      buf[counter++] = data[i].size();
      for (unsigned int j = 0; j < data[i].size(); j++)
	buf[counter++] = data[i][j];
    }
  }

  void StdMpiHelper<vector<vector<double> > >::makeFromBuffer(vector<vector<double> > &data, 
							      double *buf, 
							      int bufSize)
  {
    data.resize(static_cast<unsigned int>(buf[0]));
    int counter = 1;

    for (unsigned int i = 0; i < data.size(); i++) {
      data[i].resize(static_cast<unsigned int>(buf[counter++]));
      
      for (unsigned int j = 0; j < data[i].size(); j++)
	data[i][j] = buf[counter++];
    }
  }



  // T = vector<MeshStructure>

  int StdMpiHelper<vector<MeshStructure> >::getBufferSize(vector<MeshStructure> &data)
  {
    FUNCNAME("StdMpiHelper<vector<MeshStructure> >::getBufferSize()");

    int size = 0;
    for (unsigned int i = 0; i < data.size(); i++)
      size += data[i].getCode().size() + 2;

    return size;
  }

  void StdMpiHelper<vector<MeshStructure> >::createBuffer(vector<MeshStructure> &data, uint64_t *buf)
  {    
    FUNCNAME("StdMpiHelper<vector<MeshStructure> >::createBuffer()");

    int pos = 0;
    for (unsigned int i = 0; i < data.size(); i++) {
      buf[pos++] = data[i].getCode().size();
      buf[pos++] = data[i].getNumElements();

      for (unsigned int j = 0; j < data[i].getCode().size(); j++)
	buf[pos++] = data[i].getCode()[j];
    }
  }

  void StdMpiHelper<vector<MeshStructure> >::makeFromBuffer(vector<MeshStructure> &data, 
							    uint64_t *buf, int bufSize)
  {
    FUNCNAME("StdMpiHelper<vector<MeshStructure> >::makeFromBuffer()");

    int pos = 0;

    while (pos < bufSize) {
      int codeSize = buf[pos++];
      int nElements = buf[pos++];
      vector<uint64_t> code;
      code.resize(codeSize);
      for (int i = 0; i < codeSize; i++)
	code[i] = buf[pos++];

      MeshStructure meshCode;
      meshCode.init(code, nElements);

      data.push_back(meshCode);
    }	
  }



  // T = vector<AtomicBoundary>

  int StdMpiHelper<vector<AtomicBoundary> >::getBufferSize(vector<AtomicBoundary> &data)
  {
    return data.size() * 6;
  }

  void StdMpiHelper<vector<AtomicBoundary> >::createBuffer(vector<AtomicBoundary> &data, int *buf)
  {
    for (unsigned int i = 0; i < data.size(); i++) {
      buf[i * 6] = data[i].rankObj.elIndex;
      buf[i * 6 + 1] = data[i].rankObj.subObj;
      buf[i * 6 + 2] = data[i].rankObj.ithObj;
      buf[i * 6 + 3] = data[i].neighObj.elIndex;
      buf[i * 6 + 4] = data[i].neighObj.subObj;
      buf[i * 6 + 5] = data[i].neighObj.ithObj;
    }
  }

  void StdMpiHelper<vector<AtomicBoundary> >::makeFromBuffer(vector<AtomicBoundary> &data, int *buf, int bufSize)
  {
    if (bufSize == 0)
      return;

    TEST_EXIT(bufSize % 6 == 0)("This should not happen!\n");    

    data.resize(bufSize / 6);
    for (int i = 0; i < bufSize / 6; i++) {
      data[i].rankObj.elIndex = buf[i * 6];
      data[i].rankObj.subObj = static_cast<GeoIndex>(buf[i * 6 + 1]);
      data[i].rankObj.ithObj = buf[i * 6 + 2];
      data[i].neighObj.elIndex = buf[i * 6 + 3];
      data[i].neighObj.subObj = static_cast<GeoIndex>(buf[i * 6 + 4]);
      data[i].neighObj.ithObj = buf[i * 6 + 5];
    }
  }



  // T = map<BoundaryType, map<DegreeOfFreedom, DegreeOfFreedom> >

  int StdMpiHelper<map<BoundaryType, map<DegreeOfFreedom, DegreeOfFreedom> > >::getBufferSize(map<BoundaryType, map<DegreeOfFreedom, DegreeOfFreedom> > &data)
  {
    int size = 1;

    for (map<BoundaryType, map<DegreeOfFreedom, DegreeOfFreedom> >::iterator it = data.begin();
	 it != data.end(); ++it) {
      size += 2 + it->second.size() * 2;
    }

    return size;
  }

  void StdMpiHelper<map<BoundaryType, map<DegreeOfFreedom, DegreeOfFreedom> > >::createBuffer(map<BoundaryType, map<DegreeOfFreedom, DegreeOfFreedom> > &data, int *buf)
  {
    buf[0] = data.size();
    int counter = 1;

    for (map<BoundaryType, map<DegreeOfFreedom, DegreeOfFreedom> >::iterator it = data.begin();
	 it != data.end(); ++it) {
      buf[counter++] = it->first;
      buf[counter++] = it->second.size();

      for (map<DegreeOfFreedom, DegreeOfFreedom>::iterator it2 = it->second.begin();
	   it2 != it->second.end(); ++it2) {
	buf[counter++] = it2->first;
	buf[counter++] = it2->second;
      }
    }
  }

  void StdMpiHelper<map<BoundaryType, map<DegreeOfFreedom, DegreeOfFreedom> > >::makeFromBuffer(map<BoundaryType, map<DegreeOfFreedom, DegreeOfFreedom> > &data, int *buf, int bufSize)
  {
    data.clear();

    if (bufSize == 0)
      return;

    int counter = 1;
    
    for (int i = 0; i < buf[0]; i++) {
      BoundaryType bound = buf[counter++];
      map<DegreeOfFreedom, DegreeOfFreedom> dofs;

      int nDofs = buf[counter++];
      for (int j = 0; j < nDofs; j++) {
	DegreeOfFreedom dof0, dof1;
	dof0 = buf[counter++];
	dof1 = buf[counter++];
	dofs[dof0] = dof1;
      }

      data[bound] = dofs;
    }

    TEST_EXIT(bufSize == counter)("Should not happen!\n");
  }



  // T = vector<std::pair<int, int> >

  int StdMpiHelper<vector<std::pair<int, int> > >::getBufferSize(vector<std::pair<int, int> > &data)
  {
    return data.size() * 2;
  }

  void StdMpiHelper<vector<std::pair<int, int> > >::createBuffer(vector<std::pair<int, int> > &data, int *buf)
  {
    for (unsigned int i = 0; i < data.size(); i++) {
      buf[i * 2] = data[i].first;
      buf[i * 2 + 1] = data[i].second;
    }
  }

  void StdMpiHelper<vector<std::pair<int, int> > >::makeFromBuffer(vector<std::pair<int, int> > &data, 
									int *buf, int bufSize)
  {    
    if (bufSize == 0)
      return;

    TEST_EXIT(bufSize % 2 == 0)("This should not happen!\n");

    data.clear();
    data.reserve(bufSize / 2);

    for (int i = 0; i < (bufSize / 2); i++)
      data.push_back(std::make_pair(buf[i * 2], buf[i * 2 + 1]));
  }



  // T = vector<WorldVector<double> >

  int StdMpiHelper<vector<WorldVector<double> > >::getBufferSize(vector<WorldVector<double> > &data)
  {
    return data.size() * Global::getGeo(WORLD);
  }

  void StdMpiHelper<vector<WorldVector<double> > >::createBuffer(vector<WorldVector<double> > &data, double *buf)
  {
    int dimOfWorld = Global::getGeo(WORLD);
    int pos = 0;
    for (unsigned int i = 0; i < data.size(); i++)
      for (int j = 0; j < dimOfWorld; j++)
	buf[pos++] = data[i][j];
  }

  void StdMpiHelper<vector<WorldVector<double> > >::makeFromBuffer(vector<WorldVector<double> > &data, double *buf, int bufSize)
  {
    int dimOfWorld = Global::getGeo(WORLD);
    TEST_EXIT(bufSize % Global::getGeo(WORLD) == 0)("This should not happen!\n");

    int pos = 0;
    data.resize(bufSize / Global::getGeo(WORLD));
    for (unsigned int i = 0; i < data.size(); i++)
      for (int j = 0; j < dimOfWorld; j++)
	data[i][j] = buf[pos++];
  }



  // T = vector<WorldVector<int> >

  int StdMpiHelper<vector<WorldVector<int> > >::getBufferSize(vector<WorldVector<int> > &data)
  {
    return data.size() * Global::getGeo(WORLD);
  }

  void StdMpiHelper<vector<WorldVector<int> > >::createBuffer(vector<WorldVector<int> > &data, int *buf)
  {
    int dimOfWorld = Global::getGeo(WORLD);
    int pos = 0;
    for (unsigned int i = 0; i < data.size(); i++)
      for (int j = 0; j < dimOfWorld; j++)
	buf[pos++] = data[i][j];
  }

  void StdMpiHelper<vector<WorldVector<int> > >::makeFromBuffer(vector<WorldVector<int> > &data, int *buf, int bufSize)
  {
    int dimOfWorld = Global::getGeo(WORLD);
    TEST_EXIT(bufSize % Global::getGeo(WORLD) == 0)("This should not happen!\n");

    int pos = 0;
    data.resize(bufSize / Global::getGeo(WORLD));
    for (unsigned int i = 0; i < data.size(); i++)
      for (int j = 0; j < dimOfWorld; j++)
	data[i][j] = buf[pos++];
  }



  // T = vector<Vector<int> >

  int StdMpiHelper<vector<Vector<int> > >::getBufferSize(vector<Vector<int> > &data)
  {
    return data.size() > 0 ? data.size() * data[0].getSize() : 0;
  }

  void StdMpiHelper<vector<Vector<int> > >::createBuffer(vector<Vector<int> > &data, int *buf)
  {
    int pos = 0;
    buf[pos++] = data.size() > 0 ? data[0].getSize() : 0;
    for (unsigned int i = 0; i < data.size(); i++)
      for (int j = 0; j < data[i].getSize(); j++)
	buf[pos++] = data[i][j];
  }

  void StdMpiHelper<vector<Vector<int> > >::makeFromBuffer(vector<Vector<int> > &data, int *buf, int bufSize)
  {
    int pos = 0;
    int size = buf[pos++];
    data.resize(bufSize / size);
    for (unsigned int i = 0; i < data.size(); i++)
      for (int j = 0; j < size; j++)
	data[i][j] = buf[pos++];
  }


  // T = vector<WorldMatrix<double> >

  int StdMpiHelper<vector<WorldMatrix<double> > >::getBufferSize(vector<WorldMatrix<double> > &data)
  {
    return data.size() * sqr(Global::getGeo(WORLD));
  }

  void StdMpiHelper<vector<WorldMatrix<double> > >::createBuffer(vector<WorldMatrix<double> > &data, double *buf)
  {
    int dimOfWorld = Global::getGeo(WORLD);
    int pos = 0;
    for (unsigned int i = 0; i < data.size(); i++)
      for (int j = 0; j < dimOfWorld; j++)
	for (int k = 0; k < dimOfWorld; k++)
	  buf[pos++] = data[i][j][k];
  }

  void StdMpiHelper<vector<WorldMatrix<double> > >::makeFromBuffer(vector<WorldMatrix<double> > &data, double *buf, int bufSize)
  {
    int dimOfWorld = Global::getGeo(WORLD);
    TEST_EXIT(bufSize % sqr(Global::getGeo(WORLD)) == 0)("This should not happen!\n");

    int pos = 0;
    data.resize(bufSize / sqr(dimOfWorld));
    for (unsigned int i = 0; i < data.size(); i++)
      for (int j = 0; j < dimOfWorld; j++)
	for (int k = 0; k < dimOfWorld; k++)
	  data[i][j][k] = buf[pos++];
  }


  // T = vector<WorldMatrix<double> >

  int StdMpiHelper<vector<WorldVector<WorldVector<double> > > >::getBufferSize(vector<WorldVector<WorldVector<double> > > &data)
  {
    return data.size() * sqr(Global::getGeo(WORLD));
  }

  void StdMpiHelper<vector<WorldVector<WorldVector<double> > > >::createBuffer(vector<WorldVector<WorldVector<double> > > &data, double *buf)
  {
    int dimOfWorld = Global::getGeo(WORLD);
    int pos = 0;
    for (unsigned int i = 0; i < data.size(); i++)
      for (int j = 0; j < dimOfWorld; j++)
	for (int k = 0; k < dimOfWorld; k++)
	  buf[pos++] = data[i][j][k];
  }

  void StdMpiHelper<vector<WorldVector<WorldVector<double> > > >::makeFromBuffer(vector<WorldVector<WorldVector<double> > > &data, double *buf, int bufSize)
  {
    int dimOfWorld = Global::getGeo(WORLD);
    TEST_EXIT(bufSize % sqr(Global::getGeo(WORLD)) == 0)("This should not happen!\n");

    int pos = 0;
    data.resize(bufSize / sqr(dimOfWorld));
    for (unsigned int i = 0; i < data.size(); i++)
      for (int j = 0; j < dimOfWorld; j++)
	for (int k = 0; k < dimOfWorld; k++)
	  data[i][j][k] = buf[pos++];
  }



  // T = vector<WorldMatrix<int> >

  int StdMpiHelper<vector<WorldMatrix<int> > >::getBufferSize(vector<WorldMatrix<int> > &data)
  {
    return data.size() * sqr(Global::getGeo(WORLD));
  }

  void StdMpiHelper<vector<WorldMatrix<int> > >::createBuffer(vector<WorldMatrix<int> > &data, int *buf)
  {
    int dimOfWorld = Global::getGeo(WORLD);
    int pos = 0;
    for (unsigned int i = 0; i < data.size(); i++)
      for (int j = 0; j < dimOfWorld; j++)
	for (int k = 0; k < dimOfWorld; k++)
	  buf[pos++] = data[i][j][k];
  }

  void StdMpiHelper<vector<WorldMatrix<int> > >::makeFromBuffer(vector<WorldMatrix<int> > &data, int *buf, int bufSize)
  {
    int dimOfWorld = Global::getGeo(WORLD);
    TEST_EXIT(bufSize % sqr(Global::getGeo(WORLD)) == 0)("This should not happen!\n");

    int pos = 0;
    data.resize(bufSize / sqr(dimOfWorld));
    for (unsigned int i = 0; i < data.size(); i++)
      for (int j = 0; j < dimOfWorld; j++)
	for (int k = 0; k < dimOfWorld; k++)
	  data[i][j][k] = buf[pos++];
  }


  // T = vector<Matrix<int> >

  int StdMpiHelper<vector<Matrix<int> > >::getBufferSize(vector<Matrix<int> > &data)
  {
    return data.size() > 0 ? data.size() * data[0].getNumRows() * data[0].getNumCols() + 2 : 0;
  }

  void StdMpiHelper<vector<Matrix<int> > >::createBuffer(vector<Matrix<int> > &data, int *buf)
  {
//     int dimOfWorld = Global::getGeo(WORLD);
    int pos = 0;
    buf[pos++] = data.size() > 0 ? data[0].getNumRows() : 0;
    buf[pos++] = data.size() > 0 ? data[0].getNumCols() : 0;
    for (unsigned int i = 0; i < data.size(); i++)
      for (int j = 0; j < data[i].getNumRows(); j++)
	for (int k = 0; k < data[i].getNumCols(); k++)
	  buf[pos++] = data[i][j][k];
  }

  void StdMpiHelper<vector<Matrix<int> > >::makeFromBuffer(vector<Matrix<int> > &data, int *buf, int bufSize)
  {
    int pos = 0;
    int n_rows = buf[pos++];
    int n_cols = buf[pos++];
    data.resize(bufSize / (n_rows * n_cols));
    for (unsigned int i = 0; i < data.size(); i++) {
      data[i].resize(n_rows, n_cols);
      for (int j = 0; j < n_rows; j++)
	for (int k = 0; k < n_cols; k++)
	  data[i][j][k] = buf[pos++];
    }
  }



  // T = map<WorldVector<double>, int>

  int StdMpiHelper<map<WorldVector<double>, int> >::getBufferSize(map<WorldVector<double>, int> &data)
  {
    return data.size() * (Global::getGeo(WORLD) + 1);
  }

  void StdMpiHelper<map<WorldVector<double>, int> >::createBuffer(map<WorldVector<double>, int> &data, double* buf)
  {
    int i = 0;
    for (map<WorldVector<double>, int>::iterator it = data.begin();
	 it != data.end(); ++it) {
      for (int j = 0; j < Global::getGeo(WORLD); j++)
	buf[i++] = it->first[j];
      buf[i++] = static_cast<double>(it->second);      
    }
  }

  void StdMpiHelper<map<WorldVector<double>, int> >::makeFromBuffer(map<WorldVector<double>, int> &data, double* buf, int bufSize)
  {
    if (bufSize == 0)
      return;

    int oneEntrySize = Global::getGeo(WORLD) + 1;
    int nEntry = bufSize / oneEntrySize;

    TEST_EXIT(bufSize % oneEntrySize == 0)("This should not happen!\n");

    data.clear();
    int i = 0;
    WorldVector<double> coords;

    for (int j = 0; j < nEntry; j++) {
      for (int k = 0; k < Global::getGeo(WORLD); k++)
	coords[k] = buf[i++];
      data[coords] = static_cast<int>(buf[i++]);
    }
  }

  // T = map<WorldVector<double>, int>

  int StdMpiHelper<map<int, WorldVector<double> > >::getBufferSize(map<int, WorldVector<double> > &data)
  {
    return data.size() * (Global::getGeo(WORLD) + 1);
  }

  void StdMpiHelper<map<int, WorldVector<double> > >::createBuffer(map<int, WorldVector<double> > &data, double* buf)
  {
    int i = 0;
    for (map<int, WorldVector<double> >::iterator it = data.begin();
	 it != data.end(); ++it) {
      for (int j = 0; j < Global::getGeo(WORLD); j++)
	buf[i++] = it->second[j];
      buf[i++] = static_cast<double>(it->first);      
    }
  }

  void StdMpiHelper<map<int, WorldVector<double> > >::makeFromBuffer(map<int, WorldVector<double> > &data, double* buf, int bufSize)
  {
    if (bufSize == 0)
      return;

    int oneEntrySize = Global::getGeo(WORLD) + 1;
    int nEntry = bufSize / oneEntrySize;

    TEST_EXIT(bufSize % oneEntrySize == 0)("This should not happen!\n");

    data.clear();
    int i = 0;
    WorldVector<double> coords;

    for (int j = 0; j < nEntry; j++) {
      for (int k = 0; k < Global::getGeo(WORLD); k++)
	coords[k] = buf[i++];
      data[static_cast<int>(buf[i++])] = coords;
    }
  }

  // T = vector<vector<WorldVector<double> > >

  int StdMpiHelper<vector<vector<WorldVector<double> > > >::getBufferSize(vector<vector<WorldVector<double> > > &data)
  {
    int result = 1;

    for (unsigned int i = 0; i < data.size(); i++)
      result += 1 + (data[i].size() * Global::getGeo(WORLD));

    return result;
  }

  void StdMpiHelper<vector<vector<WorldVector<double> > > >::createBuffer(vector<vector<WorldVector<double> > > &data, 
									  double* buf)
  {
    int counter = 0;

    buf[counter++] = static_cast<double>(data.size());

    for (unsigned int i = 0; i < data.size(); i++) {
      buf[counter++] = data[i].size();
      for (unsigned int j = 0; j < data[i].size(); j++)
	for (int k = 0; k < Global::getGeo(WORLD); k++) 
	  buf[counter++] = data[i][j][k];
    }
  }

  void StdMpiHelper<vector<vector<WorldVector<double> > > >::makeFromBuffer(vector<vector<WorldVector<double> > > &data, 
									    double* buf, int bufSize)
  {
    int counter = 0;
    data.resize(static_cast<int>(buf[counter++]));

    for (unsigned int i = 0; i < data.size(); i++) {
      data[i].resize(static_cast<int>(buf[counter++]));
      for (unsigned int j = 0; j < data[i].size(); j++) 
	for (int k = 0; k < Global::getGeo(WORLD); k++) 
	  data[i][j][k] = buf[counter++];
    }

    TEST_EXIT_DBG(counter == bufSize)("There is something very wrong!\n");
  }

} }
