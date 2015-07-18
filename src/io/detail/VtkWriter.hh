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



/** \file VtkVector.hh */

#ifndef AMDIS_VTKWRITER_DETAIL_HH
#define AMDIS_VTKWRITER_DETAIL_HH

#include "BasisFunction.h"
#include "io/DataCollector.h"

namespace AMDiS  { namespace io {
  
    namespace VtkWriter
    {
	
      template<typename T>
      void Aux::writeFileToStream(T &file)
      {
	int nVertices = (*dataCollector)[0]->getNumberVertices();
	int nElements = (*dataCollector)[0]->getNumberElements();
	int vertices = (*dataCollector)[0]->getMesh()->getGeo(VERTEX);

	if (dim == 2 && degree == 2) {
	  nVertices += (*dataCollector)[0]->getNumberInterpPoints();
	  nElements *= 4;
	} else if (dim == 2 && degree == 3) {
	  nVertices += (*dataCollector)[0]->getNumberInterpPoints();
	  nElements *= 9;
	} else if (dim == 2 && degree == 4) {
	  nVertices += (*dataCollector)[0]->getNumberInterpPoints();
	  nElements *= 16;
	}
	file << "<?xml version=\"1.0\"?>\n";
	file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
	file << "  <UnstructuredGrid>\n";
	file << "    <Piece NumberOfPoints=\"" << nVertices 
	    << "\" NumberOfCells=\"" << nElements << "\">\n";
	file << "      <Points>\n";
	file << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	
	writeVertexCoords(file);
	
	file << "        </DataArray>\n";
	file << "      </Points>\n";
	file << "      <Cells>\n";
	file << "        <DataArray type=\"Int32\" Name=\"offsets\">\n";
	
	for (int i = 0; i < nElements; i++)
	  file << " " << (i + 1) * vertices << "\n";

	file << "        </DataArray>\n";
	file << "        <DataArray type=\"UInt8\" Name=\"types\">\n";

	for (int i = 0; i < nElements; i++) {
	  switch (vertices) {
	  case 2:
	    file << " 3\n";
	    break;
	  case 3: 
	    file << " 5\n";
	    break;
	  case 4:
	    file << " 10\n";
	    break;		    
	  default:
	    break;
	  }	   
	}

	file << "        </DataArray>\n";
	file << "        <DataArray type=\"Int32\" Name=\"connectivity\">\n";

	writeConnectivity(file);
	
	file << "        </DataArray>\n";    
	file << "      </Cells>\n";
	file << "      <PointData>\n";
	
	int nValues = static_cast<int>(dataCollector->size());
	nValues = writeAsVector ? std::min(nValues,1) : nValues;	
	for (int i = 0; i < nValues; i++) {
	  file << "        <DataArray type=\"Float32\" Name=\"" 
	      << componentNames[i]
	      << "\" format=\"ascii\"";
	  if (writeAsVector && dataCollector->size() > 1)
	    file << " NumberOfComponents=\"" << std::max(3,static_cast<int>(dataCollector->size())) << "\"";
	  file << ">\n"; 

	  writeVertexValues(file, i);
	  
	  file << "        </DataArray>\n";
	}
	
	file << "      </PointData>\n";
	file << "    </Piece>\n";
	file << "  </UnstructuredGrid>\n";
	file << "</VTKFile>\n";
      }
      
      template<typename T>
      void Aux::writeFileToStreamAppended(T &file)
      {
	using namespace std;
	
	bstream.str("");
	
	int valueSize = writeAsVector ? 1 : static_cast<int>(dataCollector->size());
	valueSize = std::min(valueSize, static_cast<int>(dataCollector->size()));
	int nDataArray = 4 + valueSize;
	vector<string> dataBase64Vec(nDataArray, "");
	vector<size_t> dataOffsetVec(nDataArray, 0);
	string finalData = "";
	
	int nVertices = (*dataCollector)[0]->getNumberVertices();
	int nElements = (*dataCollector)[0]->getNumberElements();
	int vertices = (*dataCollector)[0]->getMesh()->getGeo(VERTEX);

	if (dim == 2 && degree == 2) {
	  nVertices += (*dataCollector)[0]->getNumberInterpPoints();
	  nElements *= 4;
	} else if (dim == 2 && degree == 3) {
	  nVertices += (*dataCollector)[0]->getNumberInterpPoints();
	  nElements *= 9;
	} else if (dim == 2 && degree == 4) {
	  nVertices += (*dataCollector)[0]->getNumberInterpPoints();
	  nElements *= 16;
	}
	//DataArray
	for (int i = 0; i < valueSize; i++) {
	  writeVertexValues(bstream, i);
	  dataBase64Vec[i] = getStreamData();
	}
	//DataArray
	writeVertexCoords(bstream);
 	dataBase64Vec[valueSize] = getStreamData();
	//DataArray
	writeConnectivity(bstream);
 	dataBase64Vec[valueSize + 1] = getStreamData();
	//DataArray
	uint32_t tmp32 = 0;
	for (int i = 1; i <= nElements; i++) {
	  tmp32 = i * vertices;
	  bstream << tmp32;
	}
 	dataBase64Vec[valueSize + 2] = getStreamData();
	//DataArray
	uint8_t tmp8 = 0;
	for (int i = 0; i < nElements; i++) {
	  switch (vertices) {
	  case 2:
	    tmp8 = 3;
	    break;
	  case 3: 
	    tmp8 = 5;
	    break;
	  case 4:
	    tmp8 = 10;
	    break;		    
	  default:
	    tmp8 = 0;
	    break;
	  }	   
	  if(tmp8)
	    bstream << tmp8;
	}
 	dataBase64Vec[valueSize + 3] = getStreamData();
	
	finalData = dataBase64Vec[0];
	for(int i = 1; i < nDataArray; i++)
	{
	  dataOffsetVec[i] = dataBase64Vec[i - 1].length() + dataOffsetVec[i - 1];
	  finalData += dataBase64Vec[i];
	}
	
	file << "<?xml version=\"1.0\"?>\n";
	file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"";
	if(format == APPENDED_COMPRESSED)
	  file << " compressor=\"vtkZLibDataCompressor\"";
	file << ">\n";
	file << "  <UnstructuredGrid>\n";
	file << "    <Piece NumberOfPoints=\"" << nVertices 
	    << "\" NumberOfCells=\"" << nElements << "\">\n";
	file << "      <PointData>\n";
	
	for (int i = 0; i < valueSize; i++) {
	  if(bstream.isHighPrecision())
	    file << "        <DataArray type=\"Float64\" Name=\"";
	  else
	    file << "        <DataArray type=\"Float32\" Name=\"";
	  
	  file << componentNames[i] 
	      << "\" format=\"appended\" offset=\"" << dataOffsetVec[i] << "\"";
	  if (writeAsVector && dataCollector->size() > 1)
	    file << " NumberOfComponents=\"" << std::max(3,static_cast<int>(dataCollector->size())) << "\"";
	  file << "/>\n"; 
	}
	
	file << "      </PointData>\n";
	file << "      <CellData>\n";
	file << "      </CellData>\n";
	file << "      <Points>\n";
	if(bstream.isHighPrecision())
	  file << "        <DataArray type=\"Float64\"";
	else
	  file << "        <DataArray type=\"Float32\"";
	file << " Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << dataOffsetVec[valueSize] << "\"/>\n";
	file << "      </Points>\n";
	file << "      <Cells>\n";
	file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"" << dataOffsetVec[valueSize + 1] << "\"/>\n";
	file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"" << dataOffsetVec[valueSize + 2] << "\"/>\n";
	file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << dataOffsetVec[valueSize + 3] << "\"/>\n";
	file << "      </Cells>\n";
	file << "    </Piece>\n";
	file << "  </UnstructuredGrid>\n";
	file << "  <AppendedData encoding=\"base64\">\n";
	file << "   _" << finalData << "\n";
	file << "  </AppendedData>\n";
	file << "</VTKFile>\n";
	
      }

      
      template<typename T>
      void Aux::writeVertexCoords(T &file)
      {
	DOFVector<std::list<VertexInfo> > *vertexInfos = 
	  (*dataCollector)[0]->getVertexInfos();
	DOFVector<std::list<VertexInfo> >::Iterator it(vertexInfos, USED_DOFS);
	int counter = 0;

	// For all DOFs of vertices, write the coordinates.
	for (it.reset(); !it.end(); ++it) {
	  // for all vertex infos of this DOF
	  for (std::list<VertexInfo>::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
	    it2->outputIndex = counter++;
	    detail::writeCoord(file, it2->coords);
	  }     
	}

	// For the second dim case, write also the interpolation points.
	if ((dim == 2) && (degree > 1)) {
	  DOFVector<std::list< WorldVector<double> > > *interpPointCoords = (*dataCollector)[0]->getInterpPointCoords();
	  DOFVector<std::list< WorldVector<double> > >::Iterator pointIt(interpPointCoords, USED_DOFS);
	  
	  counter = 0;
	  for (pointIt.reset(); !pointIt.end(); ++pointIt) {
	    for (std::list< WorldVector<double> >::iterator it2 = pointIt->begin(); 
		it2 != pointIt->end(); ++it2) {
	      counter++;
	      detail::writeCoord(file, *it2);
	    }
	  }
	}
      }


      template<typename T>
      void Aux::writeVertexValues(T &file, int componentNo)
      {
	if (writeAsVector) {
	  writeVertexValues(file);
	  return;
	}
	
	DOFVector<int> *interpPointInd = (*dataCollector)[componentNo]->getInterpPointInd();
	DOFVector<double> *values = (*dataCollector)[componentNo]->getValues();
	DOFVector<std::list<WorldVector<double> > > *dofCoords = (*dataCollector)[componentNo]->getDofCoords();

	DOFVector<int>::Iterator intPointIt(interpPointInd, USED_DOFS);
	DOFVector<double>::Iterator valueIt(values, USED_DOFS);
	DOFVector<std::list<WorldVector<double> > >::Iterator coordIt(dofCoords, USED_DOFS);

	file << std::scientific;
	file.precision(15);
	    
	// Write the values for all vertex DOFs.
	for (intPointIt.reset(), valueIt.reset(), coordIt.reset();
	    !intPointIt.end(); 
	    ++intPointIt, ++valueIt, ++coordIt) {

	  if (*intPointIt == -2) 
	    for (int i = 0; i < static_cast<int>(coordIt->size()); i++) 
	      file << " " << (fabs(*valueIt) < 1e-40 ? 0.0 : *valueIt) << "\n";      
	}        

	// For the second dim case, write also the values of the interpolation points.
	if ((dim == 2) && (degree > 1)) {
	  DOFVector<std::list<WorldVector<double> > >::Iterator 
	    interpCoordIt((*dataCollector)[componentNo]->getInterpPointCoords(), USED_DOFS);
	
	  for (intPointIt.reset(), valueIt.reset(), interpCoordIt.reset();
	      !intPointIt.end(); 
	      ++intPointIt, ++valueIt, ++interpCoordIt) {      
	    
	    if (*intPointIt >= 0) {
	      for (unsigned int i = 0; i < interpCoordIt->size(); i++)
		file << " " << (fabs(*valueIt) < 1e-40 ? 0.0 : *valueIt) << "\n";
	    }
	  }
	}    
      }
      
      


      template<typename T>
      void Aux::writeVertexValues(T &file)
      {	
	DOFVector<int> *interpPointInd = (*dataCollector)[0]->getInterpPointInd();
	std::vector<DOFVector<double>*> values(dataCollector->size());	
	DOFVector<std::list<WorldVector<double> > > *dofCoords = (*dataCollector)[0]->getDofCoords();
	
	for (size_t i = 0; i < dataCollector->size(); i++)
	  values[i] = (*dataCollector)[i]->getValues();
	
	DOFIterator<int> intPointIt(interpPointInd, USED_DOFS);
	DOFVectorIterator<double> valueIt(values, USED_DOFS);
	DOFIterator<std::list<WorldVector<double> > > coordIt(dofCoords, USED_DOFS);
	
	file << std::scientific;
	file.precision(15);
	    
	// Write the values for all vertex DOFs.
	for (intPointIt.reset(), valueIt.reset(), coordIt.reset(); 
	     !intPointIt.end(); 
	     ++intPointIt, ++valueIt, ++coordIt) {

	  if (*intPointIt == -2) 
	    for (int i = 0; i < static_cast<int>(coordIt->size()); i++) {
	      file << " " ;
	      print(*valueIt, file);
	      file << "\n";
	    }
	}        

	// For the second dim case, write also the values of the interpolation points.
	if ((dim == 2) && (degree > 1)) {
	  DOFVector< std::list<WorldVector<double> > >::Iterator 
	    interpCoordIt((*dataCollector)[0]->getInterpPointCoords(), USED_DOFS);
	
	  for (intPointIt.reset(), valueIt.reset(), interpCoordIt.reset();
	      !intPointIt.end(); 
	      ++intPointIt, ++valueIt, ++interpCoordIt) {      
	    
	    if (*intPointIt >= 0) {
	      for (unsigned int i = 0; i < interpCoordIt->size(); i++) {
		file << " ";
		print(*valueIt, file);
		file << "\n";
	      }
	    }
	  }
	}
      }

      template<typename T>
      void Aux::writeConnectivity(T &file)
      {	
	// For the second dim case, and if higher order Lagrange elements are used,
	// write the connectivity by extra functions.
	if ((dim == 2) && (degree == 2)) {
	  detail::writeConnectivity_dim2_degree2((*dataCollector)[0], file);
	} else if ((dim == 2) && (degree == 3)) {
	  detail::writeConnectivity_dim2_degree3((*dataCollector)[0], file);
	} else if ((dim == 2) && (degree == 4)) {
	  detail::writeConnectivity_dim2_degree4((*dataCollector)[0], file);   
	} else {
	  std::list<ElementInfo> *elements = (*dataCollector)[0]->getElementInfos();
	  std::list<ElementInfo>::iterator elementIt;
	  int vertices = (*dataCollector)[0]->getMesh()->getGeo(VERTEX);
	  
	  for (elementIt = elements->begin(); elementIt != elements->end(); ++elementIt) {
	    // for all vertices
	    for (int i = 0; i < vertices; i++) {
	      file << " " << elementIt->vertexInfo[i]->outputIndex;
	    }
	    file << "\n";
	  } 
	}
      }
      
// __________________________________________________________________________ //
      
      namespace detail
      {
	template<typename S, typename T>
	void writeConnectivity_dim2_degree2(DataCollector<S>* dataCollector, T &file)
	{
	  std::list<ElementInfo> *elements = dataCollector->getElementInfos();
	  std::list<ElementInfo>::iterator elementIt;

	  std::vector< std::vector<int> > *interpPoints = dataCollector->getInterpPoints();
	  std::vector< std::vector<int> >::iterator pointIt;

	  int nVertices = dataCollector->getNumberVertices();

	  for (pointIt = interpPoints->begin(), elementIt = elements->begin(); 
	      pointIt != interpPoints->end(); 
	      ++pointIt, ++elementIt) {
	  
	    file << " " << elementIt->vertexInfo[0]->outputIndex
		<< " " << (*pointIt)[1] + nVertices
		<< " " << (*pointIt)[2] + nVertices << "\n";

	    file << " " << elementIt->vertexInfo[2]->outputIndex 
		<< " " << (*pointIt)[0] + nVertices
		<< " " << (*pointIt)[1] + nVertices << "\n";

	    file << " " << elementIt->vertexInfo[1]->outputIndex
		<< " " << (*pointIt)[0] + nVertices
		<< " " << (*pointIt)[2] + nVertices << "\n";

	    file << " " << (*pointIt)[0] + nVertices
		<< " " << (*pointIt)[1] + nVertices 
		<< " " << (*pointIt)[2] + nVertices << "\n";
	  }
	}


	template<typename S, typename T>
	void writeConnectivity_dim2_degree3(DataCollector<S>* dataCollector, T &file)
	{
	  std::list<ElementInfo> *elements = dataCollector->getElementInfos();
	  std::list<ElementInfo>::iterator elementIt;

	  std::vector< std::vector<int> > *interpPoints = dataCollector->getInterpPoints();
	  std::vector< std::vector<int> >::iterator pointIt;

	  int nVertices = dataCollector->getNumberVertices();

	  for (pointIt = interpPoints->begin(), elementIt = elements->begin(); 
	      pointIt != interpPoints->end(); 
	      ++pointIt, ++elementIt) {
				
	    file << " " << elementIt->vertexInfo[0]->outputIndex
		<< " " << (*pointIt)[3] + nVertices
		<< " " << (*pointIt)[4] + nVertices << "\n";

	    file << " " << (*pointIt)[4] + nVertices
		<< " " << (*pointIt)[5] + nVertices
		<< " " << (*pointIt)[6] + nVertices << "\n";

	    file << " " << (*pointIt)[3] + nVertices
		<< " " << (*pointIt)[4] + nVertices
		<< " " << (*pointIt)[6] + nVertices << "\n";

	    file << " " << (*pointIt)[2] + nVertices
		<< " " << (*pointIt)[3] + nVertices
		<< " " << (*pointIt)[6] + nVertices << "\n";

	    file << " " << elementIt->vertexInfo[1]->outputIndex
		<< " " << (*pointIt)[0] + nVertices
		<< " " << (*pointIt)[5] + nVertices << "\n";

	    file << " " << (*pointIt)[0] + nVertices
		<< " " << (*pointIt)[6] + nVertices
		<< " " << (*pointIt)[5] + nVertices << "\n";

	    file << " " << (*pointIt)[0] + nVertices
		<< " " << (*pointIt)[1] + nVertices
		<< " " << (*pointIt)[6] + nVertices << "\n";

	    file << " " << (*pointIt)[1] + nVertices
		<< " " << (*pointIt)[2] + nVertices
		<< " " << (*pointIt)[6] + nVertices << "\n";

	    file << " " << elementIt->vertexInfo[2]->outputIndex 
		<< " " << (*pointIt)[1] + nVertices
		<< " " << (*pointIt)[2] + nVertices << "\n";
	  
	    
	  }
	}


	template<typename S, typename T>
	void writeConnectivity_dim2_degree4(DataCollector<S>* dataCollector, T &file)
	{
	  std::list<ElementInfo> *elements = dataCollector->getElementInfos();
	  std::list<ElementInfo>::iterator elementIt;

	  std::vector< std::vector<int> > *interpPoints = dataCollector->getInterpPoints();
	  std::vector< std::vector<int> >::iterator pointIt;

	  int nVertices = dataCollector->getNumberVertices();

	  for (pointIt = interpPoints->begin(), elementIt = elements->begin(); 
	      pointIt != interpPoints->end(); 
	      ++pointIt, ++elementIt) {

	    file << " " << elementIt->vertexInfo[0]->outputIndex 
		<< " " << (*pointIt)[5] + nVertices
		<< " " << (*pointIt)[6] + nVertices << "\n";

	    file << " " << (*pointIt)[5] + nVertices
		<< " " << (*pointIt)[9] + nVertices
		<< " " << (*pointIt)[6] + nVertices << "\n";

	    file << " " << (*pointIt)[6] + nVertices
		<< " " << (*pointIt)[7] + nVertices
		<< " " << (*pointIt)[9] + nVertices << "\n";

	    file << " " << (*pointIt)[7] + nVertices
		<< " " << (*pointIt)[9] + nVertices
		<< " " << (*pointIt)[10] + nVertices << "\n";

	    file << " " << (*pointIt)[7] + nVertices
		<< " " << (*pointIt)[8] + nVertices
		<< " " << (*pointIt)[10] + nVertices << "\n";

	    file << " " << (*pointIt)[0] + nVertices
		<< " " << (*pointIt)[8] + nVertices
		<< " " << (*pointIt)[10] + nVertices << "\n";

	    file << " " << elementIt->vertexInfo[1]->outputIndex 
		<< " " << (*pointIt)[0] + nVertices
		<< " " << (*pointIt)[8] + nVertices << "\n";

	    file << " " << (*pointIt)[4] + nVertices
		<< " " << (*pointIt)[5] + nVertices
		<< " " << (*pointIt)[9] + nVertices << "\n";

	    file << " " << (*pointIt)[4] + nVertices
		<< " " << (*pointIt)[9] + nVertices
		<< " " << (*pointIt)[11] + nVertices << "\n";

	    file << " " << (*pointIt)[9] + nVertices
		<< " " << (*pointIt)[10] + nVertices
		<< " " << (*pointIt)[11] + nVertices << "\n";

	    file << " " << (*pointIt)[1] + nVertices
		<< " " << (*pointIt)[10] + nVertices
		<< " " << (*pointIt)[11] + nVertices << "\n";

	    file << " " << (*pointIt)[0] + nVertices
		<< " " << (*pointIt)[1] + nVertices
		<< " " << (*pointIt)[10] + nVertices << "\n";

	    file << " " << (*pointIt)[3] + nVertices
		<< " " << (*pointIt)[4] + nVertices
		<< " " << (*pointIt)[11] + nVertices << "\n";

	    file << " " << (*pointIt)[2] + nVertices
		<< " " << (*pointIt)[3] + nVertices
		<< " " << (*pointIt)[11] + nVertices << "\n";

	    file << " " << (*pointIt)[1] + nVertices
		<< " " << (*pointIt)[2] + nVertices
		<< " " << (*pointIt)[11] + nVertices << "\n";

	    file << " " << elementIt->vertexInfo[2]->outputIndex 
		<< " " << (*pointIt)[2] + nVertices
		<< " " << (*pointIt)[3] + nVertices << "\n";

	  }
	}
	
      } // end namespace detail
    } // end namespace VtkWriter
} } // end namespace io, AMDiS

#endif // AMDIS_VTKWRITER_DETAIL_HH
