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


/** \file VtkReader.h */

#ifndef AMDIS_VTKREADER_DETAIL_H
#define AMDIS_VTKREADER_DETAIL_H

// need some extension-methods/libs (pugixml, nanoflann)
#ifdef HAVE_EXTENSIONS

#include <cstring>
#include <boost/archive/iterators/binary_from_base64.hpp>
#include <boost/archive/iterators/transform_width.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#ifdef HAVE_COMPRESSION
#include <boost/iostreams/filter/zlib.hpp>
#endif
#include "DOFVector.h"
#include "Mesh.h"
#include "boost/filesystem.hpp"
#include "pugixml.hpp"
// extensions {
#include "kdtree_nanoflann.h"
#include "VectorOperations.h"
// }

#define BLOCKSIZE 32768

namespace AMDiS 
{
  namespace io 
  {
    namespace VtkReader 
    {
      namespace detail
      {
	typedef boost::archive::iterators::transform_width<
	          boost::archive::iterators::binary_from_base64<
	              std::string::const_iterator
	          >, 8, 6 
	        > base64_binary;
	
	inline void valueVector2type(std::vector<double> p, 
				    double& value)
	{
	  TEST_EXIT(p.size() != 0)("Not enough data for assignment!\n");
	  value = p[0];
	}

	
	inline void valueVector2type(std::vector<double> p, 
				    WorldVector<double>& value)
	{
	  TEST_EXIT(static_cast<int>(p.size()) == Global::getGeo(WORLD))("Not enough data for assignment!\n");
	  for (int i = 0; i < Global::getGeo(WORLD); i++)
	    value[i] = p[i];
	}

	
	inline void valueVector2type(std::vector<double> p, 
				    std::vector<double>& value)
	{
	  TEST_EXIT(p.size() != 0)("Not enough data for assignment!\n");
	  for (size_t i = 0; i < p.size(); i++)
	    value.push_back(p[i]);
	}
	
	inline std::string base64ToStr(std::string base64)
	{
	  using namespace std;
	  unsigned int paddChars = count(base64.begin(), base64.end(), '=');
	  std::replace(base64.begin(),base64.end(),'=','A'); 
	  string result(base64_binary(base64.begin()), base64_binary(base64.end()));
	  result.erase(result.end()-paddChars,result.end());
	  return result;
	}

#ifdef HAVE_COMPRESSION
	inline std::string decompress(std::string text)
	{
	  std::stringstream tmp1, tmp2;
	  tmp1.str(text);
	  boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
	  in.push(boost::iostreams::zlib_decompressor());
	  in.push(tmp1);
	  boost::iostreams::copy(in, tmp2);
	  return tmp2.str();
	}
#endif

	inline std::string getInnerDataArray(std::string& input, bool zlib)
	{
	  FUNCNAME("VtkReader::detail::getInnerDataArray()");
	  
	  using namespace std;
	  
	  char* tmp = NULL;
	  int* ptr = NULL;
	  int nBytes = 0;
	  string header = "", body = "", data = "";
	  
	  if(zlib) {
#ifdef HAVE_COMPRESSION
	    string s = input.substr(0, 8);
	    s = detail::base64ToStr(s);
	    tmp =  const_cast<char*>(s.c_str());
	    ptr = reinterpret_cast<int*>(tmp);
	    
	    int nBlocks = *ptr;
	    int headerSize = (((4 * nBlocks + 12) % 3) > 0) ? 4 * ((4 * nBlocks + 12) / 3 + 1) : 4 * ((4 * nBlocks + 12) / 3);
	    header = input.substr(0, headerSize);
	    body = input.substr(headerSize);
	    header = detail::base64ToStr(header);
	    body = detail::base64ToStr(body);
	    
	    int blockSize, finalSize, offset = 0;
	    tmp =  const_cast<char*>(header.c_str());
	    ptr = reinterpret_cast<int*>(tmp);
	    ptr++;
	    blockSize = *ptr++;
	    finalSize = *ptr++;
	    nBytes = (nBlocks - 1) * blockSize + finalSize;
	    
	    for(int i = 0; i < nBlocks; i++, ptr++) {
	      string subBody = body.substr(offset, *ptr);
	      data += decompress(subBody);
	      offset += *ptr;
	    }
#else
            ERROR_EXIT("HAVE_COMPRESSION OFF. VtkReader cannot read APPENDED_COMPRESSED vtu files.\n");
#endif
	  } else {
	    header = detail::base64ToStr(input);
	    tmp =  const_cast<char*>(header.c_str());
	    ptr = reinterpret_cast<int*>(tmp);
	    nBytes = *ptr;
	    data = header.substr(4);
	  }
	  TEST_EXIT((unsigned)nBytes == data.length())("Compressed DataArray is wrong. nBytes: %i and data length: %i.\n", nBytes, data.length());
	  return data;
	}
	
	inline void binary2pointList(std::string& input,
				     std::string type,
			             bool zlib,
				     std::vector<WorldVector<double> >& pointList)
	{
	  int dow = Global::getGeo(WORLD);
	  int dowExplicit = 3;
	  
	  std::string inner = getInnerDataArray(input, zlib);
	  int nBytes = inner.length();
	  char* tmp = const_cast<char*>(inner.c_str());
	  
	  pointList.clear();

	  if(type == "Float32") {
	    int size = nBytes / 4;
	    float* ptr = reinterpret_cast<float*>(tmp); 
	    
	    for(int i = 0; i < size;) {
	      WorldVector<double> p;
	      for (int j = 0; j < dowExplicit && i < size; i++, j++, ptr++)
		if (j < dow)
		  p[j] = boost::lexical_cast<double>(*ptr);
	      pointList.push_back(p);
	    }
	  }
	  else if(type == "Float64") {
	    int size = nBytes / 8;
	    double* ptr = reinterpret_cast<double*>(tmp); 
	    
	    for(int i = 0; i < size;) {
	      WorldVector<double> p;
	      for (int j = 0; j < dowExplicit && i < size; i++, j++, ptr++)
		if (j < dow)
		  p[j] = *ptr;
	      pointList.push_back(p);
	    }
	  } else {
	    ERROR_EXIT("Currently only supports Float32 and Float64 data type.\n");
	  }
	}

	
	inline void string2pointList(std::string& input, 
				    std::vector<WorldVector<double> >& pointList)
	{
	  int dow = Global::getGeo(WORLD);
	  int dowExplicit = 3;

	  std::vector<std::string> valueStringList;
	  std::stringstream ss(input);
	  std::string buf;
	  while (ss >> buf)
	    valueStringList.push_back(buf);

	  pointList.clear();

	  std::vector<std::string>::iterator it;
	  size_t j = 0;
	  for (it = valueStringList.begin(); it != valueStringList.end();j++) {
	    WorldVector<double> p;
	    for (int i = 0; i < dowExplicit && it != valueStringList.end(); i++, it++)
	      if (i < dow)
		p[i] = boost::lexical_cast<double>(*it);
	    pointList.push_back(p);
	  }
	}

	template<typename T>
	void binary2valueList(std::string& input, 
			            std::string type,
				    bool zlib,
				    std::vector<T>& valueList, 
				    int numComponent = 1, 
				    int numComponentMax = -1)
	{
	  if (numComponentMax < 0)
	    numComponentMax = numComponent;
	  
	  std::string inner = getInnerDataArray(input, zlib);
	  int nBytes = inner.length();
	  char* tmp = const_cast<char*>(inner.c_str());
	  
	  valueList.clear();
	  
	  if(type == "Float32") {
	    int size = nBytes / 4;
	    float* ptr = reinterpret_cast<float*>(tmp); 
	    
	    for(int i = 0; i < size;) {
	      std::vector<double> p;
	      for (int j = 0; j < numComponentMax && i < size; i++, j++, ptr++) 
		if (j < numComponent)
   		  p.push_back(boost::lexical_cast<double>(*ptr));
		
	      T value; valueVector2type(p, value);
	      valueList.push_back(value);
	    }
	  } else if(type == "Float64") {
	    int size = nBytes / 8;
	    double* ptr = reinterpret_cast<double*>(tmp); 
	    
	    for(int i = 0; i < size;) {
	      std::vector<double> p;
	      for (int j = 0; j < numComponentMax && i < size; i++, j++, ptr++) 
		if (j < numComponent) 
		  p.push_back(*ptr);
	      T value; valueVector2type(p, value);
	      valueList.push_back(value);
	    }
	  }
	  else {
	    ERROR_EXIT("Currently only supports Float32 and Float64 data type.\n");
	  }
	}
	
	template<typename T>
	void string2valueList(std::string& input, 
				    std::vector<T>& valueList, 
				    int numComponent = 1, 
				    int numComponentMax = -1)
	{
	  if (numComponentMax < 0)
	    numComponentMax = numComponent;

	  std::vector<std::string> valueStringList;
	  std::stringstream ss(input);
	  std::string buf;
	  while (ss >> buf)
	    valueStringList.push_back(buf);

	  valueList.clear();

	  std::vector<std::string>::iterator it;
	  for (it = valueStringList.begin(); it != valueStringList.end();) {
	    std::vector<double> p;
	    for (int i = 0; i < numComponentMax && it != valueStringList.end(); i++, it++)
	      if (i < numComponent)
		p.push_back(boost::lexical_cast<double>(*it));
	      
	    T value; valueVector2type(p, value);
	    valueList.push_back(value);
	  }
	}


	// find point in mesh using KD-tree structure
	inline size_t getNearestIndex(extensions::KD_Tree& tree, 
				      WorldVector<double>& x)
	{
	  // do a knn search
	  const size_t num_nnp = 1;
	  std::vector<size_t> ret_indexes(num_nnp);
	  std::vector<double> out_dists_sqr(num_nnp);

	  nanoflann::KNNResultSet<double> resultSet(num_nnp);
	  resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );

	  tree.index->findNeighbors(resultSet, x.begin(), nanoflann::SearchParams(10));
	  return ret_indexes[0];
	}

      } // end namespace detail    
    } // end namespace VtkReader
  } // end namespace io
} // end namespace AMDiS

#endif // HAVE_EXTENSIONS

#endif // AMDIS_VTKREADER_DETAIL_H
