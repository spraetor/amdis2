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



/** \file Serializer.h */

#ifndef AMDIS_SERIALIZER_H
#define AMDIS_SERIALIZER_H

#if HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

#include <map>
#include <boost/lexical_cast.hpp>

#include "Global.h"
#include "Initfile.h"
#include "AdaptInfo.h"
#include "io/FileWriterInterface.h"

namespace AMDiS {

  template<typename ProblemType>
  class Serializer : public FileWriterInterface
  {
  public:
    Serializer(ProblemType *prob) 
      : name(""), 
	problem(prob),
	tsModulo(1), 
	appendIndex(0),
	indexLength(5),
	indexDecimals(3),
	timestepNumber(-1)
    {
      FUNCNAME("Serializer::Serializer()");

      Parameters::get(problem->getName() + "->output->serialization filename", 
		      name);
      Parameters::get(problem->getName() + "->output->write every i-th timestep", 
		      tsModulo);
      TEST_EXIT(name != "")("No filename!\n");

      Parameters::get(problem->getName() + "->output->append serialization index", appendIndex);
      Parameters::get(problem->getName() + "->output->index length", indexLength);
      Parameters::get(problem->getName() + "->output->index decimals", indexDecimals);
    }


    Serializer(ProblemType *prob, std::string filename, int writeEveryIth)
      : name(filename),
	problem(prob),
	tsModulo(writeEveryIth),
	appendIndex(0),
	indexLength(5),
	indexDecimals(3),
	timestepNumber(-1)
    {
      FUNCNAME("Serializer::Serializer()");

      TEST_EXIT(name != "")("No filename!\n");

      Parameters::get(problem->getName() + "->output->append serialization index", appendIndex);
      Parameters::get(problem->getName() + "->output->index length", indexLength);
      Parameters::get(problem->getName() + "->output->index decimals", indexDecimals);
    }


    virtual ~Serializer() {}

    virtual void writeFiles(AdaptInfo *adaptInfo, 
			    bool force,
			    int level = -1,
			    Flag traverseFlag = Mesh::CALL_LEAF_EL,
			    bool (*writeElem)(ElInfo*) = NULL) 
    {
      FUNCNAME("Serializer::writeFiles()");

      TEST_EXIT(tsModulo > 0)
	("Parameter 'write every ith timestep' must be larger than zero!\n");

      timestepNumber++;
      timestepNumber %= tsModulo;
      if ((timestepNumber != 0) && !force)
	return;

      TEST_EXIT(adaptInfo)("No AdaptInfo\n");

      std::string filename = name;
      if (appendIndex) {
	TEST_EXIT(indexLength <= 99)("index lenght > 99\n");
	TEST_EXIT(indexDecimals <= 97)("index decimals > 97\n");
	TEST_EXIT(indexDecimals < indexLength)("index length <= index decimals\n");
	
	char formatStr[9];
	char timeStr[20];
	
	std::sprintf(formatStr, "%%0%d.%df", indexLength, indexDecimals);
	std::sprintf(timeStr, formatStr, adaptInfo ? adaptInfo->getTime() : 0.0);
	
	filename += std::string(timeStr);
      }

#if HAVE_PARALLEL_DOMAIN_AMDIS
      filename += ".p" + boost::lexical_cast<std::string>(MPI::COMM_WORLD.Get_rank());
#endif

      std::ofstream out(filename.c_str());
      TEST_EXIT(out.is_open())("Cannot open serialization file!\n");
      out.write(reinterpret_cast<const char*>(&amdisRevisionNumber), sizeof(int));
      problem->serialize(out);
      adaptInfo->serialize(out);
      out.close();

      MSG("Problem serialized to %s \n", filename.c_str());
    }

  protected:
    /// Name of file to which the problem is serialized.
    std::string name;

    /// Pointer to the problem.
    ProblemType *problem;

    /// The problem is serialized every tsModulo-th timestep.
    int tsModulo;

    /// 0: Don't append time index to filename prefix.
    /// 1: Append time index to filename prefix.
    int appendIndex;

    /// Total length of appended time index.
    int indexLength;

    /// Number of decimals in time index.
    int indexDecimals;

    /// Current timestep number.
    int timestepNumber;
  };

  namespace SerUtil {

    template<typename T>
    void serialize(std::ostream& out, T& data)
    {
      out.write(reinterpret_cast<const char*>(&data), sizeof(T));
    }   

    template<typename T>
    void deserialize(std::istream& in, T& data)
    {
      in.read(reinterpret_cast<char*>(&data), sizeof(T));
    }   



    void serialize(std::ostream& out, DofEdge& data);

    void deserialize(std::istream& in, DofEdge& data);



    void serialize(std::ostream& out, DofFace& data);

    void deserialize(std::istream& in, DofFace& data);



    template<typename T, typename U>
    void serialize(std::ostream& out, std::pair<T, U>& data)
    {
      serialize(out, data.first);
      serialize(out, data.second);
    }

    template<typename T, typename U>
    void deserialize(std::istream& in, std::pair<T, U>& data)
    {
      deserialize(in, data.first);
      deserialize(in, data.second);
    }



    template<typename T>
    void serialize(std::ostream& out, std::vector<T>& data)
    {
      int vecSize = data.size();
      serialize(out, vecSize);
      for (typename std::vector<T>::iterator it = data.begin(); 
	   it != data.end(); ++it) {
	T v = *it;
	serialize(out, v);
      }
    }

    template<typename T>
    void deserialize(std::istream& in, std::vector<T>& data)
    {
      data.clear();

      int vecSize = 0;
      deserialize(in, vecSize);
      data.resize(vecSize);

      for (int i = 0; i < vecSize; i++) {
	T v;
	deserialize(in, v);
	data[i] = v;
      }
    }



    template<typename T>
    void serialize(std::ostream& out, std::set<T>& data)
    {
      int setSize = data.size();
      serialize(out, setSize);
      for (typename std::set<T>::iterator it = data.begin(); 
	   it != data.end(); ++it) {
	T v = *it;
	serialize(out, v);
      }
    }

    template<typename T>
    void deserialize(std::istream& in, std::set<T>& data)
    {
      data.clear();

      int setSize = 0;
      deserialize(in, setSize);

      for (int i = 0; i < setSize; i++) {
	T v;
	deserialize(in, v);
	data.insert(v);
      }
    }



    template<typename T1, typename T2>
    void serialize(std::ostream& out, std::map<T1, T2>& data)
    {
      int mapSize = data.size();
      serialize(out, mapSize);

      for (typename std::map<T1,T2>::iterator it = data.begin(); 
	   it != data.end(); ++it) {
	T1 v1 = it->first;
	T2 v2 = it->second;
	serialize(out, v1);
	serialize(out, v2);
      }
    }

    template<typename T1, typename T2>
    void deserialize(std::istream& in, std::map<T1, T2>& data)
    {
      data.clear();

      int mapSize = 0;
      deserialize(in, mapSize);

      for (int i = 0; i < mapSize; i++) {
	T1 v1;
	T2 v2;
	deserialize(in, v1);
	deserialize(in, v2);
	data[v1] = v2;
      }
    }

  }

}
#endif
