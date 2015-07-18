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


#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/lexical_cast.hpp>

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

#include "VtkVectorWriter.h"
#include "DataCollector.h"
#include "DOFVector.h"
// #include "SurfaceRegion_ED.h"
// #include "ElementRegion_ED.h"

namespace AMDiS 
{
  namespace io
  {
    namespace VtkVectorWriter
    {
      
      void writeFile(std::vector<DOFVector<double>* > &values,
		     std::string filename,
		     bool writeParallel,
		     bool writeAs3dVector)
      {
	DOFVector<std::vector<double> > *newValues = new DOFVector<std::vector<double> >(values[0]->getFeSpace(), values[0]->getName());
	std::vector<DOFIterator<double>* > iterators;
	for (size_t i = 0; i < values.size(); i++)
	  iterators.push_back(new DOFIterator<double>(values[i],USED_DOFS));
	for (size_t i = 0; i < iterators.size(); i++)
	  iterators[i]->reset();
	DOFIterator<std::vector<double> > resultIter(newValues, USED_DOFS);

	for(resultIter.reset(); !resultIter.end(); resultIter++)
	{
	  std::vector<double> val(0);
	  for (size_t i = 0; i < static_cast<size_t>(iterators.size()); i++)
	    val.push_back(*(*(iterators[i])));

	  *resultIter = val;
	  
	  for (size_t i = 0; i < static_cast<size_t>(iterators.size()); i++)
	    (*(iterators[i]))++;
	}

	writeFile(newValues, filename, writeParallel, writeAs3dVector);
	for (size_t i = 0; i < iterators.size(); i++)
	  delete iterators[i];
	delete newValues;
      }


      void writeFile(WorldVector<DOFVector<double>* > &values,
		     std::string filename,
		     bool writeParallel,
		     bool writeAs3dVector)
      {
	DOFVector<WorldVector<double> > *newValues = new DOFVector<WorldVector<double> >(values[0]->getFeSpace(), values[0]->getName());
	WorldVector<DOFIterator<double>* > iterators;
	for (int i = 0; i < values.getSize(); i++)
	  iterators[i] = new DOFIterator<double>(values[i],USED_DOFS);
	for (int i = 0; i < iterators.getSize(); i++)
	  iterators[i]->reset();
	DOFIterator<WorldVector<double> > resultIter(newValues, USED_DOFS);

	for(resultIter.reset(); !resultIter.end(); resultIter++)
	{
	  for (int i = 0; i < iterators.getSize(); i++)
	    (*resultIter)[i] = *(*(iterators[i]));

	  for (int i = 0; i < iterators.getSize(); i++)
	    (*(iterators[i]))++;
	}

	writeFile(newValues, filename, writeParallel, writeAs3dVector);
	for (size_t i = 0; i < static_cast<size_t>(iterators.getSize()); i++)
	  delete iterators[i];
	delete newValues;
      }


      void writeFile(SystemVector *values,
		     std::string filename,
		     bool writeParallel,
		     bool writeAs3dVector)
      {
	DOFVector<std::vector<double> > *newValues = new DOFVector<std::vector<double> >(values->getDOFVector(0)->getFeSpace(), values->getName());
	std::vector<DOFIterator<double>* > iterators;
	for (size_t i = 0; i < static_cast<size_t>(values->getSize()); i++)
	  iterators.push_back(new DOFIterator<double>(values->getDOFVector(i),USED_DOFS));
	for (size_t i = 0; i < iterators.size(); i++)
	  iterators[i]->reset();
	DOFIterator<std::vector<double> > resultIter(newValues, USED_DOFS);

	for(resultIter.reset(); !resultIter.end(); resultIter++)
	{
	  std::vector<double> val(0);
	  for (size_t i = 0; i < static_cast<size_t>(iterators.size()); i++)
	    val.push_back(*(*(iterators[i])));

	  *resultIter = val;

	  for (size_t i = 0; i < static_cast<size_t>(iterators.size()); i++)
	    (*(iterators[i]))++;
	}

	writeFile(newValues, filename, writeParallel, writeAs3dVector);
	for (size_t i = 0; i < iterators.size(); i++)
	  delete iterators[i];
	delete newValues;
      }
      
    } // end namespace VtkVectorWriter
  } // end namespace io
} // end namespace AMDiS
