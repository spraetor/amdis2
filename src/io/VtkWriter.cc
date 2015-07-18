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

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include <mpi.h>
#endif

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/convenience.hpp>

#include "VtkWriter.h"
#include "DataCollector.h"
#include "DOFVector.h"
#include "SurfaceRegion_ED.h"
#include "ElementRegion_ED.h"
#include "AdaptInfo.h"



namespace AMDiS 
{ 
  namespace io
  {
    namespace VtkWriter
    {
      using namespace std;

      
      void writeFile(DOFVector<double> *values, 
		     string filename,
		     Vtuformat format,
		     bool highPrecision, 
		     bool writeParallel
		    )
      {
	DataCollector<> dc(values->getFeSpace(), values);
	vector<DataCollector<>*> dcList(0);
	dcList.push_back(&dc);
	writeFile(dcList, filename, format, highPrecision, writeParallel);
      }


      void writeFile(vector<DOFVector<double>* > &values,
		     string filename,
		     Vtuformat format,
		     bool highPrecision, 
		     bool writeParallel
		    )
      {
	vector<DataCollector<>*> dcList(0);
	for (size_t i = 0; i < values.size(); i++)
	  dcList.push_back(new DataCollector<>(values[i]->getFeSpace(), values[i]));	
	writeFile(dcList, filename, format, highPrecision, writeParallel);
	for (size_t i = 0; i < values.size(); i++)
	  delete dcList[i];	
      }


      void writeFile(WorldVector<DOFVector<double>* > &values,
		     string filename,
		     Vtuformat format,
		     bool highPrecision, 
		     bool writeParallel
		    )
      {
	vector<DataCollector<>*> dcList(0);
	for (int i = 0; i < values.getSize(); i++)
	  dcList.push_back(new DataCollector<>(values[i]->getFeSpace(), values[i]));	
	writeFile(dcList, filename, format, highPrecision, writeParallel);
	for (int i = 0; i < values.getSize(); i++)
	  delete dcList[i];	
      }


      void writeFile(DOFVector<WorldVector<double> > *values,
		     string filename,
		     Vtuformat format,
		     bool highPrecision, 
		     bool writeParallel
		    )
      {
	WorldVector<DOFVector<double>*> valuesWV;
	for (int i =0 ; i < valuesWV.getSize(); i++)
	  valuesWV[i] = new DOFVector<double>(values->getFeSpace(), 
					      "values["+to_string(i)+"]");
	transform(values, &valuesWV);
	writeFile(valuesWV, filename, format, highPrecision, writeParallel);
	for (int i = 0; i < valuesWV.getSize(); i++)
	  delete valuesWV[i];
      }

      
      void writeFile(SystemVector *values, 
		     string filename,
		     Vtuformat format,
		     bool highPrecision, 
		     bool writeParallel
		    )
      {
	vector<DataCollector<>*> dcList(0);
	for (int i = 0; i < values->getSize(); i++)
	  dcList.push_back(new DataCollector<>(values->getDOFVector(i)->getFeSpace(), 
					    values->getDOFVector(i)));    
	writeFile(dcList, filename, format, highPrecision, writeParallel);
	for (size_t i = 0; i < dcList.size(); i++)
	  delete dcList[i];    
      }
      
      
      void writeFile(vector<DataCollector<>*> &dcList,
		     string filename,
		     Vtuformat format,
		     bool highPrecision, 
		     bool writeParallel
		    )
      {
	vector<string> componentNames;
	for (size_t i = 0; i < dcList.size(); i++)
	  componentNames.push_back(dcList[i]->getValues()->getName());
	::AMDiS::io::VtkWriter::Aux writer(&dcList, componentNames, format, highPrecision);
	
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
	if (writeParallel) {
	  int sPos = filename.find(".vtu");
	  TEST_EXIT(sPos >= 0)("Failed to find file postfix!\n");
	  string name = filename.substr(0, sPos);
	  
	  if (MPI::COMM_WORLD.Get_rank() == 0) {

	    detail::writeParallelFile(name + ".pvtu", MPI::COMM_WORLD.Get_size(), 
				      name, ".vtu", componentNames, format, highPrecision);
	  }
	  
	  filename = name + "-p" + std::to_string(MPI::COMM_WORLD.Get_rank()) + "-.vtu";
	}
#endif
	writer.writeFile(filename);
      }

    } // end namespace VtkWriter
  } // end namespace io
} // end namespace AMDiS
