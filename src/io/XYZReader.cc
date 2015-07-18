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

#ifdef HAVE_EXTENSIONS

#include <cstring>
#include <fstream>
#include <vector>

#include "ValueReader.h"
#include "MacroInfo.h"
#include "detail/VtkReader.h"

namespace AMDiS { namespace io {

  namespace XYZReader
  {
    using namespace std;
      
    /** \brief
    * Copies the values of a value file to a DOF vector.
    *
    * The DOF vector must have been created by a corresponding mesh file. The
    * information about this file and the macro triangulation are stored in
    * macroFileInfo. The function now reads the corresponding value file and
    * copies the values to the correct positions in the DOF vector.
    */
    void readFile(std::string filename,
		  std::pair<std::vector<WorldVector<double> >, 
			    std::vector<std::vector<double> > > &data)
    {
      FUNCNAME("readFile()");
      using VtkReader::detail::string2valueList;

      TEST_EXIT(filename != "")("Filename not specified!\n");

      string line;

      // open the file and read the data to the vector values.
      ifstream file(filename.c_str(), ios_base::in);
      getline(file, line);
      int numRows = boost::lexical_cast<int>(line);
      if (numRows == 0 || file.eof())
	return;
      
      getline(file, line); // comments
      
      TEST_EXIT(!file.eof())("Not enough values stored in file!\n");
      getline(file, line); // first line
      vector<double> values_tmp;
      string2valueList(line, values_tmp);
      
      int numValues = values_tmp.size() - 1 - 3; // 1st value=nr, 2-4th value=coordinates

      int nr;
      WorldVector<double> p;
      vector<double> values(numValues);
      double tmp;
      for (int i = 0; i < Global::getGeo(WORLD); i++)
	p[i] = values_tmp[i+1];
      for (int i = 0; i < numValues; i++)
	values[i] = values_tmp[i+4];
      data.first.push_back(p);
      data.second.push_back(values);

      while (!file.eof()) {
	file >> nr;
	for (int i = 0; i < Global::getGeo(WORLD); i++)
	  file >> p[i];
	for (int i = Global::getGeo(WORLD); i < 3; i++)
	  file >> tmp;	
	for (int i = 0; i < numValues; i++)
	  file >> values[i];
	if (!file.good())
	  break;
	data.first.push_back(p);
	data.second.push_back(values);
      }

      file.close();
    }
    
  } // end namespace XYZReader
} } // end namespace io, AMDiS

#endif


