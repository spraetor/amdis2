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


#include <algorithm>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "Spreadsheet.h"
#include "Global.h"

namespace AMDiS { namespace io {
  
  using namespace std;
  using namespace boost;

  using boost::lexical_cast;
  
  void Spreadsheet::write(string filename)
  {
    int nRows = static_cast<int>(data.size());
    int maxDataLength = 0;
    
    for (int i = 0; i < nRows; i++)
      maxDataLength = std::max(maxDataLength, static_cast<int>(data[i].size()));

    for (int i = 0; i < nRows; i++)
      data[i].resize(maxDataLength, 0.0);

    ofstream file;
    file.open(filename.c_str());

    for (int i = 0; i < nRows; i++) {
      for (int j = 0; j < maxDataLength; j++) {
	file << data[i][j];
	if (j + 1 < maxDataLength)
	  file << " ";
	else 
	  file << "\n";
      }
    }

    file.close();
  }


  void Spreadsheet::read(string filename)
  {
    data.clear();

    string line;
    vector<string> lineSplit;

    ifstream file;
    file.open(filename.c_str());

    while (!file.eof()) {
      getline(file, line);
      split(lineSplit, line, is_any_of(" ,"));

      if (lineSplit.size() == 0)
	continue;

      if ((lineSplit.size() == 1) && (lineSplit[0] == ""))
	continue;

      if (lineSplit[0].find("#") == 0)
	continue;

      vector<double> lineData(lineSplit.size());
      for (unsigned int i = 0; i < lineSplit.size(); i++)
	lineData[i] = lexical_cast<double>(lineSplit[i]);

      data.push_back(lineData);
    }

    file.close();
  }

} } // end namespace io, AMDiS
