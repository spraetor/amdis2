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


#include <cstring>
#include <fstream>
#include <vector>

#include "ValueReader.h"
#include "MacroInfo.h"

namespace AMDiS
{
  namespace io
  {

    namespace ValueReader
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
      void readValue(string filename, Mesh* mesh,
                     DOFVector<double>* dofVector,
                     MacroInfo* macroFileInfo)
      {
        FUNCNAME("readValue()");

        TEST_EXIT(filename != "")("Filename not specified!\n");
        TEST_EXIT(mesh)("no mesh specified\n");
        TEST_EXIT(dofVector)("no DOF vector specified\n");
        TEST_EXIT(macroFileInfo)("no MacroInfo specified\n");

        string line;

        // open the file and read the data to the vector values.
        ifstream file(filename.c_str(), ios_base::in);

        while (!file.eof())
        {
          getline(file, line);
          if (line.find("vertex values:") != string::npos)
            break;
        }
        // A vertex value cannot be the last line of the file, if the file is correct.
        TEST_EXIT(!file.eof())("Value file does not contain vertex values\n");

        double value;
        vector<double> values;
        values.clear();

        while (!file.eof())
        {
          file >> value;
          if (!file.good())
            break;
          values.push_back(value);
        }

        file.close();

        // Using the macroFileInfo structure, copy the values to the
        // correct positions in the DOF vector.
        for (int i = 0; i < mesh->getNumberOfMacros(); i++)
        {
          for (int j = 0; j < mesh->getGeo(VERTEX); j++)
          {
            int fileIndex = macroFileInfo->mel_vertex[i][j];
            DegreeOfFreedom dofIndex = *(macroFileInfo->dof[fileIndex]);
            (*dofVector)[dofIndex] = values[fileIndex];
          }
        }
      }

    } // end namespace ValueReader
  }
} // end namespace io, AMDiS
