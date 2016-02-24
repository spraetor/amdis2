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


#include <list>
#include <vector>
#include <fstream>
#include <algorithm>
#include <cstring>

#include "ValueWriter.h"
#include "DOFVector.h"
#include "DOFAdmin.h"
#include "ElInfo.h"
#include "BasisFunction.h"
#include "DataCollector.h"

namespace AMDiS
{
  namespace io
  {

    namespace ValueWriter
    {
      using namespace std;

      void writeValues(DataCollector<>* dc, std::string filename,
                       double time, int level,
                       Flag traverseFlag,
                       bool (*writeElem)(ElInfo*))
      {
        FUNCNAME("writeValues()");

        TEST_EXIT(dc)("no data collector\n");

        ofstream file(filename.c_str(), ios_base::out);

        file << "mesh name: " << dc->getMesh()->getName() << endl << endl;
        file << "time: " << std::scientific << time << endl << endl;
        file << "number of values: 1" << endl << endl;
        file << "value description: " << dc->getValues()->getName() << endl;
        file << "number of interpolation points: " << dc->getNumberInterpPoints()
             << endl;
        file << "type: scalar" << endl;
        file << "interpolation type: lagrange" << endl;
        file << "interpolation degree: " << dc->getFeSpace()->getBasisFcts()->getDegree()
             << endl;
        file << "end of description: " << dc->getValues()->getName()
             << endl << endl;

        /* ----- write vertex values -----*/
        DOFVector<int>::Iterator intPointIt(dc->getInterpPointInd(), USED_DOFS);
        DOFVector<double>::Iterator valueIt(dc->getValues(), USED_DOFS);
        DOFVector<list<WorldVector<double>>>::Iterator
        coordIt(dc->getDofCoords(), USED_DOFS);


        file << "vertex values: " << dc->getValues()->getName() << endl;

        for (intPointIt.reset(), valueIt.reset(), coordIt.reset();
             !intPointIt.end();
             ++intPointIt, ++valueIt, ++coordIt)
        {

          if (*intPointIt == -2)
          {
            for (int i = 0; i < (int) coordIt->size(); i++)
            {
              file << std::scientific << *valueIt << endl;
            }
          }
        }

        file << endl << endl;

        /* ----- write interpolation values ----- */
        file << "interpolation values: " << dc->getValues()->getName() << endl;

        for (intPointIt.reset(), valueIt.reset();
             !intPointIt.end();
             ++intPointIt, ++valueIt)
        {

          if (*intPointIt >= 0)
          {
            file << std::scientific << *valueIt << endl;
          }
        }

        file << endl << endl;

        /* ----- write interpolation points for each simplex */
        file << "element interpolation points: " << dc->getValues()->getName() << endl;

        vector<vector<int>>* interpPoints = dc->getInterpPoints();
        vector<vector<int>>::iterator it1;
        vector<int>::iterator it2;

        for (it1 = interpPoints->begin(); it1 != interpPoints->end(); ++it1)
        {
          for (it2 = it1->begin(); it2 != it1->end(); ++it2)
          {
            file << (*it2) << " ";
          }
          file << endl;
        }

        file << endl;

        file.close();
      }

    } // end namespace ValueWriter
  }
} // end namespace io, AMDiS
