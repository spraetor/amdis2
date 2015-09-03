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



/** \file Spreadsheet.h */

#ifndef AMDIS_SPREADSHEET_H
#define AMDIS_SPREADSHEET_H

#include <vector>
#include <string>

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Input
     *  \ingroup Output
     *
     * \brief
     * Implements basic support to read and write tables of data (all are
     * assumed to be of type double. The files are written in ASCII mode and
     * allow to use comments (lines that start with #).
     */
    class Spreadsheet
    {
    public:
      typedef std::vector<std::vector<double>> DataType;

      Spreadsheet()
        : data(0)
      {}

      void addData(std::vector<double>& d)
      {
        data.push_back(d);
      }

      void addData(DataType& d)
      {
        data.insert(data.end(), d.begin(), d.end());
      }

      void setData(DataType& d)
      {
        data = d;
      }

      DataType& getData()
      {
        return data;
      }

      /// writes a spreadsheet file
      void write(std::string filename);

      /// reads from a spreadsheet file
      void read(std::string filename);

    protected:
      DataType data;
    };

  }
} // end namespace io, AMDiS

#endif
