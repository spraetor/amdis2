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



/** \file XYZReader.h */

#ifndef AMDIS_XYZREADER_H
#define AMDIS_XYZREADER_H

#ifdef HAVE_EXTENSIONS

#include <cstring>

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Input
     *
     * \brief
     * Namespace of methods which read a XYZ-file
     */
    namespace XYZReader
    {
      template<typename Container>
      void readFile(std::string filename, Container& data)
      {
        ERROR_EXIT("Can not read xyz-file to given container!\n");
      }

      /// Copies the values of a value file to a DOF vector.
      void readFile(std::string filename,
                    std::pair<std::vector<WorldVector<double>>,
                    std::vector<std::vector<double>>>& data);

    } // end namespace XYZReader
  }
} // end namespace io, AMDiS

#endif

#endif

