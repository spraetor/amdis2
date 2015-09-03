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



/** \file ValueReader.h */

#ifndef AMDIS_PNGREADER_H
#define AMDIS_PNGREADER_H

#if defined HAVE_PNG

#include <cstring>

#include "DOFVector.h"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Input
     *
     * \brief
     * Namespace which provieds readers to read a png file and gets values for each pixel
     */
    namespace PngReader
    {

      /// Interface for general containers not implemented. Specializations below.
      template<typename Container>
      void readFile(std::string filename, Container& vec)
      {
        ERROR_EXIT("PngReader not implemented for this container type!\n");
      }

      /// Copies the values of a png file to a DOF vector. Using container pointer.
      void readFile(std::string filename, DOFVector<double>* dofVector);


      /// Copies the values of a png file to a DOF vector. Using container reference.
      inline
      void readFile(std::string filename, DOFVector<double>& dofVector)
      {
        readFile(filename, &dofVector);
      }

    }

  }
} // end namespace io, AMDiS

#endif // HAVE_PNG

#endif // AMDIS_PNGREADER_H
