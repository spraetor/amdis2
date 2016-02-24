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

#ifndef AMDIS_VALUEREADER_H
#define AMDIS_VALUEREADER_H

#include <cstring>
#include "DOFVector.h"
#include "Mesh.h"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Input
     *
     * \brief
     * Namespace of methods which read a value file in AMDiS format and copies
     * the data to a DOF vector.
     */
    namespace ValueReader
    {

      template<typename Container>
      void readValue(std::string /*filename*/,
                     Mesh* /*mesh*/,
                     Container& /*vec*/,
                     MacroInfo* /*macroFileInfo*/)
      {
        ERROR_EXIT("ValueReader not implemented for this container type!\n");
      }

      /// Copies the values of a value file to a DOF vector.
      void readValue(std::string filename,
                     Mesh* mesh,
                     DOFVector<double>* dofVector,
                     MacroInfo* macroFileInfo);

      inline void readValue(std::string filename,
                     Mesh* mesh,
                     DOFVector<double>& dofVector,
                     MacroInfo* macroFileInfo)
      {
        readValue(filename, mesh, &dofVector, macroFileInfo);
      }

    } // end namespace ValueReader
  }
} // end namespace io, AMDiS

#endif
