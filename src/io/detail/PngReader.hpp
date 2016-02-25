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



/** \file PngReader.h */

#ifndef AMDIS_PNGREADER_DETAIL_H
#define AMDIS_PNGREADER_DETAIL_H

#if defined HAVE_PNG

#include "ElInfo.hpp"
#include "FixVec.hpp"
#include "Mesh.hpp"
#include "Traverse.hpp"

namespace AMDiS
{
  namespace io
  {

    /** \ingroup Input
     *
     * \brief
     * Static class which reads a png file and gets values for each pixel
     */
    namespace PngReader
    {
      namespace detail
      {
        inline
        void getMeshDimension(Mesh* mesh, double& xMin, double& xMax, double& yMin, double& yMax)
        {
          WorldVector<double> minDim;
          minDim.set(1.e10);
          WorldVector<double> maxDim;
          maxDim.set(-1.e10);

          TraverseStack stack;
          ElInfo* elInfo = stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL | Mesh::FILL_COORDS);
          while (elInfo)
          {
            for (int i = 0; i <= mesh->getDim(); i++)
            {
              WorldVector<double>& coords = elInfo->getMacroElement()->getCoord(i);
              for (int j = 0; j < coords.getSize(); ++j)
              {
                minDim[j] = std::min(minDim[j], coords[j]);
                maxDim[j] = std::max(maxDim[j], coords[j]);
              }
            }
            elInfo = stack.traverseNext(elInfo);
          }

          xMin = minDim[0];
          xMax = maxDim[0];
          yMin = minDim[1];
          yMax = maxDim[1];
        }

      } // end namespace detail
    } // end namespace PngReader
  }
} // end namespace io, AMDiS

#endif

#endif // AMDIS_PNGREADER_DETAIL_H