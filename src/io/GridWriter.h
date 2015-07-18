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



/** \file GridWriter.h */

#ifndef AMDIS_GRID_WRITER_H
#define AMDIS_GRID_WRITER_H

#include "Global.h"
#include "AMDiS_fwd.h"

namespace AMDiS { namespace io {

    /** \ingroup Output
    *  \brief
    * Produces a grid-based output of a dof-vector
    */
    namespace GridWriter {
      
      /** \brief
      * Produces a grid-based output of a dof-vector or vector of dof-vectors 
      * \param p array of world coordinates. 
      * - p[0] defines origin of grid-system.
      * - p[1] defines first basis vector by p[1] relatively to p[0].
      * - p[2] defines second basis vector by p[2] relatively to p[0]. (dim >= 2).
      * - p[3] defines third basis vector by p[3] relatively to p[0]. (dim == 3).
      * \param numPoints number of grid points for the different dirctions
      * - numPoints[0] number of grid points in the first direction
      * - numPoints[1] number of grid points in the second direction
      * - numPoints[2] number of grid points in the third direction
      * \param dist distance between two points in the different directions
      * - dist[0] distance in the first direction
      * - dist[1] distance in the second direction
      * - dist[2] distance in the third direction
      * \param vec DOFVector which is interpolated to the grid.
      * \param filename file to which the grid will be written
      * \param outFilePrecision precision used to write file
      */
    template<typename T>
    void writeGrid(const WorldVector<double> *p,
		  int               *numPoints, 
		  double            *dist,
		  DOFVector<T>      *vec,
		  const char*              filename,
		  int outFilePrecision=6); 

    template<typename T>
    void writeGrid(const WorldVector<double> *p,
		  int               *numPoints, 
		  double            *dist,
		  std::vector< DOFVector<T> * >     vec,
		  const char*              filename,
		  int outFilePrecision=6); 
  }
} } // end namespace io, AMDiS

#include "GridWriter.hh"

#endif
