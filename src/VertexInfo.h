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



/** \file VertexInfo.h */

#ifndef AMDIS_VERTEXINFO_H
#define AMDIS_VERTEXINFO_H

#include "FixVec.h"
#include "MatrixVectorOperations.h"

namespace AMDiS {
  
  /// Stores coordinates and output index for one vertex.
  struct VertexInfo 
  {
    /// Coordinates for this vertex.
    WorldVector<double> coords;
    
    /// Index for the output file.
    int outputIndex;
    
    /// Used to check, whether coords are already stored for a given dof.
    bool operator==(const WorldVector<double>& c) 
    {
      return (c == coords);
    }
    
    /// Used to check, whether coords are already stored for a given dof.
    bool operator!=(const WorldVector<double>& c) 
    {
      return (c != coords);
    }
  };
  
} // end namespace AMDiS

#endif
