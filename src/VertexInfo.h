/** \file VertexInfo.h */

#pragma once

#include "FixVec.h"
#include "MatrixVectorOperations.h"

namespace AMDiS 
{
  
  /// Stores coordinates and output index for one vertex.
  struct VertexInfo 
  {
    /// Coordinates for this vertex.
    WorldVector<double> coords;
    
    /// Index for the output file.
    int outputIndex;
    
    /// Used to check, whether coords are already stored for a given dof.
    bool operator==(const WorldVector<double>& c) const
    {
      return (c == coords);
    }
    
    /// Used to check, whether coords are already stored for a given dof.
    bool operator!=(const WorldVector<double>& c) const 
    {
      return (c != coords);
    }
  };
  
} // end namespace AMDiS
