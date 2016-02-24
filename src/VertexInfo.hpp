#pragma once

#include "FixVec.hpp"
#include "MatrixVectorOperations.hpp"

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
    bool operator==(WorldVector<double> const& c) const
    {
      return (c == coords);
    }

    /// Used to check, whether coords are already stored for a given dof.
    bool operator!=(WorldVector<double> const& c) const
    {
      return (c != coords);
    }
  };

} // end namespace AMDiS
