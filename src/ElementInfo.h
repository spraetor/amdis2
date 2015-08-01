/** \file ElementInfo.h */

#pragma once

#include <list>
#include <vector>

#include "VertexInfo.h"
#include "Boundary.h"
#include "Projection.h"
#include "AMDiS_fwd.h"

namespace AMDiS 
{
  /// Stores information for one element.
  class ElementInfo
  {
  public:
    /// Dimension specific constructor for DimVec creation.
    ElementInfo(int dim) 
      : vertices(dim),
	vertexInfo(dim),
	boundary(dim),
	projection(dim),
	neighbour(dim),
	surfaceRegions(dim)
    {}
      
    int vertices;
      
    /// Vertex infos for each element vertex.
    DimVec<std::list<VertexInfo>::iterator> vertexInfo;

    /// Boundary type for each side.
    DimVec<BoundaryType> boundary;
  
    /// Boundary projector for each side.
    DimVec<Projection*> projection;
      
    /// Neighbour output index for each side.
    DimVec<int> neighbour;
      
    /// Element type. Used in 3d.
    unsigned char type;

    ///
    int elementRegion;

    ///
    DimVec<int> surfaceRegions;
  };
  
} // end namespace AMDiS
