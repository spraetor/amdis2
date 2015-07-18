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



/** \file ElementInfo.h */

#ifndef AMDIS_ELEMENTINFO_H
#define AMDIS_ELEMENTINFO_H

#include <list>
#include <vector>
#include "VertexInfo.h"
#include "Boundary.h"
#include "Projection.h"
#include "AMDiS_fwd.h"

namespace AMDiS {

  /// Stores information for one element.
  class ElementInfo
  {
  public:
    /// Dimension specific constructor for DimVec creation.
    ElementInfo(int dim) 
      : vertices(dim),
	vertexInfo(dim, NO_INIT),
	boundary(dim, NO_INIT),
	projection(dim, NO_INIT),
	neighbour(dim, NO_INIT),
	surfaceRegions(dim, NO_INIT)
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
  
}

#endif
  
