/** \file Boundary.h */

#pragma once

namespace AMDiS {

  /// Flag to denote interior boundaryies.
  typedef enum {
    INTERIOR = 0
  } BoundaryConstants;

  /// Type specifier for the different boundary types 
  typedef signed int BoundaryType;

  BoundaryType newBound(BoundaryType oldBound, BoundaryType newBound);
  
  struct BoundaryTypeContainer
  {
    BoundaryTypeContainer(BoundaryType b) : b(b) {}
    BoundaryType b;
  };

} // end namespace AMDiS
