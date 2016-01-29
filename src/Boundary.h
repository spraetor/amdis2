/** \file Boundary.h */

#pragma once

namespace AMDiS
{

  /// Flag to denote interior boundaryies.
  typedef enum
  {
    INTERIOR = 0
  } BoundaryConstants;

  /// Type specifier for the different boundary types
  using BoundaryType = signed int;

  BoundaryType newBound(BoundaryType oldBound, BoundaryType newBound);

  struct BoundaryTypeContainer
  {
    explicit BoundaryTypeContainer(BoundaryType b) : b(b) {}
    BoundaryType b;
  };

} // end namespace AMDiS
