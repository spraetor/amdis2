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



/** \file Boundary.h */

#ifndef AMDIS_BOUNDARY_H
#define AMDIS_BOUNDARY_H

#include "Global.h"

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

#endif
