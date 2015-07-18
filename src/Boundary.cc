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


#include "Boundary.h"
#include "FixVec.h"
#include "Initfile.h"

namespace AMDiS {

  BoundaryType newBound(BoundaryType oldBound, BoundaryType newBound)
  {
    if (newBound <= INTERIOR) {
      // Face on NEUMANN-boundary or interior face; weak type.

      return oldBound;
    } else {
      // Node is already node on the DIRICHLET boundary.

      if (oldBound > newBound)
	return oldBound;
      else
	return newBound;
    }

    // New face is interior face; node type is always stronger.

    return newBound;
  }

}
