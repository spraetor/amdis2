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

} // end namespace AMDiS
