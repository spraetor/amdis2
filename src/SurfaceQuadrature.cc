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


#include <algorithm>

#include "SurfaceQuadrature.h"
#include "Quadrature.h"

namespace AMDiS {

  SurfaceQuadrature::SurfaceQuadrature(Quadrature *q,
				       VectorOfFixVecs<DimVec<double> > &c)
    : Quadrature((q->getName() + " surface").c_str(),
		 q->getDegree(),
		 q->getDim() + 1,
		 q->getNumPoints(),
		 NULL,
		 q->getWeight()), 
      quad(q),
      coords(c)
  {
    lambda = new VectorOfFixVecs<DimVec<double> >(dim, n_points, NO_INIT);

    // for each integration point
    for (int i = 0; i < n_points; i++) {
      // get coords of quadrature point in dim-1
      DimVec<double> const& origin = quad->getLambda(i);

      for (int j = 0; j < dim + 1; j++)
	(*lambda)[i][j] = 0.0;

      for (int j = 0; j < dim; j++)
	for (int k = 0; k < dim + 1; k++)
	  (*lambda)[i][k] += origin[j] * coords[j][k];
    }
  }

  void SurfaceQuadrature::scaleSurfaceQuadrature(VectorOfFixVecs<DimVec<double> >&c)
  {
    // copy coords NOTE: why???
    for (int i = 0; i < dim; i++)
      coords[i] = c[i];
  
    // for each integration point
    for (int i = 0; i < n_points; i++) {
      // get coords of quadrature point in dim-1
      DimVec<double> const& origin = quad->getLambda(i);

      for (int j = 0; j < dim + 1; j++)
	(*lambda)[i][j] = 0.0;

      for (int j = 0; j < dim; j++)
	for (int k = 0; k < dim+1; k++)
	  (*lambda)[i][k] += origin[j] * coords[j][k];
    }
  }

}
