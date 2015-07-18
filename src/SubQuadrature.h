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



/** \file SubQuadrature.h */

#ifndef AMDIS_SUBQUADRATURE_H
#define AMDIS_SUBQUADRATURE_H

#include "AMDiS_fwd.h"
#include "Quadrature.h"

namespace AMDiS {

  class SubQuadrature : public Quadrature
  {
  public:
    SubQuadrature(Quadrature *quad, int dim_);

    void scaleQuadrature(VectorOfFixVecs<DimVec<double> > &coords);

    int getSubDim() const
    { 
      return subDim_;
    }

  protected:
    Quadrature *quad_;

    int subDim_;
  };

}

#endif
