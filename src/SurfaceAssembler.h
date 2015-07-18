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



/** \file SurfaceAssembler.h */

#ifndef AMDIS_SURFACEASSEMBLER_H
#define AMDIS_SURFACEASSEMBLER_H

#include "AMDiS_fwd.h"
#include "Assembler.h"
#include "FixVec.h"

namespace AMDiS {

  /** 
   * \ingroup Integration
   *
   * \brief
   */
  class SurfaceAssembler : public Assembler
  {
  public:
    /// Creates a SurfaceAssembler conforming to operate for the given \ref coords.
    SurfaceAssembler(Operator *operat,
		     const FiniteElemSpace *rowFeSpace,
		     const FiniteElemSpace *colFeSpace,
		     VectorOfFixVecs<DimVec<double> > &coords);

    /// Destructor
    ~SurfaceAssembler();

    /// Adapt surface quadratures to \ref coords.
    void adaptSurfaceAssembler(VectorOfFixVecs<DimVec<double> > &coords);

    ///
    bool initElementVector(const ElInfo *elInfo);

  protected:
    VectorOfFixVecs<DimVec<double> > coords_;
  };

}

#endif
