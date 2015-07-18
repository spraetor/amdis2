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



/** \file ProblemInterpol.h */

#ifndef AMDIS_PROBLEM_INTERPOL_H
#define AMDIS_PROBLEM_INTERPOL_H

#include <vector>
#include "ProblemStat.h"

// TODO: remove this file!

namespace AMDiS {

  /// Interpolates a given function adaptive on spaceProblems mesh.
  class ProblemInterpol : public ProblemStatSeq
  {
  public:
    /** \brief
     * Constructor. fct will be interpolated on the mesh of spaceProblem.
     * grdFct is used, if H1 error should be used for estimation. It points
     * to the gradient of fct.
     */
    ProblemInterpol(const char *name,
		    ProblemStatSeq *spaceProblem,
		    std::vector<AbstractFunction<double, WorldVector<double> >*> *fct,
		    std::vector<AbstractFunction<WorldVector<double>, WorldVector<double> >*> *grdFct);

    /// No system assemblage.
    virtual void buildbeforeRefine(AdaptInfo *adaptInfo, Flag) {}

    /// No system assemblage.
    virtual void buildbeforeCoarsen(AdaptInfo *adaptInfo, Flag) {}

    /// No system assemblage.
    virtual void buildAfterCoarsen(AdaptInfo *adaptInfo, Flag,
				   bool assembleMatrix = true,
				   bool assembleVector = true)
    {}

    /// No equation system ins solved. Instead fct is interpolated to uh.
    virtual void solve(AdaptInfo *adaptInfo,
		       bool createMatrixData = true,
		       bool storeMatrixData = false);

    /// True H1 or L2 error is calculated.
    virtual void estimate(AdaptInfo *adaptInfo);

  protected:
    /// Function to interpolate.
    std::vector<AbstractFunction<double, WorldVector<double> >*> *interpolFct;

    /// Gradient of \ref interpolFct_. Used for H1 error in estimate().
    std::vector<AbstractFunction<WorldVector<double>, WorldVector<double> >*> *grdInterpolFct;
  };

}

#endif
