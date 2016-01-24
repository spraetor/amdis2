/** \file ProblemInterpol.h */

#pragma once

#include <vector>

#include "ProblemStat.h"

// TODO: move to deprecated

namespace AMDiS
{

  /// Interpolates a given function adaptive on spaceProblems mesh.
  class ProblemInterpol : public ProblemStatSeq
  {
  public:
    /** \brief
     * Constructor. fct will be interpolated on the mesh of spaceProblem.
     * grdFct is used, if H1 error should be used for estimation. It points
     * to the gradient of fct.
     */
    ProblemInterpol(const char* name,
                    ProblemStatSeq* spaceProblem,
                    std::vector<AbstractFunction<double, WorldVector<double>>*>* fct,
                    std::vector<AbstractFunction<WorldVector<double>, WorldVector<double>>*>* grdFct);

    /// No system assemblage.
    virtual void buildbeforeRefine(AdaptInfo& adaptInfo, Flag) {}

    /// No system assemblage.
    virtual void buildbeforeCoarsen(AdaptInfo& adaptInfo, Flag) {}

    /// No system assemblage.
    virtual void buildAfterCoarsen(AdaptInfo& adaptInfo, Flag,
                                   bool assembleMatrix = true,
                                   bool assembleVector = true)
    {}

    /// No equation system ins solved. Instead fct is interpolated to uh.
    virtual void solve(AdaptInfo& adaptInfo,
                       bool createMatrixData = true,
                       bool storeMatrixData = false);

    /// True H1 or L2 error is calculated.
    virtual void estimate(AdaptInfo& adaptInfo);

  protected:
    /// Function to interpolate.
    std::vector<AbstractFunction<double, WorldVector<double>>*>* interpolFct;

    /// Gradient of \ref interpolFct_. Used for H1 error in estimate().
    std::vector<AbstractFunction<WorldVector<double>, WorldVector<double>>*>* grdInterpolFct;
  };

} // end namespace AMDiS
