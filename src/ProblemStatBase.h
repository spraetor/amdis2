/** \file ProblemStatBase.h */

/**
 * \defgroup Problem Problem module
 *  @{ <img src="problem.png"> @}
 */

#pragma once

#include <string>

#include "Flag.h"

namespace AMDiS
{
  // Flags for controling which part of the problem should be initialized

  // For all problems
  const Flag INIT_FE_SPACE         = 0x01L;
  const Flag INIT_MESH             = 0x02L;
  const Flag CREATE_MESH           = 0x04L;
  const Flag INIT_SYSTEM           = 0x08L;
  const Flag INIT_SOLVER           = 0x10L;
  const Flag INIT_ESTIMATOR        = 0x20L;
  const Flag INIT_MARKER           = 0x40L;
  const Flag INIT_ADAPT            = 0x80L;
  const Flag INIT_FILEWRITER       = 0x100L;
  const Flag INIT_GLOBAL_REFINES   = 0x1000L;

  // For time dependent problems
  const Flag INIT_INITIAL_PROBLEM  = 0x200L;
  const Flag INIT_UH_OLD           = 0x400L;

  // For non linear problems
  const Flag INIT_UPDATER          = 0x800L;
  const Flag INIT_NONLIN_SOLVER    = 0x1000L;

  // Combined Flags
  const Flag INIT_NOTHING          = 0x00L;
  const Flag INIT_ALL        = INIT_FE_SPACE | INIT_MESH | CREATE_MESH | INIT_SYSTEM |
                               INIT_SOLVER | INIT_ESTIMATOR | INIT_MARKER |
                               INIT_ADAPT | INIT_FILEWRITER | INIT_INITIAL_PROBLEM |
                               INIT_UH_OLD | INIT_UPDATER | INIT_NONLIN_SOLVER ;

  const Flag MESH_REFINED   = 1;
  const Flag MESH_COARSENED = 2;

  // forward declaration
  class AdaptInfo;

  /**
   * \ingroup Problem
   *
   * \brief
   * Interface for time independent problems. Concrete problems must override
   * all pure virtual methods. The method \ref adaptMethodStat() should
   * initiate the adaption loop which in turn uses the other pure virtual
   * functions. The default stationary adaption loop is implemented in the class
   * AdaptStationary.
   */
  class ProblemStatBase
  {
  public:
    virtual ~ProblemStatBase() {}

    /// Marks mesh elements for refinement and coarsening.
    virtual Flag markElements(AdaptInfo& adaptInfo) = 0;

    /// Assembling of system matrices and vectors before refinement.
    virtual void buildBeforeRefine(AdaptInfo& adaptInfo, Flag flag) = 0;

    /// Assembling of system matrices and vectors before coarsening.
    virtual void buildBeforeCoarsen(AdaptInfo& adaptInfo, Flag flag) = 0;

    /** \brief
     * Assembling of system matrices and vectors after coarsening.
     * By the last two parameters, assembling can be restricted to either
     * matrices or vectors only.
     */
    virtual void buildAfterCoarsen(AdaptInfo& adaptInfo, Flag flag,
                                   bool assembleMatrix, bool assembleVector) = 0;

    /// Refinement of the mesh.
    virtual Flag refineMesh(AdaptInfo& adaptInfo) = 0;

    /// Coarsening of the mesh.
    virtual Flag coarsenMesh(AdaptInfo& adaptInfo) = 0;

    /** \brief
     * Solves the assembled system. The result is an approximative solution.
     * The last two boolean arguments can be used to controll successive
     * solutions of systems with the same matrix.
     *
     * \param  adaptInfo          Pointer to an \ref AdaptInfo object.
     * \param  createMatrixData   If false, the solver assumes that all of its
     *                            internal data structures for the system
     *                            matrix are already created. This is the case,
     *                            if we solve different systems but with the
     *                            same matrix. After the first call to this
     *                            function (with this parameter set to true),
     *                            all other calls may set it to false.
     * \param  storeMatrixData    If true, all internal data structures for the
     *                            system matrix are not deleted such that they
     *                            can be used for next solutions with the same
     *                            system matrix.
     */
    virtual void solve(AdaptInfo& adaptInfo,
                       bool createMatrixData = true,
                       bool storeMatrixData = false) = 0;

    /** \brief
     * A posteriori error estimation of the calculated solution. Should store
     * a local error estimation at each elements leaf data and return the
     * total error sum.
     */
    virtual void estimate(AdaptInfo& adaptInfo) = 0;

    /// Returns the name of the problem.
    virtual std::string getName() const = 0;
  };

} // end namespace AMDiS
