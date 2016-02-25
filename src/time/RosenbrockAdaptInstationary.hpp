/** \file RosenbrockAdaptInstationary.h */

#pragma once

// AMDiS headers
#include "AMDiS_fwd.hpp"
#include "AdaptInstationary.hpp"
#include "CreatorMap.hpp"

#include "time/RosenbrockStationary.hpp"
#include "time/RosenbrockMethod.hpp"

namespace AMDiS
{

  class RosenbrockAdaptInstationary : public AdaptInstationary
  {
  public:
    /** \brief
     * Creates a AdaptInstationary object for Rosenbrock method
     * with the given name for the time
     * dependent problem problemInstat.
     */
    RosenbrockAdaptInstationary(std::string name,
                                RosenbrockStationary& problemStat,
                                AdaptInfo& info,
                                ProblemTimeInterface& problemInstat,
                                AdaptInfo& initialInfo,
                                std::time_t initialTimestamp = 0);

    /// Runs the Rosenbrock loop until one timestep is accepted.
    void oneTimestep();

    virtual double getTimeEst(AdaptInfo& adaptInfo);

    /** \brief
     * This funciton is used only to avoid double code in both constructors. If the
     * obsolte constructure, which uses pointers instead of references, will be
     * removed, remove also this function.
     * TODO: Remove if obsolete constructor will be removed.
     */
    void initConstructor(RosenbrockStationary* problemStat);

    void reset();

  protected:
    /// Pointer to the Rosenbrock method that should be used.
    RosenbrockMethod* rosenbrockMethod;

    /// Pointer to the stationary problem;
    RosenbrockStationary* rosenbrockStat;

    /// Indicates, if this is the very first timestep.
    bool firstTimestep;

    /// If true, the last timestep was rejected.
    bool lastTimestepRejected;


    /// If true, more than one of the last timesteps were rejected.
    bool succRejection;

    /// Maximal nr. of rejected timesteps due to a solver error
    int maxRejectedSolverError;

    /// If greater than 0, than for the first given number of timesteps the timestep
    /// will be not changed and is set to the very first one.
    int fixFirstTimesteps;

    /// Timestep of the last accepted time iteration.
    double tauAcc;

    /// Error estimation of the last accepted time iteration.
    double estAcc;

    /// Timestep of the last rejected time iteration.
    double tauRej;

    /// Error estimation of the last rejected time iteration.
    double estRej;

    /// Current timestep.
    double tau;

    /// The value tau * gamma, where gamma is a value of the used
    /// Rosenbrock method.
    double tauGamma, minusTauGamma, invTauGamma, minusInvTauGamma;

    /// If true, the first timestep is calculated with different timesteps.
    /// This is usually used to make a study how the time error estimator
    /// behavous for different timesteps.
    bool dbgTimestepStudy;

    /// If \ref dbgTimestepStudy is set to true, then this array contains the
    /// timesteps for which the first timestep should be calculated.
    std::vector<double> dbgTimesteps;
  };

} // end namespace AMDiS
