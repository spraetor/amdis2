/** \file AdaptInstationary.h */

#pragma once

#include <string>
#include <ctime>
#include <queue>

#include "AMDiS_fwd.h"
#include "AdaptBase.h"

namespace AMDiS {

  /** \ingroup Adaption  
   * \brief
   * AdaptInstationary implements the adaptive procdure for time dependent 
   * problems (see ProblemInstat). It contains a pointer to a ProblemInstat
   * object.
   */
  class AdaptInstationary : public AdaptBase
  {
  public:
    /// Creates a AdaptInstationary object with the given name for the time 
    /// dependent problem problemInstat.
    AdaptInstationary(std::string name, 
            		      ProblemIterationInterface &problemStat,
            		      AdaptInfo &info,
            		      ProblemTimeInterface &problemInstat,
            		      AdaptInfo &initialInfo,
                      std::time_t initialTimestampSet = 0);

    /// Sets \ref strategy to aStrategy
    void setStrategy(int aStrategy) 
    { 
      strategy = aStrategy; 
    }

    /// Returns \ref strategy
    int getStrategy() const 
    {
      return strategy;
    }

    /// Implementation of AdaptBase::adapt()
    virtual int adapt() override;

  protected:
    /** \brief
     * Implements one (maybe adaptive) timestep. Both the explicit and the 
     * implicit time strategy are implemented. The semi-implicit strategy 
     * is only a special case of the implicit strategy with a limited number of 
     * iterations (exactly one).
     * The routine uses the parameter \ref strategy to select the strategy:
     * strategy 0: Explicit strategy, 
     * strategy 1: Implicit strategy.
     */
    virtual void oneTimestep();

    /// Initialisation of this AdaptInstationary object
    void initialize(std::string aName);

    /// Implements the explit time strategy. Used by \ref oneTimestep().
    virtual void explicitTimeStrategy();

    /// Implements the implicit time strategy. Used by \ref oneTimestep().
    virtual void implicitTimeStrategy();

    /** \brief
     * This iteration strategy allows the timestep and the mesh to be adapted 
     * after each timestep solution. There are no inner loops for mesh adaption and
     * no refused timesteps.
     */
    void simpleAdaptiveTimeStrategy();

  protected:
    /// Strategy for choosing one timestep
    int strategy; 

    /// Parameter \f$ \delta_1 \f$ used in time step reduction
    double timeDelta1;

    /// Parameter \f$ \delta_2 \f$ used in time step enlargement
    double timeDelta2;

    /// If this parameter is 1 and the instationary problem is stable, hence the number
    /// of solver iterations to solve the problem is zero, the adaption loop will stop.
    int breakWhenStable;

    ///
    bool fixedTimestep;
  };

} // end namespace AMDiS
