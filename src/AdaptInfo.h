/** \file AdaptInfo.h */

#pragma once

#include <vector>
#include <string>
#include <cmath>

#include "Log.h"
#include "Math.hpp"

namespace AMDiS
{

  /**
   * \ingroup Adaption
   *
   * \brief
   * Holds adapt parameters and infos about the problem. Base class
   * for AdaptInfoScal and AdaptInfoVec.
   */
  class AdaptInfo
  {
  protected:
    /** \brief
     * Stores adapt infos for a scalar problem or for one component of a
     * vector valued problem.
     */
    class ScalContent
    {
    public:
      /// Constructor.
      explicit ScalContent(std::string prefix);

      /// Sum of all error estimates
      double est_sum = 0.0;

      /// Sum of all time error estimates
      double est_t_sum = 0.0;

      /// maximal local error estimate.
      double est_max = 0.0;

      /// Maximum of all time error estimates
      double est_t_max = 0.0;

      /// factors to combine max and integral time estimate
      double fac_max = 0.0, fac_sum = 1.0;

      /// Tolerance for the (absolute or relative) error
      double spaceTolerance = 0.0;

      /// Time tolerance.
      double timeTolerance = 0.0;

      /// Relative time tolerance
      double timeRelativeTolerance = 0.0;

      /// Lower bound for the time error.
      double timeErrLow = 0.0;

      /// true if coarsening is allowed, false otherwise.
      int coarsenAllowed = 0;

      /// true if refinement is allowed, false otherwise.
      int refinementAllowed = 1;

      /** \brief
       * parameter to tell the marking strategy how many bisections should be
       * performed when an element is marked for refinement; usually the value is
       * 1 or DIM
       */
      int refineBisections = 1;

      /** \brief
       * parameter to tell the marking strategy how many bisections should
       * be undone when an element is marked for coarsening; usually the value is
       * 1 or DIM
       */
      int coarseBisections = 1;
    };

  public:
    /// Constructor.
    AdaptInfo(std::string name_, int size = 1);

    /// Destructor.
    virtual ~AdaptInfo()
    {
      for (size_t i = 0;  i < scalContents.size(); i++)
        delete scalContents[i];
    }

    /// Resets all variables to zero (or something equivalent)
    void reset();

    /// Returns whether space tolerance is reached.
    virtual bool spaceToleranceReached() const
    {
      for (size_t i = 0; i < scalContents.size(); i++)
      {
        if (!(scalContents[i]->est_sum < scalContents[i]->spaceTolerance))
          return false;
      }

      return true;
    }

    /// Returns whether space tolerance of component i is reached.
    virtual bool spaceToleranceReached(int i) const
    {
      if (!(scalContents[i]->est_sum < scalContents[i]->spaceTolerance))
        return false;
      else
        return true;
    }

    /// Returns whether time tolerance is reached.
    virtual bool timeToleranceReached() const
    {
      for (size_t i = 0; i < scalContents.size(); i++)
        if (!(getTimeEstCombined(i) < scalContents[i]->timeTolerance))
          return false;

      return true;
    }

    /// Returns whether time tolerance of component i is reached.
    virtual bool timeToleranceReached(int i) const
    {
      if (!(getTimeEstCombined(i) < scalContents[i]->timeTolerance))
        return false;
      else
        return true;
    }

    /// Returns whether time error is under its lower bound.
    virtual bool timeErrorLow() const
    {
      for (size_t i = 0; i < scalContents.size(); i++)
        if (!(getTimeEstCombined(i) < scalContents[i]->timeErrLow))
          return false;

      return true;
    }
    /// Returns the time estimation as a combination
    /// of maximal and integral time error
    double getTimeEstCombined(unsigned i) const
    {
      return
        scalContents[i]->est_t_max * scalContents[i]->fac_max +
        scalContents[i]->est_t_sum * scalContents[i]->fac_sum;
    }


    /// Print debug information about time error and its bound.
    void printTimeErrorLowInfo() const;

    /// Returns \ref spaceIteration.
    int getSpaceIteration() const
    {
      return spaceIteration;
    }

    /// Sets \ref spaceIteration.
    void setSpaceIteration(int it)
    {
      spaceIteration = it;
    }

    /// Returns \ref maxSpaceIteration.
    int getMaxSpaceIteration() const
    {
      return maxSpaceIteration;
    }

    /// Sets \ref maxSpaceIteration.
    void setMaxSpaceIteration(int it)
    {
      maxSpaceIteration = it;
    }

    /// Increments \ref spaceIteration by 1;
    void incSpaceIteration()
    {
      spaceIteration++;
    }

    /// Sets \ref timestepIteration.
    void setTimestepIteration(int it)
    {
      timestepIteration = it;
    }

    /// Returns \ref timestepIteration.
    int getTimestepIteration() const
    {
      return timestepIteration;
    }

    /// Increments \ref timestepIteration by 1;
    void incTimestepIteration()
    {
      timestepIteration++;
    }

    /// Returns \ref maxTimestepIteration.
    int getMaxTimestepIteration() const
    {
      return maxTimestepIteration;
    }

    /// Sets \ref maxTimestepIteration.
    void setMaxTimestepIteration(int it)
    {
      maxTimestepIteration = it;
    }

    /// Sets \ref timeIteration.
    void setTimeIteration(int it)
    {
      timeIteration = it;
    }

    /// Returns \ref timeIteration.
    int getTimeIteration() const
    {
      return timeIteration;
    }

    /// Increments \ref timesIteration by 1;
    void incTimeIteration()
    {
      timeIteration++;
    }

    /// Returns \ref maxTimeIteration.
    int getMaxTimeIteration() const
    {
      return maxTimeIteration;
    }

    /// Sets \ref maxTimeIteration.
    void setMaxTimeIteration(int it)
    {
      maxTimeIteration = it;
    }

    /// Returns \ref timestepNumber.
    int getTimestepNumber() const
    {
      return timestepNumber;
    }

    /// Sets \ref timestepNumber.
    void setTimestepNumber(int num)
    {
      timestepNumber = std::min(nTimesteps, num);
    }

    /// Returns \ref nTimesteps.
    int getNumberOfTimesteps() const
    {
      return nTimesteps;
    }

    /// Sets \ref nTimesteps.
    void setNumberOfTimesteps(int num)
    {
      nTimesteps = std::max(0, num);
    }

    /// Increments \ref timestepNumber by 1;
    void incTimestepNumber()
    {
      timestepNumber++;
    }

    /// Sets \ref est_sum.
    void setEstSum(double e, int index)
    {
      scalContents[index]->est_sum = e;
    }

    /// Sets \ref est_max.
    void setEstMax(double e, int index)
    {
      scalContents[index]->est_max = e;
    }

    /// Sets \ref est_max.
    void setTimeEstMax(double e, int index)
    {
      scalContents[index]->est_t_max = e;
    }

    /// Sets \ref est_t_sum.
    void setTimeEstSum(double e, int index)
    {
      scalContents[index]->est_t_sum = e;
    }

    /// Returns \ref est_sum.
    double getEstSum(int index) const
    {
      FUNCNAME_DBG("AdaptInfo::getEstSum()");
      TEST_EXIT_DBG(static_cast<size_t>(index) < scalContents.size())
      ("Wrong index for adaptInfo!\n");

      return scalContents[index]->est_sum;
    }

    /// Returns \ref est_t_sum.
    double getEstTSum(int index) const
    {
      return scalContents[index]->est_t_sum;
    }

    /// Returns \ref est_max.
    double getEstMax(int index) const
    {
      FUNCNAME_DBG("AdaptInfo::getEstSum()");
      TEST_EXIT_DBG(static_cast<size_t>(index) < scalContents.size())
      ("Wrong index for adaptInfo!\n");

      return scalContents[index]->est_max;
    }

    /// Returns \ref est_max.
    double getTimeEstMax(int index) const
    {
      return scalContents[index]->est_t_max;
    }

    /// Returns \ref est_t_sum.
    double getTimeEstSum(int index) const
    {
      return scalContents[index]->est_t_sum;
    }

    /// Returns \ref est_t the estimated overall time error
    double getTimeEst() const
    {
      return est_t;
    }

    void setTimeEst(double value)
    {
      est_t = value;
    }

    /// Returns \ref spaceTolerance.
    double getSpaceTolerance(int index) const
    {
      return scalContents[index]->spaceTolerance;
    }

    /// Sets \ref spaceTolerance.
    void setSpaceTolerance(int index, double tol)
    {
      scalContents[index]->spaceTolerance = tol;
    }

    /// Returns \ref timeTolerance.
    double getTimeTolerance(int index) const
    {
      return scalContents[index]->timeTolerance;
    }

    /// Returns \ref timeRelativeTolerance.
    double getTimeRelativeTolerance(int index) const
    {
      return scalContents[index]->timeRelativeTolerance;
    }

    /// Sets \ref time
    double setTime(double t)
    {
      time = t;
      if (time > endTime)
        time = endTime;
      if (time < startTime)
        time = startTime;

      return time;
    }

    /// Gets \ref time
    double getTime() const
    {
      return time;
    }

    /// Gets \ref &time
    double* getTimePtr()
    {
      return &time;
    }

    /// Sets \ref timestep
    double setTimestep(double t)
    {
      timestep = t;
      if (timestep > maxTimestep)
        timestep = maxTimestep;
      if (timestep < minTimestep)
        timestep = minTimestep;
      if (time + timestep > endTime)
        timestep = endTime - time;

      return timestep;
    }
    /// Gets \ref timestep
    double getTimestep() const
    {
      return timestep;
    }

    void setLastProcessedTimestep(double t)
    {
      lastProcessedTimestep = t;
    }

    double getLastProcessedTimestep() const
    {
      return lastProcessedTimestep;
    }

    /// Returns true, if the end time is reached and no more timestep
    /// computations must be done.
    bool reachedEndTime() const
    {
      if (nTimesteps > 0)
        return !(timestepNumber < nTimesteps);

      return !(std::abs(time - endTime) > DBL_TOL);
    }


    /// Sets \ref minTimestep
    void setMinTimestep(double t)
    {
      minTimestep = t;
    }

    /// Gets \ref minTimestep
    double getMinTimestep() const
    {
      return minTimestep;
    }

    /// Sets \ref maxTimestep
    void setMaxTimestep(double t)
    {
      maxTimestep = t;
    }

    /// Gets \ref maxTimestep
    double getMaxTimestep() const
    {
      return maxTimestep;
    }

    /// Gets \ref &timestep
    double* getTimestepPtr()
    {
      return &timestep;
    }

    /// Sets \ref startTime = time
    void setStartTime(double time)
    {
      startTime = time;
    }

    /// Sets \ref endTime = time
    void setEndTime(double time)
    {
      endTime = time;
    }

    /// Returns \ref startTime
    double getStartTime() const
    {
      return startTime;
    }

    /// Returns \ref endTime
    double getEndTime() const
    {
      return endTime;
    }

    /// Returns \ref timeErrLow.
    double getTimeErrLow(int index) const
    {
      return scalContents[index]->timeErrLow;
    }

    /// Returns whether coarsening is allowed or not.
    bool isCoarseningAllowed(int index) const
    {
      return (scalContents[index]->coarsenAllowed == 1);
    }

    /// Returns whether coarsening is allowed or not.
    bool isRefinementAllowed(int index) const
    {
      return (scalContents[index]->refinementAllowed == 1);
    }

    ///
    void allowRefinement(bool allow, int index)
    {
      scalContents[index]->refinementAllowed = allow;
    }

    ///
    void allowCoarsening(bool allow, int index)
    {
      scalContents[index]->coarsenAllowed = allow;
    }

    /// Returns \ref refineBisections
    int getRefineBisections(int index) const
    {
      return scalContents[index]->refineBisections;
    }

    /// Returns \ref coarseBisections
    int getCoarseBisections(int index) const
    {
      return scalContents[index]->coarseBisections;
    }

    int getSize() const
    {
      return static_cast<int>(scalContents.size());
    }

    // TODO: remove from AdaptInfo
    bool getRosenbrockMode() const
    {
      return rosenbrockMode;
    }

    void setSolverIterations(int it)
    {
      solverIterations = it;
    }

    int getSolverIterations() const
    {
      return solverIterations;
    }

    void setMaxSolverIterations(int it)
    {
      maxSolverIterations = it;
    }

    int getMaxSolverIterations() const
    {
      return maxSolverIterations;
    }

    void setSolverTolerance(double tol)
    {
      solverTolerance = tol;
    }

    double getSolverTolerance() const
    {
      return solverTolerance;
    }

    void setSolverResidual(double res)
    {
      solverResidual = res;
    }

    double getSolverResidual() const
    {
      return solverResidual;
    }

    void setGlobalTimeTolerance(double tol)
    {
      globalTimeTolerance = tol;
    }

    double getGlobalTimeTolerance() const
    {
      return globalTimeTolerance;
    }

    // TODO: remove from AdaptInfo
    void setRosenbrockMode(bool b)
    {
      rosenbrockMode = b;
    }

    /// Creates new scalContents with the given size.
    void setScalContents(int newSize);

    /** \brief
     * Resets timestep, current time and time boundaries without
     * any check. Is used by the parareal algorithm.
     */
    void resetTimeValues(double newTimeStep,
                         double newStartTime,
                         double newEndTime)
    {
      time = newStartTime;
      startTime = newStartTime;
      endTime = newEndTime;
      timestep = newTimeStep;
      timestepNumber = 0;
    }

  protected:
    /// Name.
    std::string name;

    /// Current space iteration
    int spaceIteration = -1;

    /** \brief
     * maximal allowed number of iterations of the adaptive procedure; if
     * maxIteration <= 0, no iteration bound is used
     */
    int maxSpaceIteration = -1;

    /// Current timestep iteration
    int timestepIteration = 0;

    /// Maximal number of iterations for choosing a timestep
    int maxTimestepIteration = 30;

    /// Current time iteration
    int timeIteration = 0;

    /// Maximal number of time iterations
    int maxTimeIteration = 30;

    /// Actual time, end of time interval for current time step
    double time = 0.0;

    /// Initial time
    double startTime = 0.0;

    /// Final time
    double endTime = 1.0;

    ///Time step size to be used
    double timestep = 0.0;

    /// Last processed time step size of finished iteration
    double lastProcessedTimestep = 0.0;

    /// Minimal step size
    double minTimestep = 0.0;

    /// Maximal step size
    double maxTimestep = 1.0;

    /// Number of current time step
    int timestepNumber = 0;

    /** \brief
     * Per default this value is 0 and not used. If it is set to a non-zero value,
     * the computation of the stationary problem is done nTimesteps times with a
     * fixed timestep.
     */
    int nTimesteps = 0;

    /// number of iterations needed of linear or nonlinear solver
    int solverIterations = 0;

    /// maximal number of iterations needed of linear or nonlinear solver
    int maxSolverIterations = 0;

    ///
    double solverTolerance = 1.e-8;

    ///
    double solverResidual = 0.0;

    /// tolerance for the overall time error
    double globalTimeTolerance = 1.0;

    /// Scalar adapt infos.
    std::vector<ScalContent*> scalContents;

    /// Is true, if the adaptive procedure was deserialized from a file. TODO: remove deserialization
    bool deserialized = false;

    /// Is true, if the time adaption is controlled by a Rosenbrock Method.
    bool rosenbrockMode = false;

    /// overall time error estimate
    double est_t = 0.0;
  };

} // end namespace AMDiS
