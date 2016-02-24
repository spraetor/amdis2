#pragma once

#include "AdaptInstationary.hpp"
#include "ProblemStat.hpp"
#include "ProblemTimeInterface.hpp"

namespace AMDiS
{
  /**
   * \ingroup Problem
   *
   * \brief
   * Base class for \ref ProblemInstat.
   */
  class ProblemInstatBase : public ProblemTimeInterface,
			    public ProblemStatBase   // NOTE: Why is this derived from ProblemStatBase
  {
  public:
    /// Constructor.
    ProblemInstatBase(std::string probName,
                      ProblemStatBase* initialProb)
      : name(probName),
        initialProblem(initialProb ? initialProb : this)
    {}

    /// Destructor.
    virtual ~ProblemInstatBase() {}

    /// Implementation of \ref ProblemTimeInterface::setTime().
    virtual void setTime(AdaptInfo& adaptInfo) override
    {
      cTime = adaptInfo.getTime();
      tau = adaptInfo.getTimestep();
      invTau = 1.0 / tau;
    }

    void solve(AdaptInfo& adaptInfo) {}

    /// Implementation of \ref ProblemStatBase::solve().
    virtual void solve(AdaptInfo& adaptInfo, bool, bool) override
    {
      solve(adaptInfo);
    }

    /// Implementation of \ref ProblemStatBase::estimate().
    virtual void estimate(AdaptInfo& adaptInfo) override {}

    /// Implementation of \ref ProblemStatBase::buildBeforeRefine().
    virtual void buildBeforeRefine(AdaptInfo& adaptInfo, Flag) override {}

    /// Implementation of \ref ProblemStatBase::buildBeforeCoarsen().
    virtual void buildBeforeCoarsen(AdaptInfo& adaptInfo, Flag) override {}

    /// Implementation of \ref ProblemStatBase::buildAfterCoarsen().
    virtual void buildAfterCoarsen(AdaptInfo& adaptInfo, Flag, bool, bool) override {}

    /// Implementation of \ref ProblemStatBase::markElements().
    virtual Flag markElements(AdaptInfo& adaptInfo) override
    {
      return {0};
    }

    /// Implementation of \ref ProblemStatBase::refineMesh().
    virtual Flag refineMesh(AdaptInfo& adaptInfo) override
    {
      return {0};
    }

    /// Implementation of \ref ProblemStatBase::coarsenMesh().
    virtual Flag coarsenMesh(AdaptInfo& adaptInfo) override
    {
      return {0};
    }

    /// Implementation of \ref ProblemTimeInterface::closeTimestep().
    virtual void closeTimestep(AdaptInfo& adaptInfo) override {}

    /// Implementation of \ref ProblemStatBase::getName().
    std::string getName() const
    {
      return name;
    }

    /// Implementation of \ref ProblemTimeInterface::solveInitialProblem().
    virtual void solveInitialProblem(AdaptInfo& adaptInfo) override;

    double* getTime()
    {
      return &cTime;
    }

    double* getTau()
    {
      return &tau;
    }

    double* getInvTau()
    {
      return &invTau;
    }

  protected:
    /// Name of the problem.
    std::string name;

    ProblemStatBase* initialProblem;

    /// Time
    double cTime = 0.0;

    /// Timestep
    double tau = 1.0;

    /// 1 / timestep
    double invTau = 1.0;
  };


  /**
   * \ingroup Problem
   *
   * \brief
   * Standard implementation of ProblemTimeInterface for a time
   * dependent problems.
   */
  class ProblemInstat : public ProblemInstatBase
  {
  public:
    /// Constructs a ProblemInstatVec with prob as its stationary problem.
    ProblemInstat(std::string name, 
		  ProblemStatSeq& prob)
      : ProblemInstatBase(name, NULL),
	problemStat(prob)
    {}
    
    ProblemInstat(std::string name, 
		  ProblemStatSeq& prob,
                  ProblemStatBase& initialProb)
      : ProblemInstatBase(name, &initialProb),
	problemStat(prob)
    {}

    /// Destructor.
    virtual ~ProblemInstat();

    /// Initialisation of the problem.
    virtual void initialize(Flag initFlag,
                            ProblemInstat* adoptProblem = NULL,
                            Flag adoptFlag = INIT_NOTHING);

    /// Used in \ref initialize().
    virtual void createUhOld();

    /// Implementation of \ref ProblemTimeInterface::initTimestep().
    virtual void initTimestep(AdaptInfo& adaptInfo) override;

    /// Implementation of \ref ProblemTimeInterface::closeTimestep().
    virtual void closeTimestep(AdaptInfo& adaptInfo) override;

    /// Returns \ref problemStat.
    ProblemStatSeq& getStatProblem()
    {
      return problemStat;
    }

    /// Returns \ref oldSolution.
    SystemVector* getOldSolution()
    {
      return oldSolution;
    }

    DOFVector<double>* getOldSolution(int i)
    {
      return oldSolution->getDOFVector(i);
    }

    /// Implementation of \ref ProblemTimeInterface::transferInitialSolution().
    virtual void transferInitialSolution(AdaptInfo& adaptInfo);

  protected:
    /// Space problem solved in each timestep.
    ProblemStatSeq& problemStat;

    /// Solution of the last timestep.
    SystemVector* oldSolution = NULL;

    /// In parallel computations, we want to print the overall computational time
    /// that is used for one timestep.
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    double lastTimepoint = 0.0;
#endif
  };

} // end namespace AMDiS
