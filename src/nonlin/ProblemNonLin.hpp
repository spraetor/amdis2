#pragma once

// std c++ headers
#include <vector>

// AMDiS includes
#include "DOFVector.hpp"
#include "NonLinSolver.hpp"
#include "ProblemStat.hpp"
#include "SystemVector.hpp"

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
#include "parallel/ParallelProblemStat.hpp"
#endif

namespace AMDiS
{

  /**
   * \ingroup Problem
   *
   * \brief
   * Standard implementation for a vector valued non linear problem.
   */
  class ProblemNonLin : public ProblemStat
  {
  public:
    /// Constructs a ProblemNonLin with given name.
    ProblemNonLin(const std::string& name_)
      : ProblemStat(name_),
        nonLinSolver(NULL)
    {
      u0.resize(nComponents);
      for (int i = 0; i < nComponents; i++)
        u0[i] = NULL;
    }

    /// Sets \ref u0 and interpolates it to \ref solution.
    inline void setU0(AbstractFunction<double, WorldVector<double>>* u0Fct,
                      int index)
    {
      FUNCNAME("ProblemNonLinVec::setU0()");

      TEST_EXIT(index < nComponents)("Invalid index!\n");
      u0[index] = u0Fct;
      solution->getDOFVector(index)->interpol(u0Fct);
    }

    /// Destructor
    virtual ~ProblemNonLin()
    {}

    /// Initialization of the problem.
    virtual void initialize(Flag initFlag,
                            ProblemStatSeq* adoptProblem = NULL,
                            Flag adoptFlag = INIT_NOTHING);

    /// Used in \ref initialize().
    virtual void createNonLinSolver();

    /** \brief
     * Overrides ProblemStat::solve(). Uses the non linear solver
     * \ref nonLinSolver.
     */
    void solve(AdaptInfo& adaptInfo,
               bool createMatrixData = true,
               bool storeMatrixData = false);

    /// Returns \ref nonLinSolver.
    inline NonLinSolver* getNonLinSolver()
    {
      return nonLinSolver;
    }

    /// Sets \ref nonLinSolver.
    inline void setNonLinSolver(NonLinSolver* s)
    {
      nonLinSolver = s;
    }

  protected:
    /// Initial guess for the solution.
    std::vector<AbstractFunction<double, WorldVector<double>>*> u0;

    /// Non linear solver object. Used in \ref solve().
    NonLinSolver* nonLinSolver;
  };

} // end namespace AMDiS

#endif

