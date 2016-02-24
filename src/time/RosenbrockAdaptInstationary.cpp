#include "time/RosenbrockAdaptInstationary.hpp"

// AMDiS headers
#include "AdaptInfo.hpp"
#include "ProblemTimeInterface.hpp"

namespace AMDiS
{
  RosenbrockAdaptInstationary::RosenbrockAdaptInstationary(std::string name,
							   RosenbrockStationary& problemStat,
							   AdaptInfo& info,
							   ProblemTimeInterface& problemInstat,
							   AdaptInfo& initialInfo,
							   std::time_t initialTimestamp)
    : AdaptInstationary(name, problemStat, info, problemInstat, initialInfo, initialTimestamp),
      rosenbrockStat(&problemStat),
      dbgTimestepStudy(false)
  {
    FUNCNAME_DBG("RosenbrockAdaptInstationary::initConstructor()");

    std::string initFileStr(name + "->rosenbrock method");
    std::string str = Parameters::get<std::string>(initFileStr, "");
    RosenbrockMethodCreator* creator =
      dynamic_cast<RosenbrockMethodCreator*>(CreatorMap<RosenbrockMethod>::getCreator(str, initFileStr));
    rosenbrockMethod = creator->create();

    TEST_EXIT_DBG(rosenbrockMethod)("Rosenbrock method could not be created!\n");
    
    // set initial parameters
    reset();

    Parameters::get(name + "->fix first timesteps", fixFirstTimesteps);
    problemStat.setRosenbrockMethod(rosenbrockMethod);

    adaptInfo.setRosenbrockMode(true);

    problemStat.setTauGamma(&tauGamma, &minusTauGamma,
                            &invTauGamma, &minusInvTauGamma);
    problemStat.setTau(&tau);

    Parameters::get(name + "->rosenbrock->timestep study", dbgTimestepStudy);
    Parameters::get(name + "->rosenbrock->timestep study steps", dbgTimesteps);

    Parameters::get(name + "->rosenbrock->max rejected timesteps", maxRejectedSolverError);
    
    rosenbrockStat->setOldTime(adaptInfo.getStartTime());
  }


  void RosenbrockAdaptInstationary::reset()
  {
    firstTimestep = true;
    lastTimestepRejected = false;
    succRejection = false;
    fixFirstTimesteps = 0;
    tau = 1.0;
    tauGamma = 1.0;
    minusTauGamma = -1.0;
    invTauGamma = 1.0;
    minusInvTauGamma = -1.0;
  }


  void RosenbrockAdaptInstationary::oneTimestep()
  {
    FUNCNAME("RosenbrockAdaptInstationary::oneTimestep()");

    // estimate before first adaption
    if (adaptInfo.getTime() <= adaptInfo.getStartTime())
      problemIteration->oneIteration(adaptInfo, ESTIMATE);

    bool rejected = false;
    double timeTol = adaptInfo.getGlobalTimeTolerance();

    int studyTimestep = -1;
    int rejectedByError = 0;

    do
    {
      if (dbgTimestepStudy)
      {
        if (studyTimestep >= 0)
          adaptInfo.setTime(adaptInfo.getTime() - dbgTimesteps[studyTimestep]);

        studyTimestep++;
        TEST_EXIT(studyTimestep < static_cast<int>(dbgTimesteps.size()))("Should not happen!\n");

        adaptInfo.setTimestep(dbgTimesteps[studyTimestep]);
      }

      // increment time
      adaptInfo.setTime(adaptInfo.getTime() + adaptInfo.getTimestep());
      problemTime->setTime(adaptInfo);
      tau = adaptInfo.getTimestep();

      tauGamma = tau * rosenbrockMethod->getGamma();
      minusTauGamma = -tauGamma;
      invTauGamma = 1.0 / (tau * rosenbrockMethod->getGamma());
      minusInvTauGamma = -1.0;

      INFO(info, 6)("time = %e   timestep = %e\n",
                    adaptInfo.getTime(), adaptInfo.getTimestep());

      adaptInfo.setSpaceIteration(0);

      // do the iteration
      problemIteration->beginIteration(adaptInfo);
      bool solverError = false;
      try
      {
        problemIteration->oneIteration(adaptInfo, FULL_ITERATION);
      }
      catch(...)
      {
        rejectedByError++;
        solverError = true;
      }
      problemIteration->endIteration(adaptInfo);

      TEST_EXIT(rejectedByError < maxRejectedSolverError)
      ("Number of failed timestep iterations, because of solver errors, reached a limit!\n");
      if (!solverError)
        rejectedByError = 0;

      double errorEst = getTimeEst(adaptInfo);

      double newTimestep = 0.0;
      double order = rosenbrockMethod->getOrder();

      if (errorEst < timeTol || fixFirstTimesteps > 0 || firstTimestep || solverError)
      {
        double fac = 1.0;

        if (solverError)
        {
          newTimestep = adaptInfo.getTimestep() * 0.5;
        }
        else if (fixFirstTimesteps > 0)
        {
          newTimestep = adaptInfo.getTimestep();
        }
        else
        {
          if (firstTimestep || succRejection)
          {
            fac = std::pow((timeTol / errorEst), 1.0 / order);
          }
          else
          {
            fac = adaptInfo.getTimestep() / tauAcc *
                  std::pow((timeTol * estAcc / (errorEst * errorEst)), 1.0 / order);
          }
          fac = std::max(0.1, std::min(fac, 3.0));
          newTimestep = fac * adaptInfo.getTimestep();
          newTimestep *= 0.95;
        }

        tauAcc = adaptInfo.getTimestep();
        estAcc = errorEst;

        lastTimestepRejected = false;
        succRejection = false;
      }
      else
      {
        if (lastTimestepRejected)
        {
          succRejection = true;

          double reducedP =
            std::log(errorEst / estRej) / std::log(adaptInfo.getTimestep() / tauRej);

          if (reducedP < order && reducedP > 0.0)
            order = reducedP;
        }

        newTimestep = std::pow((timeTol / errorEst), 1.0 / order) * adaptInfo.getTimestep();
        newTimestep *= 0.95;

        tauRej = adaptInfo.getTimestep();
        estRej = errorEst;

        lastTimestepRejected = true;
      }


      if (errorEst < timeTol || fixFirstTimesteps || solverError
          || !(adaptInfo.getTimestep() > adaptInfo.getMinTimestep()) )
      {
        INFO(info, 6)("Accepted timestep at time = %e   with timestep = %e\n",
                      adaptInfo.getTime()  - adaptInfo.getTimestep(),
                      adaptInfo.getTimestep());

        adaptInfo.setTimestep(newTimestep);

        rejected = false;

        for (int i = 0; i < adaptInfo.getSize(); i++)
          INFO(info, 6)("time estimate for component %d = %e\n",
                        i, adaptInfo.getTimeEstSum(i));

        INFO(info, 6)("Combined time estimate = %e\n",
                      errorEst);

        if (errorEst > timeTol)
          MSG("Accepted timestep but tolerance TOL = %e not reached \n", timeTol);



      }
      else
      {
        for (int i = 0; i < adaptInfo.getSize(); i++)
          INFO(info, 6)("time estimate for component %d = %e\n",
                        i, adaptInfo.getTimeEstSum(i));

        INFO(info, 6)("Combined time estimate = %e\n",
                      errorEst);

        INFO(info, 6)("Rejected timestep at time = %e   with timestep = %e\n",
                      adaptInfo.getTime()  - adaptInfo.getTimestep(),
                      adaptInfo.getTimestep());

        adaptInfo.setTime(adaptInfo.getTime() - adaptInfo.getTimestep());
        adaptInfo.setTimestep(newTimestep);

        rejected = true;
      }

      if (firstTimestep)
        firstTimestep = false;

      if (fixFirstTimesteps > 0)
        fixFirstTimesteps--;
    }
    while (rejected ||
           (dbgTimestepStudy && (studyTimestep + 1 < static_cast<int>(dbgTimesteps.size()))));

    rosenbrockStat->acceptTimestep(adaptInfo);

    adaptInfo.setLastProcessedTimestep(adaptInfo.getTimestep());
    adaptInfo.incTimestepNumber();
  }


  double RosenbrockAdaptInstationary::getTimeEst(AdaptInfo& adaptInfo)
  {
    double errorEst = 0.0;
    for (int i = 0; i < adaptInfo.getSize(); i++)
    {
      double weight = adaptInfo.getTimeTolerance(i) + rosenbrockStat->getSolution(i)->L2Norm() * adaptInfo.getTimeRelativeTolerance(i);
      errorEst += sqr(adaptInfo.getTimeEstCombined(i)) / sqr(weight);
    }

    errorEst = std::sqrt(errorEst / adaptInfo.getSize());
    adaptInfo.setTimeEst(errorEst);
    return errorEst;
  }

}
