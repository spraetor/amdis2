#include <string>
#include <iostream>

#include "AdaptInfo.h"
#include "Initfile.h"

namespace AMDiS
{

  AdaptInfo::ScalContent::ScalContent(std::string prefix)
  {
    Parameters::get(prefix + "->tolerance", spaceTolerance);
    Parameters::get(prefix + "->time tolerance", timeTolerance);
    Parameters::get(prefix + "->time relative tolerance", timeRelativeTolerance);
    Parameters::get(prefix + "->coarsen allowed", coarsenAllowed);
    Parameters::get(prefix + "->refinement allowed", refinementAllowed);
    Parameters::get(prefix + "->refine bisections", refineBisections);
    Parameters::get(prefix + "->coarsen bisections", coarseBisections);
    Parameters::get(prefix + "->sum factor", fac_sum);
    Parameters::get(prefix + "->max factor", fac_max);

    if (timeTolerance == 0.0 && timeRelativeTolerance == 0.0)
      timeTolerance = 1.0;
    timeErrLow = timeTolerance * 0.3;
  }


  AdaptInfo::AdaptInfo(std::string name_, int size)
    : name(name_),
      scalContents(size)
  {
    // init();
    Parameters::get(name + "->start time", startTime);
    time = startTime;
    Parameters::get(name + "->timestep", timestep);
    Parameters::get(name + "->end time", endTime);
    Parameters::get(name + "->max iteration", maxSpaceIteration);
    Parameters::get(name + "->max timestep iteration", maxTimestepIteration);
    Parameters::get(name + "->max time iteration", maxTimeIteration);
    Parameters::get(name + "->min timestep", minTimestep);
    Parameters::get(name + "->max timestep", maxTimestep);
    Parameters::get(name + "->number of timesteps", nTimesteps);
    Parameters::get(name + "->time tolerance", globalTimeTolerance);

    for (int i = 0; i < size; i++)
    {
      scalContents[i] = new ScalContent(name + "[" + std::to_string(i) + "]");
    }
  }


  void AdaptInfo::setScalContents(int newSize)
  {
    int oldSize = static_cast<int>(scalContents.size());

    if (newSize > oldSize)
    {
      scalContents.resize(newSize);

      for (int i = oldSize; i < newSize; i++)
        scalContents[i] =
          new ScalContent(name + "[" + std::to_string(i) + "]");
    }
  }


  void AdaptInfo::printTimeErrorLowInfo() const
  {
    for (size_t i = 0; i < scalContents.size(); i++)
    {
      std::cout << "    Time error estimate     ["<<i<<"] = "
                << getTimeEstCombined(i) << "\n"
                << "    Time error estimate sum ["<<i<<"] = "
                << scalContents[i]->est_t_sum << "\n"
                << "    Time error estimate max ["<<i<<"] = "
                << scalContents[i]->est_t_max << "\n"
                << "    Time error low bound    ["<<i<<"] = "
                << scalContents[i]->timeErrLow << "\n"
                << "    Time error high bound   ["<<i<<"] = "
                << scalContents[i]->timeTolerance << "\n";
    }
  }

  void AdaptInfo::reset()
  {
    spaceIteration = -1;
    timestepIteration = 0;
    timeIteration = 0;
    time = 0.0;
    timestep = 0.0;
    timestepNumber = 0;
    solverIterations = 0;
    solverResidual = 0.0;

    Parameters::get(name + "->timestep", timestep);
    lastProcessedTimestep=timestep;
  }

} // end namespace AMDiS
