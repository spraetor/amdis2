#include <string>

#include "AdaptInfo.h"

namespace AMDiS 
{
  
  AdaptInfo::ScalContent::ScalContent(std::string prefix)
  	: est_sum(0.0),
  	  est_t_sum(0.0),
  	  est_max(0.0),
  	  est_t_max(0.0),
  	  fac_max(0.0),
  	  fac_sum(1.0),
  	  spaceTolerance(0.0),
  	  timeTolerance(0.0),
  	  timeRelativeTolerance(0.0),
  	  timeErrLow(0.0),
  	  coarsenAllowed(0),
  	  refinementAllowed(1),
  	  refineBisections(1),
  	  coarseBisections(1)	  	
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
      	spaceIteration(-1),
      	maxSpaceIteration(-1),
      	timestepIteration(0),
      	maxTimestepIteration(30),
      	timeIteration(0),
      	maxTimeIteration(30),
      	time(0.0),
      	startTime(0.0),
      	endTime(1.0),
      	timestep(0.0),
      	lastProcessedTimestep(0.0),
      	minTimestep(0.0),
      	maxTimestep(1.0),
      	timestepNumber(0),
      	nTimesteps(0),
      	solverIterations(0),
      	maxSolverIterations(0),
      	solverTolerance(1e-8),
      	solverResidual(0.0),
      	globalTimeTolerance(1.0),
        scalContents(size),
      	deserialized(false),
      	rosenbrockMode(false)
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
    
    for (int i = 0; i < size; i++) {
    	scalContents[i] = new ScalContent(name + "[" + std::to_string(i) + "]");  
    }
  }
  
  
  void AdaptInfo::setScalContents(int newSize) 
  {
    int oldSize = static_cast<int>(scalContents.size());

    if (newSize > oldSize) { 
      scalContents.resize(newSize);

      for (int i = oldSize; i < newSize; i++)
	scalContents[i] = 
	  new ScalContent(name + "[" + std::to_string(i) + "]"); 
    }
  }
  
  
  void AdaptInfo::printTimeErrorLowInfo() const
  {
    for (size_t i = 0; i < scalContents.size(); i++){
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
