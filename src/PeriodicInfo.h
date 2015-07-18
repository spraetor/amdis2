/** \file PeriodicInfo.h */

#pragma once

#include <map>

namespace AMDiS {

  class PeriodicInfo 
  {
  public:
    int mode;
      
    int type;

    int outputIndex;

    int neighIndex;

    std::map<int, int> vertexMap;
  };

} // end namespace AMDiS
