/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/



/** \file PeriodicInfo.h */

#ifndef AMDIS_PERIODICINFO_H
#define AMDIS_PERIODICINFO_H

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

}

#endif
