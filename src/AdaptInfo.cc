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


#include <string>
#include "AdaptInfo.h"

namespace AMDiS {

  using boost::lexical_cast;

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

}
