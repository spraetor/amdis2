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



/** \file FixVecConvert.h */

#ifndef AMDIS_FIXVECCONVERT_H_
#define AMSID_FIXVECCONVERT_H_

#include "Global.h"

namespace AMDiS {

  template<typename T,GeoIndex d1,GeoIndex d2>
  class VecConv
  {
  public:
    static FixVec<T,d1>& convertVec(FixVec<T,d2>& rhs, Mesh* mesh) 
    {
      TEST_EXIT(mesh->getGeo(d1) == mesh->getGeo(d2))("Incompatible dimensions!\n");
      return reinterpret_cast<FixVec<T,d1>&>(rhs);
    }
  };

}

#endif
