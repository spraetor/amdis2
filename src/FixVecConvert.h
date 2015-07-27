/** \file FixVecConvert.h */

#pragma once

#include "Global.h"

namespace AMDiS 
{

  template <class T, GeoIndex d1, GeoIndex d2>
  class VecConv
  {
  public:
    static FixVec<T,d1>& convertVec(FixVec<T,d2>& rhs, Mesh* mesh) 
    {
      TEST_EXIT(mesh->getGeo(d1) == mesh->getGeo(d2))("Incompatible dimensions!\n");
      return reinterpret_cast<FixVec<T,d1>&>(rhs);
    }
  };

} // end namespace AMDiS
