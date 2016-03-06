#pragma once

namespace AMDiS
{
  struct LazyOperatorTermBase
  {
    int getDegree() const
    {
      return 0;
    }

//     virtual ~LazyOperatorTermBase() { MSG("~LazyOperatorTermBase()"); }
  };

} // end namespace AMDiS
