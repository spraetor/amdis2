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

// Written by Simon Praetorius


#ifndef ITL_PREONLY_INCLUDE
#define ITL_PREONLY_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>

namespace itl {

  /// Solver that simply applies a preconditioner to the rhs vector
  template < typename Matrix, typename Vector, typename Preconditioner, typename Iteration >
  int preonly(const Matrix &A, Vector &x, const Vector &b, const Preconditioner &P, Iteration& iter)
  {
    if (size(b) == 0)
      throw mtl::logic_error("empty rhs vector");
    
    typedef typename mtl::Collection<Vector>::value_type Scalar;
        
    // simple richardson iteration
    Vector r(b - A*x);
    Scalar res = two_norm(r);
    for (; !iter.finished(res); ++iter) {      
      x += Vector(solve(P, r));
      r = b - A*x;
      res = two_norm(r);
    }
    return iter;
  }

} // namespace itl

#endif // ITL_PREONLY_INCLUDE

