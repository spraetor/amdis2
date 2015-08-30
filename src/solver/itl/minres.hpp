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

// Written by Thomas Witkowski


#ifndef AMDIS_ITL_MINRES_INCLUDE
#define AMDIS_ITL_MINRES_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>

namespace itl
{

  /// Minimal Residual method
  template <typename Matrix, typename Vector,
            typename LeftPreconditioner, typename RightPreconditioner,
            typename Iteration>
  int minres(const Matrix& A, Vector& x, const Vector& b,
             const LeftPreconditioner& L, const RightPreconditioner& R,
             Iteration& iter)
  {
    using std::abs;
    using math::reciprocal;
    typedef typename mtl::Collection<Vector>::value_type Scalar;

    if (size(b) == 0)
      throw mtl::logic_error("empty rhs vector");

    Scalar                zero= math::zero(b[0]), one= math::one(b[0]);
    Vector v0(size(x), zero), v1(b - A * x), v2(v1), z1(solve(L, v1)), z2(size(x), zero);
    Vector w0(size(x), zero), w1(size(x), zero), w2(size(x), zero);

    Scalar s0(zero), s1(zero), c0(one), c1(one), gamma0(one);
    Scalar gamma1(sqrt(dot(z1, v1))), gamma2(zero), eta(gamma1);
    Scalar sigma1(one), alpha0(zero), alpha1(zero), alpha2(zero), alpha3(zero);

    while (!iter.finished(abs(eta)))
    {
      z1 *= reciprocal(gamma1);
      v2 = A * z1;
      sigma1 = dot(v2, z1);
      v2 += -(sigma1 / gamma1) * v1 - (gamma1 / gamma0) * v0;

      z2 = solve(L, v2);

      gamma2 = sqrt(dot(z2, v2));
      alpha0 = c1 * sigma1 - c0 * s1 * gamma1;
      alpha1 = sqrt(alpha0 * alpha0 + gamma2 * gamma2);
      alpha2 = s1 * sigma1 + c0 * c1 * gamma1;
      alpha3 = s0 * gamma1;

      c0 = c1;
      c1 = alpha0 / alpha1;
      s0 = s1;
      s1 = gamma2 / alpha1;

      w2 = z1 - alpha3 * w0 - alpha2 * w1;
      w2 *=  reciprocal(alpha1);

      x += c1 * eta * w2;
      eta *= -s1;

      w0 = w1;
      w1 = w2;
      v0 = v1;
      v1 = v2;
      z1 = z2;

      gamma0 = gamma1;
      gamma1 = gamma2;

      ++iter;
    }

    return iter;
  }

} // namespace itl;

#endif // AMDIS_ITL_MINRES_INCLUDE
