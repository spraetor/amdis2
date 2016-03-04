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


#ifndef ITL_GCR_INCLUDE
#define ITL_GCR_INCLUDE

#include <boost/numeric/mtl/concept/collection.hpp>

namespace itl
{

  template <typename Matrix, typename Vector, typename Preconditioner, typename Iteration>
  int gcr_full(const Matrix& A, Vector& x, const Vector& b, const Preconditioner& P, Iteration& iter)
  {
    using math::reciprocal;
    typedef typename mtl::Collection<Vector>::value_type Scalar;
    typedef typename mtl::Collection<Vector>::size_type Size;

    if (size(b) == 0)
      throw mtl::logic_error("empty rhs vector");

    Scalar zero= math::zero(b[0]), dbl_tol= 1.e-16; // TODO: tolerance as parameter or otherwise generalized

    Size k, kmax(std::min(size(x), Size(iter.max_iterations() - iter.iterations())));

    Vector r(b - A * x);
    Vector uk(size(x), zero), ck(size(x), zero);
    mtl::matrix::multi_vector<Vector> C(Vector(resource(x), zero), kmax);
    mtl::matrix::multi_vector<Vector> U(Vector(resource(x), zero), kmax);

    Scalar res(two_norm(r));
    Scalar alpha(zero), beta(zero), gamma(zero);

    for (k= 0; k < kmax && !iter.finished(res); ++k, ++iter)
    {
      uk = solve(P, r);
      // In order to avoid breakdown, use LSQR switch
      // Cite: C. Vuik, Further experiences with GMRESR, 1993
      //  if (two_norm(uk) < dbl_tol)
      //    uk = trans(A) * r; // requires transposed multiplication

      ck = A * uk;
      for (size_t i = 0; i < k; i++)
      {
        alpha = dot(C.vector(i), ck);

        ck -= alpha * C.vector(i);
        uk -= alpha * U.vector(i);
      }

      beta = two_norm(ck);
      if (beta < dbl_tol)
        return iter.fail(2, "search direction close to 0");

      C.vector(k) = ck * reciprocal(beta);
      U.vector(k) = uk * reciprocal(beta);

      gamma = dot(C.vector(k), r);
      x += U.vector(k) * gamma;
      r -= C.vector(k) * gamma;

      res = two_norm(r);
    }

    return iter;
  }

#if 0
  template <typename Matrix, typename Vector, typename LeftPreconditioner, typename RightPreconditioner, typename Iteration>
  int gmresr_trunc(const Matrix& A, Vector& x, const Vector& b, const LeftPreconditioner& L, const RightPreconditioner& R,
                   Iteration& iter, , typename mtl::Collection<Vector>::size_type nTrunc)
  {
TODO:
    implement truncated GCR solver, instead of restarted version
Cite:
    C. Vuik, Further experiences with GMRESR, 1993
  }
#endif

  /// Generalized Conjugate Residual method with restart
  /// Cite: S.C. Eisenstat, H.C. Elman, M.H. Schultz, Variational iterative methods for non symmetric systems of linear equations, 1983
  /// Cite: C. Vuik, New insight in GMRES-like methods with variable preconditioners
  template <typename Matrix, typename Vector, typename LeftPreconditioner,
            typename RightPreconditioner, typename Iteration>
  int gcr(const Matrix& A, Vector& x, const Vector& b, LeftPreconditioner& /*L*/, RightPreconditioner& R,
          Iteration& iter, typename mtl::Collection<Vector>::size_type restart)
  {
    do
    {
      Iteration inner(iter);
      inner.set_max_iterations(std::min(int(iter.iterations()+restart), iter.max_iterations()));
      inner.suppress_resume(true);
      gcr_full(A, x, b, R, inner);
      iter.update_progress(inner);
    }
    while (!iter.finished());

    return iter;
  }
} // namespace itl;

#endif // ITL_GMR_INCLUDE

