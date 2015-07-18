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

// Written by Simon Praetorius (adopted from previous implementation)


#ifndef ITL_FGMRES_HOUSEHOLDER_INCLUDE
#define ITL_FGMRES_HOUSEHOLDER_INCLUDE

#include <algorithm>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/multi_vector.hpp>
#include <boost/numeric/mtl/operation/givens.hpp>
#include <boost/numeric/mtl/operation/two_norm.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>

#include "solver/itl/details.hpp"

namespace itl {

/// Flexible Generalized Minimal Residual method (without restart) using householder othogonalization
/** It computes at most kmax_in iterations (or size(x) depending on what is smaller) 
    regardless on whether the termination criterion is reached or not.   **/
template < typename Matrix, typename Vector, typename RightPreconditioner, typename Iteration >
int fgmres_householder_full(const Matrix &A, Vector &x, const Vector &b,
               RightPreconditioner &R, Iteration& iter)
{
    using mtl::irange; using std::abs; using math::reciprocal; using mtl::iall; using mtl::imax;
    using mtl::signum; using mtl::vector::dot; using mtl::conj;
    typedef typename mtl::Collection<Vector>::value_type Scalar;
    typedef typename mtl::Collection<Vector>::size_type  Size;

    if (size(b) == 0) throw mtl::logic_error("empty rhs vector");

    const Scalar zero= math::zero(Scalar());
    Scalar       rho, bnrm2, temp, alpha, beta;
    Size         k, kmax(std::min(size(x), Size(iter.max_iterations() - iter.iterations())));
    Vector       w(solve(R, b)), r(b - A*x);
    mtl::matrix::multi_vector<Vector>   V(Vector(resource(x), zero), kmax+1); 
    mtl::matrix::multi_vector<Vector>   Z(Vector(resource(x), zero), kmax+1); 
    mtl::vector::dense_vector<Scalar>   sn(kmax, zero), cs(kmax, zero), s(kmax+1, zero), y(kmax, zero);  // replicated in distributed solvers 
    mtl::matrix::dense2D<Scalar>        H(kmax, kmax);

    bnrm2 = two_norm(b);
    if (bnrm2 == zero) {
      set_to_zero(x);
      // b == 0 => solution = 0
      return iter.terminate(bnrm2);
    }
    
    rho = two_norm(r);				// norm of preconditioned residual
    if (iter.finished(rho))			// initial guess is good enough solution
      return iter;
    
    // u = r + sing(r(0))*||r||*e0
    beta = signum(r[0])*rho;
    w = r;
    w[0] += beta;
    w *= reciprocal(two_norm(w));
    
    V.vector(0) = w;
    H = zero;
    s[0] = -beta;

    // GMRES iteration
    for (k= 0; k < kmax && !iter.finished(rho); ++k, ++iter) {
      
      w = (-2.0 * V.vector(k)[k])*V.vector(k);
      w[k] += 1.0;
      // v := P_0*...*P_{k-2}*(P_{k-1} * e_k)
      for (Size i= k; i > 0; i--) {
	temp = 2.0 * dot(V.vector(i-1), w);
	w -= temp * V.vector(i-1);
      }
      
      temp = two_norm(w);
      if (temp == zero)
	return iter.fail(2, "GMRES: breakdown");
      
      // Explicitly normalize v to reduce the effects of round-off.
      Z.vector(k) = solve(R, w);
      w = A * Z.vector(k);
      
      // P_{k-1}*...*P_0*Av
      for (Size i = 0; i <= k; i++) {
	temp = 2.0 * dot(V.vector(i), w);
	w -= temp * V.vector(i);
      }
      
      temp = two_norm(w);
      if (temp == zero)
	return iter.fail(3, "GMRES: breakdown");
      
      irange range_to_end(k+1,imax);
      set_to_zero(V.vector(k+1));
      V.vector(k+1)[range_to_end] = w[range_to_end];
      alpha = two_norm(V.vector(k+1));
      if (alpha != 0.0) {
	alpha *= signum(w[k+1]);
	V.vector(k+1)[k+1] += alpha;
	V.vector(k+1) *= reciprocal(two_norm(V.vector(k+1)));
	
	w[k+1] = -alpha;
      }
      
      for (Size i= 0; i < k; i++) {
	temp   =  conj(cs[i])*w[i] + conj(sn[i])*w[i+1];
	w[i+1] =- sn[i]*w[i] + cs[i]*w[i+1];
	w[i]   =  temp;
      }
      
      details::rotmat(w[k], w[k+1], cs[k], sn[k]);
      
      s[k+1] = -sn[k]*s[k];
      s[k]   = conj(cs[k])*s[k];
      w[k] = cs[k]*w[k] + sn[k]*w[k+1];
      w[k+1] = 0.0;
      
      irange range(num_rows(H));
      H[iall][k] = w[range];
      
      rho = std::abs(s[k+1]);
    }
        
    // reduce k, to get regular matrix
//     while (k > 0 && std::abs(s[k-1]) <= iter.atol()) k--;

    // iteration is finished -> compute x: solve H*y=s as far as rank of H allows
    irange range(k);
    for (; !range.empty(); --range) {
      try {
	  y[range] = upper_trisolve(H[range][range], s[range]); 
      } catch (mtl::matrix_singular) { continue; } // if singular then try with sub-matrix
      break;
    }

    if (range.finish() < k)
  	std::cerr << "GMRES orhogonalized with " << k << " vectors but matrix singular, can only use " 
		  << range.finish() << " vectors!\n";
    if (range.empty())
        return iter.fail(3, "GMRES did not find any direction to correct x");
    

    x += Z.vector(range) * y[range];
    
    r = b - A*x;
    return iter.terminate(r);
}

/// Felxible Generalized Minimal Residual method with restart
template < typename Matrix, typename Vector,
           typename RightPreconditioner, typename Iteration >
int fgmres_householder(const Matrix &A, Vector &x, const Vector &b,
          RightPreconditioner &R,
	  Iteration& iter, typename mtl::Collection<Vector>::size_type restart)
{   
  do {
    Iteration inner(iter);
    inner.set_max_iterations(std::min(int(iter.iterations()+restart), iter.max_iterations()));
    inner.suppress_resume(true);
    fgmres_householder_full(A, x, b, R, inner);
    iter.update_progress(inner);
  } while (!iter.finished());

  return iter;
}



} // namespace itl

#endif // ITL_FGMRES_HOUSEHOLDER_INCLUDE


