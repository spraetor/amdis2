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


#ifndef ITL_FGMRES_INCLUDE
#define ITL_FGMRES_INCLUDE

#include <algorithm>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/matrix/dense2D.hpp>
#include <boost/numeric/mtl/matrix/multi_vector.hpp>
#include <boost/numeric/mtl/operation/givens.hpp>
#include <boost/numeric/mtl/operation/two_norm.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/irange.hpp>

#include "details.hpp"

namespace itl {
  
/// Flexible Generalized Minimal Residual method (without restart)
/// Cite: Youcef Saad, A Flexible Inner-Outer Preconditioned GMRES Algorithm, 1993
/** It computes at most kmax_in iterations (or size(x) depending on what is smaller) 
    regardless on whether the termination criterion is reached or not.   **/
template < typename Matrix, typename Vector, typename Preconditioner, typename Iteration >
int fgmres_full(const Matrix &A, Vector &x, const Vector &b,
               Preconditioner &P, Iteration& iter)
{
    using mtl::irange; using mtl::iall; using std::abs; using std::sqrt; using math::reciprocal;
    typedef typename mtl::Collection<Vector>::value_type Scalar;
    typedef typename mtl::Collection<Vector>::size_type  Size;

    if (size(b) == 0) throw mtl::logic_error("empty rhs vector");
    
    const Scalar                zero= math::zero(Scalar()), breakdown_tol= 1.e-16, kappa = 10.0;
    Scalar                      rho, bnrm2, temp, hr;
    Size                        k, kmax(std::min(size(x), Size(iter.max_iterations() - iter.iterations())));
    Vector                      w(resource(x)), r(b - A*x);
    mtl::matrix::multi_vector<Vector>   V(Vector(resource(x), zero), kmax+1); 
    mtl::matrix::multi_vector<Vector>   Z(Vector(resource(x), zero), kmax+1); 
    mtl::matrix::dense2D<Scalar>        H(kmax+1, kmax);
    
    mtl::vector::dense_vector<Scalar>   sn(kmax, zero), cs(kmax, zero), s(kmax+1, zero), y(kmax, zero);  // replicated in distributed solvers 

    bnrm2 = two_norm(b);
    if (bnrm2 == zero) {
      set_to_zero(x);
      // b == 0 => solution = 0
      return iter.terminate(bnrm2);
    }
    
    rho = two_norm(r);				// norm of preconditioned residual
    if (iter.finished(rho))			// initial guess is good enough solution
      return iter;
    
    V.vector(0) = r * reciprocal(rho);
    H = zero;
    s[0] = rho;

    // FGMRES iteration
    for (k= 0; k < kmax && !iter.finished(rho) ; ++k, ++iter) {
      Z.vector(k) = solve(P, V.vector(k));
      w = A * Z.vector(k);
      temp = two_norm(w);
      
      for (Size j= 0; j < k+1; j++) {
	H[j][k] = dot(V.vector(j), w);
	w -= H[j][k] * V.vector(j);
      }
      H[k+1][k]= two_norm(w);
            
      // reorthogonalization, only if "heuristic" condition is fulfilled
      if (H[k+1][k] < temp * reciprocal(kappa)) {
	for (Size i= 0; i < k+1; i++) {
	    hr = dot(w, V.vector(i));
            H[i][k] += hr;
            w -= hr * V.vector(i);
	}
        temp = two_norm(w);
	
	if (temp < H[k+1][k] * reciprocal(kappa)) {
	  set_to_zero(w);
	  H[k+1][k] = 0.0;
	} else {
	  H[k+1][k] = temp;
	}
      }
      
      if (H[k+1][k] < breakdown_tol)
	return iter.fail(2, "FGMRES: Singular matrix - nearly hard breakdown");
      
      V.vector(k+1) = w * reciprocal(H[k+1][k]);

      // k Given's rotations
      for(Size i= 0; i < k; i++) {
	temp =        cs[i]*H[i][k] + sn[i]*H[i+1][k];
	H[i+1][k] = - sn[i]*H[i][k] + cs[i]*H[i+1][k];
	H[i][k] = temp;
      }
      
      details::rotmat(H[k][k], H[k+1][k], cs[k], sn[k]);
      
      s[k+1] = -sn[k]*s[k];
      s[k]   =  cs[k]*s[k];
      H[k][k] = cs[k]*H[k][k] + sn[k]*H[k+1][k];
      H[k+1][k] = 0.0;
      
      rho = std::abs(s[k+1]);
    }
    
    // reduce k, to get regular matrix
//     while (k > 0 && abs(s[k-1]) <= iter.atol()) k--;

    // iteration is finished -> compute x: solve H*y=g as far as rank of H allows
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
    
//     r = b - A*x;
    return iter.terminate(rho);
}

/// Flexible Generalized Minimal Residual method with restart
template < typename Matrix, typename Vector, typename LeftPreconditioner,
           typename RightPreconditioner, typename Iteration >
int fgmres(const Matrix &A, Vector &x, const Vector &b,
          LeftPreconditioner &L, RightPreconditioner &R,
	  Iteration& iter, typename mtl::Collection<Vector>::size_type restart)
{   
    do {
	Iteration inner(iter);
	inner.set_max_iterations(std::min(int(iter.iterations()+restart), iter.max_iterations()));
	inner.suppress_resume(true);
	fgmres_full(A, x, b, R, inner);
	iter.update_progress(inner);
    } while (!iter.finished());

    return iter;
}



} // namespace itl

#endif // ITL_FGMRES_INCLUDE


