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


#ifndef ITL_DETAIL_INCLUDE
#define ITL_DETAIL_INCLUDE

#include <boost/math/special_functions/sign.hpp>

namespace itl {
  
  namespace details {
    
    /// Compute the Givens rotation matrix parameters for a and b.
//     template<typename T>
//     void rotmat(const T& a, const T& b , T& c, T& s)
//     {
//       using std::abs; using std::sqrt; using mtl::conj;
//       
//       const T zero = math::zero(T());
//       if (a == zero) {
// 	c = 0.0;
// 	s = 1.0;
//       } else {
// 	double temp = abs(a) / sqrt( conj(a)*a + conj(b)*b );
// 	c = temp;
// 	s = temp * (b / a);
//       }
//     }
    
    inline void rotmat(const double& a, const double& b , double& c, double& s)
    {
      using std::abs; using std::sqrt;
      if ( b == 0.0 ) {
	  c = 1.0;
	  s = 0.0;
      } else if ( abs(b) > abs(a) ) {
	  double temp = a / b;
	  s = 1.0 / sqrt( 1.0 + temp*temp );
	  c = temp * s;
      } else {
	  double temp = b / a;
	  c = 1.0 / sqrt( 1.0 + temp*temp );
	  s = temp * c;
      }
    }
  }
}

#endif // ITL_DETAIL_INCLUDE
