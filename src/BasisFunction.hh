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


#include "BasisFunction.h"

namespace AMDiS {

  template <class T>
  T BasisFunction::evalUh(const DimVec<double>& lambda,
			  const mtl::dense_vector<T>& uh_loc) const
  {
    T val;
    nullify(val);
    
    for (int i = 0; i < nBasFcts; i++)
      val += uh_loc[i] * (*(*phi)[i])(lambda);

    return val;
  }

  template <class T>
  typename GradientType<T>::type& 
  BasisFunction::evalGrdUh(const DimVec<double>& lambda,
			   const DimVec<WorldVector<double> >& grd_lambda,
			   const mtl::dense_vector<T>& uh_loc,
			   typename GradientType<T>::type& val) const
  {

    mtl::dense_vector<double> grdTmp1(dim + 1);
    T null;
    nullify(null);
    mtl::dense_vector<T> grdTmp2(dim + 1, null);

    for (int i = 0; i < nBasFcts; i++) {
      (*(*grdPhi)[i])(lambda, grdTmp1);

      for (int j = 0; j < dim + 1; j++)
	grdTmp2[j] += uh_loc[i] * grdTmp1[j];
    }
    
    nullify(val);
    for (int i = 0; i < dow; i++) {
      for (int j = 0; j < dim + 1; j++)
	val[i] += grd_lambda[j][i] * grdTmp2[j];
    }

    return val;
  }
}
