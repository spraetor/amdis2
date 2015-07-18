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
  
} // end namespace AMDiS
