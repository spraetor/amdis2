#include "BasisFunction.hpp"

namespace AMDiS
{

  template <class T>
  T BasisFunction::evalUh(const DimVec<double>& lambda,
                          const DenseVector<T>& uh_loc) const
  {
    T val;
    nullify(val);

    for (int i = 0; i < nBasFcts; i++)
      val += uh_loc[i] * (*(*phi)[i])(lambda);

    return val;
  }

  template <class T> Gradient_t<T>&
  BasisFunction::evalGrdUh(const DimVec<double>& lambda,
                           const DimVec<WorldVector<double>>& grd_lambda,
                           const DenseVector<T>& uh_loc,
                           Gradient_t<T>& val) const
  {

    DenseVector<double> grdTmp1(dim + 1);
    T null;
    nullify(null);
    DenseVector<T> grdTmp2(dim + 1, null);

    for (int i = 0; i < nBasFcts; i++)
    {
      (*(*grdPhi)[i])(lambda, grdTmp1);

      for (int j = 0; j < dim + 1; j++)
        grdTmp2[j] += uh_loc[i] * grdTmp1[j];
    }

    nullify(val);
    for (int i = 0; i < dow; i++)
    {
      for (int j = 0; j < dim + 1; j++)
        val[i] += grd_lambda[j][i] * grdTmp2[j];
    }

    return val;
  }

} // end namespace AMDiS
