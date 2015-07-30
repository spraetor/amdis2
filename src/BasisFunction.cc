#include <algorithm>

#include "FixVec.h"
#include "DOFVector.h"
#include "BasisFunction.h"
#include "Lagrange.h"

namespace AMDiS {


  /****************************************************************************/
  /*  Lagrangian basisfunctions of order 0-4; these                           */
  /*  functions are evaluated in barycentric coordinates; the derivatives     */
  /*  are those corresponding to these barycentric coordinates.               */
  /****************************************************************************/

  BasisFunction::BasisFunction(std::string name_, int dim_, int degree_)
    : name(name_), 
      degree(degree_), 
      dim(dim_)
  {
    nDOF = new DimVec<int>(dim, -1);
    dow = Global::getGeo(WORLD);
  }


  BasisFunction::~BasisFunction()
  {
    delete nDOF;
  }


  const WorldMatrix<double>& 
  BasisFunction::evalD2Uh(const DimVec<double>& lambda,
        						      const DimVec<WorldVector<double> >& grd_lambda,
        						      const ElementVector& uh_loc, 
        						      WorldMatrix<double>* D2_uh) const
  {
    // TODO: REMOVE STATIC
    static WorldMatrix<double> D2(DEFAULT_SIZE, DEFAULT_SIZE, 0.0);
    DimMat<double> D2_b(dim, 0.0);
    DimMat<double> D2_tmp(dim, 0.0);
    WorldMatrix<double> *val = D2_uh ? D2_uh : &D2;
  
    for (int i = 0; i < nBasFcts; i++) {
      (*(*d2Phi)[i])(lambda, D2_b);
      for (int k = 0; k < dim + 1; k++)
      	for (int l = 0; l < dim + 1; l++)
      	  D2_tmp[k][l] += uh_loc[i] * D2_b[k][l];
    }

    for (int i = 0; i < dow; i++)
      for (int j = 0; j < dow; j++) {
      	(*val)[i][j] = 0.0;
      	for (int k = 0; k < dim + 1; k++)
      	  for (int l = 0; l < dim + 1; l++)
      	    (*val)[i][j] += grd_lambda[k][i] * grd_lambda[l][j] * D2_tmp[k][l];
      }
    
    return ((*val));
  }


  int BasisFunction::getNumberOfDofs(int i) const
  { 
    return (*nDOF)[i];
  }

} // end namespace AMDiS
