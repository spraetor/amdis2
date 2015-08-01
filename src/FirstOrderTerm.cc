#include "FirstOrderTerm.h"
#include "ElInfo.h"
#include "traits/num_rows.hpp"
#include "traits/size.hpp"
#include "MatrixVector.h"
#include "MatrixVectorOperations.h"

namespace AMDiS 
{
  
  /// Implenetation of FirstOrderTerm::evalImpl().
  void FirstOrderTerm::evalImpl(int nPoints,
                        				const DenseVector<double>& uhAtQP,
                        				const DenseVector<WorldVector<double> >& grdUhAtQP,
                        				const DenseVector<WorldMatrix<double> >& D2UhAtQP,
                        				DenseVector<double>& result,
                        				double factor)  const
  {
    int dow = Global::getGeo(WORLD);

    if (num_rows(grdUhAtQP) > 0) {
      for (int iq = 0; iq < nPoints; iq++) {
      	double resultQP = 0.0;
      	for (int i = 0; i < dow; i++)
      	  resultQP += grdUhAtQP[iq][i];
      	
      	result[iq] += resultQP * factor;
      	// result[iq] += sum(grdUhAtQP[iq]) * factor;
      }
    }
  }


  void FirstOrderTerm::l1(const DimVec<WorldVector<double> >& Lambda,
			  DenseVector<double>& Lb,
			  double factor) const
  {
    const int dim = size(Lambda);

    for (int i = 0; i < dim; i++) {
      double val = 0.0;
      for (int j = 0; j < dimOfWorld; j++)
        val += Lambda[i][j];

      Lb[i] += val * factor;
      // Lb[i] += sum(Lambda[i]) * factor;
    }    
  }
  
  
  void FirstOrderTerm::lb(const DimVec<WorldVector<double> >& Lambda,
			  const WorldVector<double>& b,
			  DenseVector<double>& Lb,
			  double factor) const
  {
    const int dim = size(Lambda);

    for (int i = 0; i < dim; i++)
      Lb[i] += (Lambda[i] * b) * factor;
  }
  
  
  void FirstOrderTerm::lb_one(const DimVec<WorldVector<double> >& Lambda,
			      DenseVector<double>& Lb,
			      double factor) const
  {
    const int dim = size(Lambda);

    for (int i = 0; i < dim; i++)
      Lb[i] += Lambda[i][bOne] * factor;
  }
  
} // end namespace AMDiS
