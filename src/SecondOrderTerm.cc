#include "SecondOrderTerm.h"
#include "ElInfo.h"
#include "MatrixVector.h"
#include "MatrixVectorOperations.h"

namespace AMDiS 
{
  void SecondOrderTerm::lalt(const DimVec<WorldVector<double> >& Lambda,
			     const WorldMatrix<double>& matrix,
			     mtl::dense2D<double>& LALt,
			     bool symm,
			     double factor) const
  {
    int dim = num_rows(LALt);
  
    if (symm) {
      for (int i = 0; i < dim; i++) {
	double val = 0.0;
	for (int k = 0; k < dimOfWorld; k++)
	  for (int l = 0; l < dimOfWorld; l++)
	    val += Lambda[i][k] * matrix[k][l] * Lambda[i][l];
	val *= factor;
	LALt[i][i] += val;
	for (int j = i + 1; j < dim; j++) {
	  val = 0.0;
	  for (int k = 0; k < dimOfWorld; k++)
	    for (int l = 0; l < dimOfWorld; l++)
	      val += Lambda[i][k] * matrix[k][l] * Lambda[j][l];
	  val *= factor;
	  LALt[i][j] += val;
	  LALt[j][i] += val;
	}
      }    
    } else {	
      for (int i = 0; i < dim; i++) {
	for (int j = 0; j < dim; j++) {
	  double val = 0.0;
	  for (int k = 0; k < dimOfWorld; k++)
	    for (int l = 0; l < dimOfWorld; l++) 
	      val += Lambda[i][k] * matrix[k][l] * Lambda[j][l];
	  val *= factor;
	  LALt[i][j] += val;
	}
      }    
    }      
  }


  void SecondOrderTerm::lalt_kl(const DimVec<WorldVector<double> >& Lambda,
				int k, int l,
				mtl::dense2D<double>& LALt,
				double factor)
  {
    int dim = num_rows(LALt);

    for (int i = 0; i < dim; i++)
      for (int j = 0; j < dim; j++)
	LALt[i][j] += factor * Lambda[i][k] * Lambda[j][l];
  }
  
  
  void SecondOrderTerm::l1lt(const DimVec<WorldVector<double> >& Lambda, 
			     mtl::dense2D<double>& LALt,
			     double factor) const
  {
    const int dim = num_rows(LALt);
    
    for (int i = 0; i < dim; i++) {
      double val = 0.0;
      for (int k = 0; k < dimOfWorld; k++)
	val += Lambda[i][k] * Lambda[i][k];
      val *= factor;
      LALt[i][i] += val;
      for (int j = i + 1; j < dim; j++) {
	val = 0.0;
	for (int k = 0; k < dimOfWorld; k++)
	  val += Lambda[i][k] * Lambda[j][k];
	val *= factor;
	LALt[i][j] += val;
	LALt[j][i] += val;
      }    
    }      
  }

} // end namespace AMDiS
