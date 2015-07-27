#include "ScalableQuadrature.h"
#include "SubElInfo.h"

namespace AMDiS 
{
  ScalableQuadrature::ScalableQuadrature(Quadrature *quadrature)
    : Quadrature(*quadrature)
  {
    // Change name of quadrature.
    name = "Scalable" + getName();

    // copy quadrature->lambda to oldLambda
    oldLambda = new VectorOfFixVecs<DimVec<double> >(*(quadrature->getLambda()));

    // First assume all quadrature points to be valid.
    valid.resize(n_points);
    for (int i = 0; i < n_points; i++)
      valid[i] = true;
  }
  

  void ScalableQuadrature::scaleQuadrature(const SubElInfo& subElInfo)
  {
    //******************************************************************************
    // Manipulates the quadrature points for the assemblage of a subelement.    
    // A quadrature point for the assemblage of the subelement is then given in
    // barycentric coordinates with respect to the element containing the
    // subelement.
    // Thus the assemblage of the subelement is done by assembling the element
    // using the quadrature with manipulated quadrature points.
    // Note: The result must be corrected by the determinant of the subelement !!!
    //******************************************************************************

    /** 
     * If (l_0, l_1, l_2) is the old quadrature point and x_0, x_1, x_2 are the
     * barycentric coordinates of the vertices of the subelement, the new quadrature
     * point (n_0, n_1, n_2) is given as follows:
     * (n_0, n_1, n_2) = l_0 * x_0  +  l_1 * x_1  +  l_2 * x_2 .
     */
    for (int iq = 0; iq < n_points; iq++) {
      for (int i = 0; i <= dim; i++) {
      	/** 
      	 * Calculate the i-th component of the iq-th new quadrature point.
      	 */
      	double l = 0.0;
      
      	for (int j = 0; j <= dim; j++)
      	  l +=  getOldLambda(iq, j) * subElInfo.getLambda(j, i);
      
      	(*lambda)[iq][i] = l;
      
      	// If one component is less than zero, the scaled quadrature point
      	// is not inside the target element and is therefore set to be invalid.	
      	if (l < 0.0)
      	  valid[iq] = false;
      }
    }
  }

  
  void ScalableQuadrature::scaleQuadrature(DimMat<double> *scalMat)
  {
    for (int iq = 0; iq < n_points; iq++) {
      (*lambda)[iq].multMatrixVec(*scalMat, *getOldLambda(iq));

      for (int i = 0; i <= dim; i++) {
      	// If one component is less than zero, the scaled quadrature point
      	// is not inside the target element and is therefore set to be invalid.
      	if ((*lambda)[iq][i] < 0.0) {
      	  valid[iq] = false;
      	  break;
      	}
      }
    }
  }

} // end namespace AMDiS
