#include "SubQuadrature.h"
#include "Global.h"
#include "FixVec.h"

namespace AMDiS 
{
  SubQuadrature::SubQuadrature(Quadrature *quad, int dim_)
    : Quadrature((quad->getName() + " sub").c_str(),
		  quad->getDegree(),
		  dim_,
		  quad->getNumPoints(),
		  NULL,
		  quad->getWeight()),
      quad_(quad),
      subDim_(quad_->getDim())
  {
    if (dim == subDim_) 
      lambda = quad_->getLambda();
  }
  

  void SubQuadrature::scaleQuadrature(VectorOfFixVecs<DimVec<double> > &coords) 
  {
    if (!lambda || (lambda == quad_->getLambda()))
      lambda = new VectorOfFixVecs<DimVec<double> >(dim, n_points, NO_INIT);

    TEST_EXIT_DBG(coords.getSize() == (subDim_ + 1))
      ("invalid number of scale coords or invalid quadrature dim\n");
    TEST_EXIT_DBG(coords.getSizeOfFixVec() == (dim + 1))
      ("invalid coord dimension\n");

    for (int i = 0; i < n_points; i++) {
      DimVec<double> const& origin = quad_->getLambda(i);
      for (int j = 0; j < dim + 1; j++)
	(*lambda)[i][j] = 0.0;
      for (int j = 0; j < dim + 1; j++)
	for (int k = 0; k < subDim_ + 1; k++)
	  (*lambda)[i][j] += origin[k] * coords[k][j];
    }  
  }
  
} // end namespace AMDiS
