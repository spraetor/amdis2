#include "SurfaceQuadrature.hpp"

#include <algorithm>

#include "Quadrature.hpp"

namespace AMDiS
{
  SurfaceQuadrature::SurfaceQuadrature(Quadrature* q,
                                       VectorOfFixVecs<DimVec<double>> const& c)
    : Quadrature((q->getName() + " surface").c_str(),
                 q->getDegree(),
                 q->getDim() + 1,
                 q->getNumPoints(),
                 NULL,
                 q->getWeight()),
      quad(q),
      coords(c)
  {
    lambda = new VectorOfFixVecs<DimVec<double>>(dim, n_points, NO_INIT);

    // for each integration point
    for (int i = 0; i < n_points; i++)
    {
      // get coords of quadrature point in dim-1
      DimVec<double> const& origin = quad->getLambda(i);

      for (int j = 0; j < dim + 1; j++)
        (*lambda)[i][j] = 0.0;

      for (int j = 0; j < dim; j++)
        for (int k = 0; k < dim + 1; k++)
          (*lambda)[i][k] += origin[j] * coords[j][k];
    }
  }

  void SurfaceQuadrature::scaleSurfaceQuadrature(VectorOfFixVecs<DimVec<double>>& c)
  {
    // copy coords NOTE: why???
    for (int i = 0; i < dim; i++)
      coords[i] = c[i];

    // for each integration point
    for (int i = 0; i < n_points; i++)
    {
      // get coords of quadrature point in dim-1
      DimVec<double> const& origin = quad->getLambda(i);

      for (int j = 0; j < dim + 1; j++)
        (*lambda)[i][j] = 0.0;

      for (int j = 0; j < dim; j++)
        for (int k = 0; k < dim+1; k++)
          (*lambda)[i][k] += origin[j] * coords[j][k];
    }
  }

} // end namespace AMDiS
