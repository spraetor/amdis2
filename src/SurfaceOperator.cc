#include "SurfaceOperator.h"
#include "SurfaceQuadrature.h"
#include "FiniteElemSpace.h"
#include "Mesh.h"

namespace AMDiS
{
  SurfaceOperator::SurfaceOperator(Operator const& operat,
                                   VectorOfFixVecs<DimVec<double>>& coords)
    : Operator(operat),
      coords_(coords),
      quad2(NULL),
      quad1GrdPsi(NULL),
      quad1GrdPhi(NULL),
      quad0(NULL)
  {
    assembler = NULL;

    int dim = rowFeSpace->getMesh()->getDim();
    int degree;

    // create surface quadratures
    if (secondOrder.size() > 0)
    {
      degree = getQuadratureDegree(2);
      quad2 =
        new SurfaceQuadrature(Quadrature::provideQuadrature(dim - 1, degree), coords);
    }

    if (firstOrderGrdPsi.size() > 0)
    {
      degree = getQuadratureDegree(1, GRD_PSI);
      quad1GrdPsi =
        new SurfaceQuadrature(Quadrature::provideQuadrature(dim - 1, degree), coords);
    }

    if (firstOrderGrdPhi.size() > 0)
    {
      degree = getQuadratureDegree(1, GRD_PHI);
      quad1GrdPhi =
        new SurfaceQuadrature(Quadrature::provideQuadrature(dim - 1, degree), coords);
    }

    if (zeroOrder.size() > 0)
    {
      degree = getQuadratureDegree(0);
      quad0 =
        new SurfaceQuadrature(Quadrature::provideQuadrature(dim - 1, degree), coords);
    }

    // initialize assembler with surface quadratures
    optimized = false;
    initAssembler(quad2, quad1GrdPsi, quad1GrdPhi, quad0);
  }


  void SurfaceOperator::adaptSurfaceOperator(VectorOfFixVecs<DimVec<double>>& coords)
  {
    coords_ = coords;

    if (quad2)
      quad2->scaleSurfaceQuadrature(coords);
    if (quad1GrdPsi)
      quad1GrdPsi->scaleSurfaceQuadrature(coords);
    if (quad1GrdPhi)
      quad1GrdPhi->scaleSurfaceQuadrature(coords);
    if (quad0)
      quad0->scaleSurfaceQuadrature(coords);
  }


  void SurfaceOperator::getElementMatrix(ElInfo const* elInfo,
                                         ElementMatrix& userMat,
                                         double factor = 1.0)
  {
    int dim = rowFeSpace->getMesh()->getDim();
    double origDet = elInfo->getDet();

    FixVec<WorldVector<double>, VERTEX> worldCoords(dim-1);

    // transform barycentric coords to world coords
    for (int i = 0; i < dim; i++)
      elInfo->coordToWorld(coords_[i], worldCoords[i]);

    // set determinant for world coords of the side
    const_cast<ElInfo*>(elInfo)->setDet(elInfo->calcDet(worldCoords));

    // calc element matrix
    Operator::getElementMatrix(elInfo, userMat, factor);

    // set determinant for world coords of the side
    const_cast<ElInfo*>(elInfo)->setDet(origDet);
  }


  void SurfaceOperator::getElementVector(ElInfo const* elInfo,
                                         DenseVector<double>& userVec,
                                         double factor = 1.0)
  {
    int dim = rowFeSpace->getMesh()->getDim();
    double origDet = elInfo->getDet();

    FixVec<WorldVector<double>, VERTEX> worldCoords(dim-1);

    // transform barycentric coords to world coords
    for (int i = 0; i < dim; i++)
      elInfo->coordToWorld(coords_[i], worldCoords[i]);

    // set determinant for world coords of the side
    const_cast<ElInfo*>(elInfo)->setDet(elInfo->calcDet(worldCoords));

    // calc element vector
    Operator::getElementVector(elInfo, userVec, factor);

    const_cast<ElInfo*>(elInfo)->setDet(origDet);
  }

} // end namespace AMDiS
