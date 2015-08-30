#include "SubElInfo.h"
#include "ElInfo.h"
#include "Mesh.h"

namespace AMDiS
{
  SubElInfo::SubElInfo(VectorOfFixVecs<DimVec<double>>* lambda_,
                       const ElInfo* elInfo_)
    : elInfo(elInfo_)
  {
    FUNCNAME("SubElInfo::SubElInfo");

    int nPoints = lambda_->getSize();
    int dim = elInfo_->getMesh()->getDim();

    TEST_EXIT(nPoints == dim + 1)
    ("invalid number of vertices of subelement\n");

    FixVec<WorldVector<double>, VERTEX> worldCoords(dim);

    lambda = new VectorOfFixVecs<DimVec<double>>(dim, dim + 1, NO_INIT);
    for (int i = 0; i <= dim; i++)
    {
      (*lambda)[i] = (*lambda_)[i];
    }

    /**
     * Get worldcoordinates of the vertices of subelement in order to
     * calculate the corresponding determinant.
     */
    for (int i = 0; i < nPoints; i++)
    {
      elInfo->coordToWorld((const DimVec<double>)((*lambda)[i]),
                           worldCoords[i]);
    }

    det = elInfo->calcDet(worldCoords);
  }

} // end namespace AMDiS
