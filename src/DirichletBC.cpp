#include "DirichletBC.hpp"
#include "ElInfo.hpp"
#include "BasisFunction.hpp"
#include "DOFVector.hpp"
#include "DOFMatrix.hpp"

namespace AMDiS
{
  namespace detail
  {
    void DirichletBC::fillBoundaryCondition(DOFMatrix* DBG_VAR( matrix ),
                                            ElInfo* /*elInfo*/,
                                            const DegreeOfFreedom* /*dofIndices*/,
                                            const BoundaryType* /*localBound*/,
                                            int /*nBasFcts*/)
    {
      FUNCNAME_DBG("DirichletBC::fillBoundaryCondition()");
      TEST_EXIT_DBG(matrix->getRowFeSpace() == rowFeSpace)("invalid row fe space\n");
    }

    void DirichletBC::initVector(DOFVectorBase<double>* vec)
    {
      if (dynamic_cast<DOFVector<double>*>(vec))
        dynamic_cast<DOFVector<double>*>(vec)->getDirichletValues().clear();
    }

  } // end namespace detail

} // end namespace AMDiS
