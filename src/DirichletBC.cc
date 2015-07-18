#include "DirichletBC.h"
#include "ElInfo.h"
#include "BasisFunction.h"
#include "DOFVector.h"
#include "DOFMatrix.h"

namespace AMDiS 
{
  namespace detail 
  {
      void DirichletBC::fillBoundaryCondition(DOFMatrix* matrix,
					      ElInfo* elInfo,
					      const DegreeOfFreedom* dofIndices,
					      const BoundaryType* localBound,
					      int nBasFcts)
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
