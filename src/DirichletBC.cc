/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/


#include "DirichletBC.h"
#include "ElInfo.h"
#include "BasisFunction.h"
#include "DOFVector.h"
#include "DOFMatrix.h"

namespace AMDiS {
  
  namespace detail {
    
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
  }
}
