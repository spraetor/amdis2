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


#include "ElInfo.h"

namespace AMDiS {

  template<typename T>
  void OperatorTerm::getVectorAtQPs(DOFVectorBase<T>* vec,
				    const ElInfo* elInfo, 
				    SubAssembler* subAssembler,
				    Quadrature *quad,
				    mtl::dense_vector<T>& vecAtQPs)
  {
    FUNCNAME_DBG("OperatorTerm::getVectorAtQPs()");

    TEST_EXIT_DBG(elInfo->getMesh() == vec->getFeSpace()->getMesh())
      ("There is something wrong!\n");
      
    subAssembler->getVectorAtQPs(vec, elInfo, quad, vecAtQPs);
  }


  template<typename T>
  void OperatorTerm::getVectorAtQPs(DOFVectorBase<T>* vec,
				    const ElInfo* smallElInfo, 
				    const ElInfo* largeElInfo, 
				    SubAssembler* subAssembler,
				    Quadrature *quad,
				    mtl::dense_vector<T>& vecAtQPs)
  {
    FUNCNAME("OperatorTerm::getVectorAtQPs()");

    TEST_EXIT(smallElInfo->getMesh() == vec->getFeSpace()->getMesh() ||
	      largeElInfo->getMesh() == vec->getFeSpace()->getMesh())
      ("There is something wrong!\n");

    if (smallElInfo->getLevel() == largeElInfo->getLevel()) {

      // Both elements have the same size, so we can use the simple procedure
      // to determine the vecAtQPs.
      
      if (vec->getFeSpace()->getMesh() == smallElInfo->getMesh())
	subAssembler->getVectorAtQPs(vec, smallElInfo, quad, vecAtQPs);
      else
	subAssembler->getVectorAtQPs(vec, largeElInfo, quad, vecAtQPs);

    } else {

      // The two elements are different. If the vector is defined on the mesh of the
      // small element, we can still use the simple procedure to determine the vecAtQPs.

      if (vec->getFeSpace()->getMesh() == largeElInfo->getMesh())
	subAssembler->getVectorAtQPs(vec, smallElInfo, largeElInfo, quad, vecAtQPs);
      else
	subAssembler->getVectorAtQPs(vec, smallElInfo, quad, vecAtQPs);
    }
  }

  /* ======== GradientAtQPs =========== */

  template<typename T>
  void OperatorTerm::getGradientsAtQPs( DOFVectorBase<T>* vec,
					const ElInfo* elInfo,
					SubAssembler* subAssembler,
					Quadrature *quad,
					mtl::dense_vector<typename GradientType<T>::type>& grdAtQPs)
  {
    FUNCNAME_DBG("OperatorTerm::getGradientsAtQPs()");

    TEST_EXIT_DBG(elInfo->getMesh() == vec->getFeSpace()->getMesh())
      ("There is something wrong!\n");

    subAssembler->getGradientsAtQPs(vec, elInfo, quad, grdAtQPs);
  }

  template<typename T>
  void OperatorTerm::getGradientsAtQPs( DOFVectorBase<T>* vec,
					const ElInfo* smallElInfo,
					const ElInfo* largeElInfo,
					SubAssembler* subAssembler,
					Quadrature *quad,
					mtl::dense_vector<typename GradientType<T>::type>& grdAtQPs)
  {
    FUNCNAME("OperatorTerm::getGradientsAtQPs()");

//     ERROR_EXIT("Not yet tested!\n");

    TEST_EXIT(smallElInfo->getMesh() == vec->getFeSpace()->getMesh() ||
	      largeElInfo->getMesh() == vec->getFeSpace()->getMesh())
      ("There is something wrong!\n");

    if (smallElInfo->getLevel() == largeElInfo->getLevel()) {

      // Both elements have the same size, so we can use the simple procedure
      // to determine the gradients.

      if (vec->getFeSpace()->getMesh() == smallElInfo->getMesh())
	subAssembler->getGradientsAtQPs(vec, smallElInfo, quad, grdAtQPs);
      else
	subAssembler->getGradientsAtQPs(vec, largeElInfo, quad, grdAtQPs);

    } else {

      // The two elements are different. If the vector is defined on the mesh of the
      // small element, we can still use the simple procedure to determine the gradients.

      if (vec->getFeSpace()->getMesh() == largeElInfo->getMesh())
	subAssembler->getGradientsAtQPs(vec, smallElInfo, largeElInfo, quad, grdAtQPs);
      else
	subAssembler->getGradientsAtQPs(vec, smallElInfo, quad, grdAtQPs);
    }
  }
}
