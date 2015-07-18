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


#include <vector>
#include <boost/numeric/mtl/mtl.hpp>
#include "SubElementAssembler.h"
#include "ScalableQuadrature.h"
#include "SubPolytope.h"

namespace compositeFEM {

  SubElementAssembler::SubElementAssembler(Operator *op, 
					   const FiniteElemSpace *rowFeSpace_,
					   const FiniteElemSpace *colFeSpace_)
    : StandardAssembler(op, NULL, NULL, NULL, NULL, rowFeSpace_, colFeSpace_)
  {
    /** 
     * Create a scalable quadrature for subassembler and replace the original 
     * quadrature of the subassembler with the scalable quadrature.
     */

    if (zeroOrderAssembler) {
      checkQuadratures();
      zeroOrderScalableQuadrature = 
	new ScalableQuadrature(zeroOrderAssembler->getQuadrature());
      zeroOrderAssembler->setQuadrature(zeroOrderScalableQuadrature);
    } else {
      zeroOrderScalableQuadrature = NULL;
    }

    if (firstOrderAssemblerGrdPsi) {
      checkQuadratures();
      firstOrderGrdPsiScalableQuadrature = 
	new ScalableQuadrature(firstOrderAssemblerGrdPsi->getQuadrature());
      firstOrderAssemblerGrdPsi->setQuadrature(firstOrderGrdPsiScalableQuadrature);
    } else {
      firstOrderGrdPsiScalableQuadrature = NULL;
    }

    if (firstOrderAssemblerGrdPhi) {
      checkQuadratures();
      firstOrderGrdPhiScalableQuadrature = 
	new ScalableQuadrature(firstOrderAssemblerGrdPhi->getQuadrature());
      firstOrderAssemblerGrdPhi->setQuadrature(firstOrderGrdPhiScalableQuadrature);
    } else {
      firstOrderGrdPhiScalableQuadrature = NULL;
    }

    if (secondOrderAssembler) {
      checkQuadratures();
      secondOrderScalableQuadrature = 
	new ScalableQuadrature(secondOrderAssembler->getQuadrature());
      secondOrderAssembler->setQuadrature(secondOrderScalableQuadrature);
    } else {
      secondOrderScalableQuadrature = NULL;
    }
  }

  void SubElementAssembler::scaleQuadratures(const SubElInfo& subElInfo)
  {
    if (zeroOrderAssembler) {
      zeroOrderScalableQuadrature->scaleQuadrature(subElInfo);
    }
    if (firstOrderAssemblerGrdPsi) {
      firstOrderGrdPsiScalableQuadrature->scaleQuadrature(subElInfo);
    }
    if (firstOrderAssemblerGrdPhi) {
      firstOrderGrdPhiScalableQuadrature->scaleQuadrature(subElInfo);
    }
    if (secondOrderAssembler) {
      secondOrderScalableQuadrature->scaleQuadrature(subElInfo);
    }
  }

  void SubElementAssembler::getSubElementVector(SubElInfo *subElInfo, 
						const ElInfo *elInfo, 
						ElementVector& userVec)
  {
    /** 
     * Manipulate the quadratures of the SubAssemblers for subelement.
     */
    scaleQuadratures(*subElInfo);

    calculateElementVector(elInfo, userVec);

    /** 
     * The integration has been performed with a quadrature living on element. The 
     * determinant of element has been used instead of the determinant of subelement. Thus
     * the result must be corrected with respect to subelement.
     */
    double corrFactor = subElInfo->getDet() / fabs(elInfo->getDet());
    for (int i = 0; i < nRow; i++)
      userVec[i] *= corrFactor;
  }

  void SubElementAssembler::getSubElementMatrix(SubElInfo *subElInfo, 
						const ElInfo *elInfo, 
						ElementMatrix& userMat)
  {
    /** 
     * Manipulate the quadratures of the SubAssemblers for subelement.
     */
    scaleQuadratures(*subElInfo);

    /** 
     * Integrate using the manipulated quadrature.
     */
    calculateElementMatrix(elInfo, userMat);

    /** 
     * The integration has been performed with a quadrature living on element. The 
     * determinant of element has been used instead of the determinant of subelement. 
     * Thus the result must be corrected with respect to subelement.
     */
    double corrFactor = subElInfo->getDet() / fabs(elInfo->getDet());
    for (int i = 0; i < nRow; i++) {
      for (int j = 0; j < nCol; j++) {
	userMat[i][j] *= corrFactor;
      }
    }
  }

  void SubElementAssembler::getSubPolytopeVector(SubPolytope *subPolytope,
						 SubElementAssembler *subElementAssembler,
						 const ElInfo *elInfo,
						 ElementVector& subPolVec)
  {
    /// Note: There is no reset of subPolVec.
    std::vector<SubElInfo *>::iterator it;
    ElementVector subElVec(nRow);

    /// Assemble for each subelement of subpolytope.
    for (it = subPolytope->getSubElementsBegin(); 
	 it != subPolytope->getSubElementsEnd();
	 it++) {
      set_to_zero(subElVec);
      subElementAssembler->getSubElementVector(*it, elInfo, subElVec);

      /// Add results for subelement to total result for subpolytope.
      subPolVec += subElVec;
    }
  }

  void SubElementAssembler::getSubPolytopeMatrix(SubPolytope *subPolytope,
						 SubElementAssembler *subElementAssembler,
						 const ElInfo *elInfo,
						 ElementMatrix& subPolMat)
  {
    /**
     * Note: There is no reset of subPolMat.
     */
    std::vector<SubElInfo *>::iterator it;
    ElementMatrix subElMat(nRow, nCol);

    /**
     * Assemble for each subelement of subpolytope.
     */
    for (it = subPolytope->getSubElementsBegin(); 
	 it != subPolytope->getSubElementsEnd();
	 it++) {
      set_to_zero(subElMat);
      subElementAssembler->getSubElementMatrix(*it, elInfo, subElMat);

      /**
       * Add results for subelement to total result for subpolytope.
       */
      for (int i = 0; i < nRow; i++)
	for (int j = 0; j < nCol; j++)
	  subPolMat[i][j] += subElMat[i][j];
    }
  }

}
