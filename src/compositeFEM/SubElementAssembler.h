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



/** \file SubElementAssembler.h */

#ifndef AMDIS_SUBELEMENTASSEMBLER_H
#define AMDIS_SUBELEMENTASSEMBLER_H

#include "Assembler.h"
#include "SubElInfo.h"
#include "ScalableQuadrature.h"


namespace compositeFEM {
  
  using namespace AMDiS;

  class SubPolytope;

  // ============================================================================
  // ===== class SubElementAssembler ============================================
  // ============================================================================
  //
  // Class Desription:
  // The class \ref SubElementAssembler holds the routines for the assemblage on
  // subpolytopes and subelements.
  // The integration on a subpolytope takes place by integrating on the 
  // subelements and afterwards summing up the results.
  // If S' is a sublement and S the element containing S', the numerical integration
  // on S' is done by integration on the element S with a manipulated quadrature 
  // formula and multiplication of the result with a correction term consisting of
  // the determinants corresponding to S and S'. That means, the quadrature points
  // are manipulated in the following way:
  // The original quadrature formula holds quadrature points which are given in
  // barycentric coordinates with respect to S'. Now we express these quadrature 
  // points in barycentric coordinates with respect to S and obtain the manipulated
  // quadrature points we need. Obviously, the corresponding quadrature points 
  // in world coordinates coincide. Thus the numerical integration on S with the
  // manipulated quadrature formula gives us the same result as the numerical
  // integration on S' with the original quadrature formula up to the determinant.
  // This method for integration on subelements allows the reuse of the routines for
  // the integration on elements.
  //
  // The manipulation of the quadrature formula takes place for the corresponding
  // assembler type (ZeroOrderAssembler, FirstOrderAssemblerGrdPsi, ...).
  //
  // Main routines:
  // SubElementAssembler()  - Creates a scalable quadrature for the appropriate
  //                          assembler type and assigns this quadrature to the 
  //                          assembler.
  // scaleQuadratures()     - Manipulates the scalable quadrature of the 
  //                          appropriate assembler type with respect to a 
  //                          subelement.
  // getSubPolytopeVector() - Calculates the righthandside vector for a polytope.
  // getSubPolytopeMatrix() - Calculates the system matrix for a polytope.
  // getSubElementVector()  - Calculates the righthandside vector for a subelement.
  // getSubElementMatrix()  - Calculates the system matrix for a subelement.
  // ============================================================================
  class SubElementAssembler : public StandardAssembler
  {
  public:
    SubElementAssembler(Operator *op, 
			const FiniteElemSpace *rowFeSpace,
			const FiniteElemSpace *colFeSpace = NULL);

    virtual ~SubElementAssembler()
    {
      if (zeroOrderScalableQuadrature)
	delete zeroOrderScalableQuadrature;
      if (firstOrderGrdPsiScalableQuadrature)
	delete firstOrderGrdPsiScalableQuadrature;
      if (firstOrderGrdPhiScalableQuadrature)
	delete firstOrderGrdPhiScalableQuadrature;
      if (secondOrderScalableQuadrature) 
	delete secondOrderScalableQuadrature;
    }

    void scaleQuadratures(const SubElInfo& subElInfo);

    void getSubElementVector(SubElInfo *subElInfo, 
			     const ElInfo *elInfo, 
			     DenseVector<double>& userVec);

    void getSubElementMatrix(SubElInfo *subElInfo, 
			     const ElInfo *elInfo, 
			     ElementMatrix& userMat);

    void getSubPolytopeVector(SubPolytope *subPolytope,
			      SubElementAssembler *subElementAssembler,
			      const ElInfo *elInfo,
			      DenseVector<double>& userVec);

    void getSubPolytopeMatrix(SubPolytope *subPolytope,
			      SubElementAssembler *subElementAssembler,
			      const ElInfo *elInfo,
			      ElementMatrix& userMat);

  protected:
    ScalableQuadrature *zeroOrderScalableQuadrature;
    ScalableQuadrature *firstOrderGrdPsiScalableQuadrature;
    ScalableQuadrature *firstOrderGrdPhiScalableQuadrature;
    ScalableQuadrature *secondOrderScalableQuadrature;
  };

}

#endif  // AMDIS_SUBELEMENTASSEMBLER_H
