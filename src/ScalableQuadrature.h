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



/** \file ScalableQuadrature.h */

#ifndef AMDIS_SCALABLEQUADRATURE_H
#define AMDIS_SCALABLEQUADRATURE_H

#include "Quadrature.h"
#include "FixVec.h"
#include "AMDiS_fwd.h"

namespace AMDiS {

  // ===============================================================================
  // Class Description:
  // The class \ref ScalableQuadrature holds the functionality for the manipulation 
  // of a quadrature formula used for the integration on subelements.
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
  // Main routines:
  // ScalableQuadrature() - Create a new quadrature formula which is a copy of the 
  //                        original quadrature formula and store the original
  //                        quadrature points in \ref oldLambda.
  // scaleQuadrature()    - Manipulate the original quadrature formula for a
  //                        subelement.
  // ===============================================================================
  class ScalableQuadrature : public Quadrature
  {
  public:
    /// Constructor
    ScalableQuadrature(Quadrature *quadrature);

    /// Destructor
    ~ScalableQuadrature() 
    {
      delete oldLambda;
    }
  
    /// Manipulates the quadrature points for the assemblage of a subelement.
    void scaleQuadrature(const SubElInfo& subElInfo);

    /** \brief
     * Scales the quadrature points using a matrix, that defines the transformation
     * from one element to another one.
     */
    void scaleQuadrature(DimMat<double> *scalMat);

    /// Get b-th coordinate of the a-th original quadrature point.
    double getOldLambda(int a, int b) const 
    {
      if (oldLambda) 
	return (*oldLambda)[a][b]; 
      else 
	return 0.0;
    }

    /// Get the a-th original quadrature point.
    DimVec<double>* getOldLambda(int a) const 
    {
      if (oldLambda) {
	return &((*oldLambda)[a]);
      } else {
	ERROR_EXIT("Should not happen!\n");
	return NULL;
      }
    }

    /** \brief
     * Returns, if the scaled quadrature point is within the target element (true),
     * or outside of the element (false).
     */
    bool isValid(int iq) const 
    {
      return valid[iq];
    }

  protected:
    /// Original quadrature points.
    VectorOfFixVecs<DimVec<double> > *oldLambda; 

    /** \brief
     * If quadrature points are scaled from on element to another one, they may
     * become invalide, i.e., there are outside of the target element. In this
     * case, the corresponding element of this field is false, otherwise true.
     */
    std::vector<bool> valid;
  };

}

#endif  // AMDIS_SCALABLEQUADRATURE_H
