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




#ifndef VELOCITYEXTFROMVELOCITYFIELD_H
#define VELOCITYEXTFROMVELOCITYFIELD_H

#include "AdaptInfo.h"
#include "BasisFunction.h"
#include "DOFVector.h"
#include "ElInfo.h"
#include "NormEps.h"
#include "VelocityExt.h"

namespace reinit {

using namespace AMDiS;

/////////////////////////////////////////////////////////////////////////////
//  c l a s s   V e l o c i t y E x t F r o m V e l o c i t y F i e l d    //
/////////////////////////////////////////////////////////////////////////////
class VelocityExtFromVelocityField : public VelocityExt
{
public:

  VelocityExtFromVelocityField(int dim_)
    : VelocityExt(dim_),
      lSFct(NULL),
      elNormalVel(dim_),
      basFcts(NULL)
  {
    lSFctVal.change_dim(dim + 1);
    
    // ===== set epsilon for norm regularization =====
    NormEps::setEps();
  }

  ~VelocityExtFromVelocityField()
  {}

  /**
   * Set velocity field.
   */
  void setVelocityField(std::vector<DOFVector<double> *> &velField_,
			const DOFVector<double> *lSFct_,
			DOFVector<double> *velDOF_)
  {
    FUNCNAME("VelocityExtFromVelocityField::setVelocityField()");

    nVelDOFs = 1;

    velField = velField_;
    lSFct = lSFct_;
    velDOF.clear();
    velDOF.push_back(velDOF_);
    origVelDOF.clear();

    TEST_EXIT(lSFct)("level set function not defined !\n");
    TEST_EXIT((int)velField.size() == dim)("illegal velocity field !\n");
    TEST_EXIT(lSFct->getFeSpace() == velDOF_->getFeSpace())
      ("different feSpaces !\n");

    basFcts = lSFct->getFeSpace()->getBasisFcts();
  };

  /**
   * Adaption of VelocityExt::calcVelocityBoundary().
   * Used to initialize normal velocity at interface elements.
   */
  void calcVelocityBoundary(DegreeOfFreedom *locInd, const int indexV);

  /**
   * Sets elInfo.
   */
  inline void setElInfo(ElInfo *elInfo_) 
  {
    elInfo = elInfo_;
  }

 protected:
  /**
   * Velocity field at interface.
   */
  std::vector<DOFVector<double> *> velField;

  /**
   * Interface is zero level set of level set function lSFct.
   */
  const DOFVector<double> *lSFct;

  /**
   * Normal velocity on a single element. Used to set normal velocity
   * in velDOF on interface elements.
   */
  DimVec<double> elNormalVel;

  /// Values of level set function in vertices of element.
  ElementVector lSFctVal;

  /// Basis functions.
  const BasisFunction *basFcts;

  /// ElInfo used in calcVelocityBoundary().
  ElInfo *elInfo;
};

}

using reinit::VelocityExtFromVelocityField;

#endif  // VELOCITYEXTFROMVELOCITYFIELD_H
