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



/** \file CompositeFEMOperator.h */

#ifndef AMDIS_COMPOSITEFEMOPERATOR_H
#define AMDIS_COMPOSITEFEMOPERATOR_H

#include "FixVec.h"
#include "Flag.h"
#include "ElementLevelSet.h"
#include "Operator.h"

#include "ElementLevelSet.h"
#include "SubElementAssembler.h"

namespace AMDiS {
class ElInfo;
class FiniteElemSpace;
}

namespace compositeFEM {

using namespace AMDiS;
using namespace std;

// ===========================================================================
// === class CompositeFEMOperator ============================================
// ===========================================================================
//
// Class Description:
// Class CompositeFEMOperator realizes the calculation of element-vectors
// and element-matrices for subelements. 
// The boundary of the integration domain is the level set zero of the level
// set function (object elLevelSet). If the boundary intersects an 
// element the integration for subelements is used (see class SubPolytope). 
// Else, integration is done as usual.
// ============================================================================
class CompositeFEMOperator : public Operator
{
public:
  /// Constructor.
  CompositeFEMOperator(ElementLevelSet *elLS_,
		       const FiniteElemSpace *rowFeSpace_,
		       const FiniteElemSpace *colFeSpace_ = NULL)
    : Operator(rowFeSpace_, colFeSpace_),
      elLS(elLS_),
      subElementAssembler(NULL),
      elStatus(ElementLevelSet::LEVEL_SET_UNDEFINED)
  {}

  /// Destructor.
  ~CompositeFEMOperator()
  {
    if (subElementAssembler)
      delete subElementAssembler;
  }

  /** 
   * Calculates the element matrix for this ElInfo and adds it multiplied by
   * factor to userMat.
   * In addition to Operator::getElementMatrix(), it distinguishes
   * between elements divided by the boundary (level set zero intersects 
   * the element) and elements lying completely inside or outside the 
   * integration domain.
   */
  void getElementMatrix(const ElInfo *elInfo, 
			ElementMatrix& userMat, 
			double factor = 1.0);

  /** 
   * Calculates the element vector for this ElInfo and adds it multiplied by
   * factor to userVec.
   * Similar to getElementMatrix it distinguishes between elements divided by 
   * the boundary and elements lying completely inside or outside the 
   * integration domain.
   */
  void getElementVector(const ElInfo *elInfo, 
			ElementVector& userVec, 
			double factor = 1.0);

protected:
  /**
   * Holds level set function and functionalities for intersection point
   * calculation.
   */
  ElementLevelSet *elLS;

  /** 
   * Calculates the element matrix and/or the element vector for subelements. 
   * It is created especially for this Operator, when getElementMatrix()
   * or getElementVector is called for the first time.
   */
  SubElementAssembler *subElementAssembler;

  /**
   * Indicator for the position of an element. It shows whether an element 
   * lies inside the integration domain, outside of the integration domain 
   * or on the boundary.
   */
  int elStatus;
};

}

using compositeFEM::CompositeFEMOperator;

#endif  // AMDIS_COMPOSITEFEMOPERATOR_H
