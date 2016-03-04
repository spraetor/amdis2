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



/** \file PenaltyOperator.h */

#ifndef AMDIS_PENALTYOPERATOR_H
#define AMDIS_PENALTYOPERATOR_H

#include "ElementLevelSet.h"
#include "Flag.h"
#include "Operator.h"
#include "ElementLevelSet.h"

namespace AMDiS
{
  class ElInfo;
  class FiniteElemSpace;
  class SurfaceOperator;
}

namespace compositeFEM
{

  using namespace AMDiS;
  using namespace std;

  class PenaltyOperator : public Operator
  {
  public:
    /// Constructor.
    PenaltyOperator(ElementLevelSet* elLS_,
                    double factor_,
                    bool penaltyCoeffFlag_,
                    FiniteElemSpace* rowFeSpace_,
                    FiniteElemSpace* colFeSpace_ = NULL)
      : Operator(rowFeSpace_, colFeSpace_),
        elLS(elLS_),
        elStatus(ElementLevelSet::LEVEL_SET_UNDEFINED),
        factor(factor_),
        penaltyCoeffFlag(penaltyCoeffFlag_),
        surfaceOp(NULL),
        dim(getRowFeSpace()->getMesh()->getDim()),
        degree(getRowFeSpace()->getBasisFcts()->getDegree())
    {
      TEST_EXIT(elLS->getLevelSetFct() && elLS->getMesh())
      ("ElementLevelSet not initialized!\n");

      tempCoords = new VectorOfFixVecs<DimVec<double>>(dim, dim, NO_INIT);
    }

    /// Destructor.
    ~PenaltyOperator();

    /**
     * Calculates the element matrix for this ElInfo and adds it multiplied by
     * factor to userMat.
     * In addition to \ref Operator::getElementMatrix(), it distinguishes
     * between elements divided by the boundary (level set zero intersects the
     * element) and elements lying completely inside or outside the integration
     * domain.
     */
    void getElementMatrix(const ElInfo* elInfo,
                          ElementMatrix& userMat,
                          double factor = 1.0);

    /**
     * Calculates the element vector for this ElInfo and adds it multiplied by
     * factor to userVec.
     * Similar to \ref getElementMatrix it distinguishes between elements
     * divided by the boundary and elements lying completely inside or outside
     * the integration domain.
     */
    void getElementVector(const ElInfo* elInfo,
                          DenseVector<double>& userVec,
                          double factor = 1.0);

  protected:
    /// Calculate coefficient of the penalty term.
    double getPenaltyCoeff(const ElInfo* elInfo);

  protected:
    /**
     * Holds level set function and functionalities for intersection point
     * calculation.
     */
    ElementLevelSet* elLS;

    /**
     * Indicator for the position of an element. It shows whether an element
     * lies inside the integration domain, outside of the integration domain
     * or on the boundary.
     */
    int elStatus;

    /**
     * Coefficient of the penalty term is 1/eps. eps depends on the size of
     * element and on factor.
     */
    double eps;

    /// Factor needed for eps.
    double factor;

    /**
     * Flag to indicate whether to use penalty coefficient.
     *   true - use penalty coefficient
     *   false - do not use
     */
    bool penaltyCoeffFlag;

    /// Surface operator for surface integration
    SurfaceOperator* surfaceOp;

    /// Problem Dimension.
    int dim;

    /// Degree of basis functions.
    int degree;

    /// Variables used for calculation.
    VectorOfFixVecs<DimVec<double>>* tempCoords;
  };

}

using compositeFEM::PenaltyOperator;

#endif  // AMDIS_PENALTYOPERATOR_H
