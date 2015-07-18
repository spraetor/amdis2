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



/** \file RobinBC.h */

#ifndef AMDIS_ROBINBC_H
#define AMDIS_ROBINBC_H

#include "BoundaryCondition.h"
#include "DOFMatrix.h"

namespace AMDiS {

  /** 
   * \ingroup Assembler
   *
   * \brief
   * Sub class of BoundaryCondition. Implements Robin and Neumann boundary conditions.
   * The flux in normal direction is given by \f$ A\nabla u\cdot \nu\phi =: j = j_0 - \alpha u \f$ where
   * \f$ j_0 \f$ and \f$ alpha \f$ are functions evaluated at world coordinates
   * and \f$ u \f$ is the problem solution.
   */
  class RobinBC : public BoundaryCondition
  {
  public:
    /// Constructor. \f$ j \f$ and \f$ alpha \f$ are given as AbstractFunction objects.
    template <class JExpr, class AlphaExpr>
    RobinBC(BoundaryType type,
	    JExpr const& j,
	    AlphaExpr const& alpha,
	    const FiniteElemSpace *rowFeSpace,
	    const FiniteElemSpace *colFeSpace = NULL)
      : BoundaryCondition(type, rowFeSpace, colFeSpace), 
	neumannOperators(NULL), 
	robinOperators(NULL)
    {
      using JExprType = typename taits::to_expr<JExpr>::to;
      using AlphaExprType = typename taits::to_expr<AlphaExpr>::to;
      
      if (boost::is_same<AlphaExpr, tag::dummy>::value)
	init(JExprType::get(j), tag::dummy());
      else
	init(JExprType::get(j), AlphaExprType::get(j));
    }
    
    // TODO: move to private section:
    
    /// Implements BoundaryCondition::fillBoundaryCondition();
    virtual void fillBoundaryCondition(DOFMatrix* matrix,
				       ElInfo* elInfo,
				       const DegreeOfFreedom* dofIndices,
				       const BoundaryType* localBound,
				       int nBasFcts);
  
    /// Implements BoundaryCondition::fillBoundaryCondition();
    virtual void fillBoundaryCondition(DOFVectorBase<double>* vector, 
				       ElInfo* elInfo,
				       const DegreeOfFreedom* dofIndices,
				       const BoundaryType* localBound,
				       int nBasFcts);

    /// Implements BoundaryCondition::boundResidual();
    virtual double boundResidual(ElInfo *elInfo, 
				 DOFMatrix *matrix,
				 const DOFVectorBase<double> *dv);

  protected:
    template <class JExpr, class AlphaExpr>
    void init(JExpr const& j, AlphaExpr const& alpha);
    
    template <class JExpr>
    void init(JExpr const& j, tag::dummy);
    
  protected:
    /// Surface operators for each element side for the Neumann part.
    DimVec<SurfaceOperator*>* neumannOperators;

    /// Surface operators for each element side for the Robin part.
    DimVec<SurfaceOperator*>* robinOperators;

    VectorOfFixVecs<DimVec<double> >**coords;
  };


  struct NeumannBC : public RobinBC
  {
    template <class JExpr>
    NeumannBC(BoundaryType type,
	      JExpr const& j,
	      const FiniteElemSpace *rowFeSpace,
	      const FiniteElemSpace *colFeSpace = NULL)
      : RobinBC(type, j, tag::dummy(), rowFeSpace, colFeSpace)
    { }
  };

}

#endif
