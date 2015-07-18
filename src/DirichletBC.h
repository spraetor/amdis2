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



/** \file DirichletBC.h */

#ifndef AMDIS_DIRICHLETBC_H
#define AMDIS_DIRICHLETBC_H

#include "AMDiS_fwd.h"
#include "BoundaryCondition.h"
#include "FixVec.h"

namespace AMDiS 
{  
  namespace detail 
  {
    class DirichletBC : public BoundaryCondition
    {
    public:
		  
      /// Constructor.
      DirichletBC(BoundaryType type,
		  const FiniteElemSpace *rowFeSpace,
		  const FiniteElemSpace *colFeSpace,
		  bool apply)
	: BoundaryCondition(type, rowFeSpace, colFeSpace), 
	  applyBC(apply) 
      { }

      /// Implementation of BoundaryCondition::fillBoundaryCondition().
      virtual void fillBoundaryCondition(DOFMatrix* matrix,
				ElInfo* elInfo,
				const DegreeOfFreedom* dofIndices,
				const BoundaryType* localBound,
				int nBasFcts) override;

      ///
      void initVector(DOFVectorBase<double>* vec);

      /// Implementation of BoundaryCondition::boundResidual().
      double boundResidual(ElInfo*, DOFMatrix*, const DOFVectorBase<double>*) 
      { 
	return 0.0; 
      }

      /// Because this is a Dirichlet boundary condition, always return true.
      bool isDirichlet() const
      { 
	return true; 
      }

      /// Returns \ref applyBC.
      bool applyBoundaryCondition() const
      {
	return applyBC;
      }

    protected:

      /// Defines, if the boundary condition must be applied to the matrix. See
      /// comment of \ref BoundaryCondition::applyBoundaryCondition.
      bool applyBC;
    };
  
  } // end namespace detail

  
  /**
  * \ingroup Assembler
  *
  * \brief
  * Sub class of BoundaryCondition. Implements Dirichlet boundary conditions.
  * A DOFVectors is set to a given value at a Dirichlet dof and in a DOFMatrix
  * the row corresponding to a Dirichlet dof is replaced by a row containing
  * only a 1.0 in the diagonal.
  */
  template <class Expr>
  class DirichletBC : public detail::DirichletBC
  {
    using Super = detail::DirichletBC;
    using ToExpr = typename traits::to_expr<Expr>::to;
    using ExprType typename ToExpr::type ExprType;
    
  public:
    /// Constructor.
    DirichletBC(BoundaryType type,
		Expr fct_,
		const FiniteElemSpace *rowFeSpace,
		const FiniteElemSpace *colFeSpace = NULL,
		bool apply = true)
      : Super(type, rowFeSpace, colFeSpace, apply),
	fct(ToExpr::get(fct_)),
	term(fct)
    { }
    
  
    /// Implementation of BoundaryCondition::fillBoundaryCondition().
    virtual void fillBoundaryCondition(DOFVectorBase<double>* vector, 
					ElInfo* elInfo,
					const DegreeOfFreedom* dofIndices,
					const BoundaryType* localBound,
					int nBasFcts) override
    {
      const BasisFunction *basFcts = rowFeSpace->getBasisFcts();
      // initialize expression on ElInfo
      fct.initElement(&term, elInfo, NULL, NULL, basFcts);
      for (int i = 0; i < nBasFcts; i++)
	if (localBound[i] == boundaryType) {
	  typename ExprType::value_type value = fct(i);
	  vector->setDirichletDofValue(dofIndices[i], value);
	  (*vector)[dofIndices[i]] = value;
	}
    }
			      
  protected:
    ExprType fct;
    GenericOperatorTerm<ExprType> term;
  };
  
} // end namespace AMDiS

#endif
