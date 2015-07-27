/** \file RobinBC.h */

#pragma once

#include "BoundaryCondition.h"
#include "expressions/expr_traits.hpp"

namespace AMDiS 
{
  class DOFMatrix;
  
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
    RobinBC(BoundaryType type, JExpr&& j_expr, AlphaExpr&& a_expr,
      	    const FiniteElemSpace *rowFeSpace,
      	    const FiniteElemSpace *colFeSpace = NULL)
      : BoundaryCondition(type, rowFeSpace, colFeSpace), 
      	neumannOperators(NULL), 
      	robinOperators(NULL)
    {
      using JExprType = traits::to_expr<JExpr>;
      using AlphaExprType = traits::to_expr<AlphaExpr>;
      
      if (traits::is_valid_arg<AlphaExpr>::value)
        init(JExprType::get(j_expr), AlphaExprType::get(a_expr));
      else
        init(JExprType::get(j_expr));
    }
    
    // TODO: move to private section:
    
    /// Implements BoundaryCondition::fillBoundaryCondition();
    virtual void fillBoundaryCondition(DOFMatrix* matrix,
                        				       ElInfo* elInfo,
                        				       const DegreeOfFreedom* dofIndices,
                        				       const BoundaryType* localBound,
                        				       int nBasFcts) override;
  
    /// Implements BoundaryCondition::fillBoundaryCondition();
    virtual void fillBoundaryCondition(DOFVectorBase<double>* vector, 
                        				       ElInfo* elInfo,
                        				       const DegreeOfFreedom* dofIndices,
                        				       const BoundaryType* localBound,
                        				       int nBasFcts) override;

    /// Implements BoundaryCondition::boundResidual();
    virtual double boundResidual(ElInfo *elInfo, 
                        				 DOFMatrix *matrix,
                        				 const DOFVectorBase<double> *dv) override;

  protected:
    template <class JExpr, class AlphaExpr>
    void init(JExpr&& j_expr, AlphaExpr&& a_expr);
    
    template <class JExpr>
    void init(JExpr&& j_expr);
    
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
      	      JExpr&& j_expr,
      	      const FiniteElemSpace *rowFeSpace,
      	      const FiniteElemSpace *colFeSpace = NULL)
      : RobinBC(type, j_expr, tag::dummy(), rowFeSpace, colFeSpace)
    {}
  };

} // end namespace AMDiS
