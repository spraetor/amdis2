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



/** \file FirstOrderTerm.h */

#ifndef AMDIS_FIRST_ORDER_TERM_H
#define AMDIS_FIRST_ORDER_TERM_H

#include <vector>

#include "AMDiS_fwd.h"
#include "OperatorTerm.h"

namespace AMDiS {

  /**
   * \ingroup Assembler
   * 
   * \brief
   * Describes the first order terms: \f$ b \cdot \nabla u(\vec{x}) \f$
   */
  class FirstOrderTerm : public OperatorTerm
  {
  public:
    /// Constructor.
    FirstOrderTerm(int deg) 
      : OperatorTerm(deg) 
    {}

    /// Destructor.
    virtual ~FirstOrderTerm() {}
    
    /// Evaluation of \f$ \Lambda b \f$.
    void getLb(const ElInfo *elInfo,
	       std::vector<mtl::dense_vector<double> >& result) const
    {
      getLbImpl(elInfo, result);
    }
    
  private:
    virtual void getLbImpl(const ElInfo *elInfo,
			   std::vector<mtl::dense_vector<double> >& result) const = 0;

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
			  const mtl::dense_vector<double>& uhAtQP,
			  const mtl::dense_vector<WorldVector<double> >& grdUhAtQP,
			  const mtl::dense_vector<WorldMatrix<double> >& D2UhAtQP,
			  mtl::dense_vector<double>& result,
			  double factor) const override;

  protected:
    /// Evaluation of \f$ \Lambda \cdot b\f$ if b contains the value 1.0
    /// in each component.
    void l1(const DimVec<WorldVector<double> >& Lambda,
	    mtl::dense_vector<double>& Lb,
	    double factor) const;

    /// Evaluation of \f$ \Lambda \cdot b\f$.
    void lb(const DimVec<WorldVector<double> >& Lambda,
	    const WorldVector<double>& b,
	    mtl::dense_vector<double>& Lb,
	    double factor) const;

    /// Evaluation of \f$ \Lambda \cdot b\f$.
    void lb_one(const DimVec<WorldVector<double> >& Lambda,
		mtl::dense_vector<double>& Lb,
		double factor) const;
  };
  
  
  /* ----- BASIC OPERATOR-TERMS USED IN EXPRESONS --------------------------- */
  
  
  /// FirstOrder OperatorTerm for expressions: < 1 * expr() * grad(u), v >
  template <class Expr>
  struct GenericFirstOrderTerm_1 : public GenericOperatorTerm<Expr, 1>
  {
    GenericFirstOrderTerm_1(const Expr& expr_)
      : GenericOperatorTerm<Expr, 1>(expr_) 
    { }

  private:
    /// Implements FirstOrderTerm::getLb().
    virtual void getLbImpl(const ElInfo *elInfo,
			   std::vector<mtl::dense_vector<double> >& Lb) const override;
  };


  /// FirstOrder OperatorTerm for expressions: < e_i * expr() * grad(u), v >
  template <int I, class Expr>
  struct GenericFirstOrderTerm_i : public GenericOperatorTerm<Expr, 1>
  {
    GenericFirstOrderTerm_i(const Expr& expr_)
      : GenericOperatorTerm<Expr, 1>(expr_) 
    {
      this->FirstOrderTerm::bOne = I;
    }
    
    GenericFirstOrderTerm_i(const Expr& expr_, int I0)
      : GenericOperatorTerm<Expr, 1>(expr_)  
    {
      this->FirstOrderTerm::bOne = I0;
      TEST_EXIT_DBG( I < 0 && I0 >= 0 )
	("You should specify eather template<int I>, or constructor(expr, int I0)\n");
    }

  private:
    /// Implements FirstOrderTerm::getLb().
    virtual void getLbImpl(const ElInfo *elInfo, 
			   std::vector<mtl::dense_vector<double> >& Lb) const override;
  };


  /// FirstOrder OperatorTerm for expressions: < Expr() * grad(u), v >
  template <class Expr>
  struct GenericFirstOrderTerm_b : public GenericOperatorTerm<Expr, 1>
  {
    GenericFirstOrderTerm_b(const Expr& expr_)
      : GenericOperatorTerm<Expr, 1>(expr_)
    { }

  private:
    /// Implements FirstOrderTerm::getLb().
    virtual void getLbImpl(const ElInfo *elInfo, 
			   std::vector<mtl::dense_vector<double> >& Lb) const override;
  };
  
  
  /* ----- IMPLEMENTATION DETAILS ------------------------------------------- */
  
  
  template <class Expr>
  void GenericFirstOrderTerm_1<Expr>::getLbImpl(
	  const ElInfo *elInfo,
	  std::vector<mtl::dense_vector<double> >& Lb) const;
  {
    const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();
    const int nPoints = static_cast<int>(Lb.size());

    for (int iq = 0; iq < nPoints; iq++)
      this->l1(grdLambda, Lb[iq], this->term(iq));
  }
  

  template <int I, class Expr>
  void GenericFirstOrderTerm_i<I, Expr>::getLbImpl(
	  const ElInfo *elInfo,
	  std::vector<mtl::dense_vector<double> >& Lb) const
  {
    const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();
    const int nPoints = static_cast<int>(Lb.size());

    for (int iq = 0; iq < nPoints; iq++)
      this->lb_one(grdLambda, Lb[iq], this->term(iq));
  }
    
    
  template <class Expr>
  void GenericFirstOrderTerm_b<Expr>::getLbImpl(
	  const ElInfo *elInfo,
	  std::vector<mtl::dense_vector<double> >& Lb) const
  {
    const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();
    const int nPoints = static_cast<int>(Lb.size());

    for (int iq = 0; iq < nPoints; iq++)
      this->lb(grdLambda, this->term(iq), Lb[iq], 1.0);
  }
  
} // end namespace AMDiS

#endif
