/** \file ZeroOrderTerm.h */

#pragma once

#include "AMDiS_fwd.h"
#include "OperatorTerm.h"

namespace AMDiS {

  /**
   * \ingroup Assembler
   *
   * \brief
   * Describes zero order terms: \f$ cu(\vec{x}) \f$.
   */
  class ZeroOrderTerm : public OperatorTerm
  {
  public:
    /// Constructor.
    ZeroOrderTerm(int deg) : OperatorTerm(deg) {}

    /// Destructor.
    virtual ~ZeroOrderTerm() {}
    
    /// Evaluates \f$ c \f$
    void getC(const ElInfo *elInfo, int nPoints, ElementVector& C) const
    {
      getCImpl(elInfo, nPoints, C);
    }

  private:
    virtual void getCImpl(const ElInfo *elInfo, int nPoints, ElementVector& C) const = 0;
  };
  
  
  /* ----- BASIC OPERATOR-TERMS USED IN EXPRESSIONS --------------------------- */
  
  
  template <class Expr>
  struct GenericZeroOrderTerm : public GenericOperatorTerm<Expr, 0>
  {
    GenericZeroOrderTerm(const Expr& expr_)
      : GenericOperatorTerm<Expr, 0>(expr_)
    { }

  private:
    /// Implemetation of ZeroOrderTerm::getC().
    virtual void getCImpl(const ElInfo *elInfo, int nPoints, ElementVector& C) const override
    {
      for (int iq = 0; iq < nPoints; iq++)
	C[iq] += this->expr(iq);
    }

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
			  const DenseVector<double>& uhAtQP,
			  const DenseVector<WorldVector<double> >& grdUhAtQP,
			  const DenseVector<WorldMatrix<double> >& D2UhAtQP,
			  DenseVector<double>& result,
			  double fac) const override;
    {
      for (int iq = 0; iq < nPoints; iq++)
	result[iq] += fac * this->expr(iq) * uhAtQP[iq];
    }
  };
  
} // end namespace AMDiS
