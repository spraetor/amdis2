/** \file ZeroOrderTerm.h */

#pragma once

#include "AMDiS_fwd.h"
#include "OperatorTerm.h"

namespace AMDiS 
{
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
    
    /// Evaluates \f$ c \f$
    void getC(ElInfo const* elInfo, int nPoints, DenseVector<double>& C) const
    {
      getCImpl(elInfo, nPoints, C);
    }

  private:
    virtual void getCImpl(ElInfo const* elInfo, int nPoints, DenseVector<double>& C) const = 0;
  };
  
  
  /* ----- BASIC OPERATOR-TERMS USED IN EXPRESSIONS --------------------------- */
  
  
  template <class Term>
  struct GenericZeroOrderTerm : public GenericOperatorTerm<Term, 0>
  {
    template <class Term_>
    GenericZeroOrderTerm(Term_&& term_)
      : GenericOperatorTerm<Term, 0>(std::forward<Term_>(term_))
    {}

  private:
    /// Implemetation of ZeroOrderTerm::getC().
    virtual void getCImpl(ElInfo const* elInfo, int nPoints, DenseVector<double>& C) const override
    {
      for (int iq = 0; iq < nPoints; iq++)
        C[iq] += this->term(iq);
    }

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
                  			  DenseVector<double> const& uhAtQP,
                  			  DenseVector<WorldVector<double> > const& grdUhAtQP,
                  			  DenseVector<WorldMatrix<double> > const& D2UhAtQP,
                  			  DenseVector<double>& result,
                  			  double fac) const override
    {
      for (int iq = 0; iq < nPoints; iq++)
        result[iq] += fac * this->term(iq) * uhAtQP[iq];
    }
  };
  
} // end namespace AMDiS
