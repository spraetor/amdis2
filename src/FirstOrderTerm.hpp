/** \file FirstOrderTerm.h */

#pragma once

#include <vector>

#include "AMDiS_fwd.hpp"
#include "OperatorTerm.hpp"

namespace AMDiS
{
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

    /// Evaluation of \f$ \Lambda b \f$.
    void getLb(ElInfo const* elInfo,
               std::vector<DenseVector<double>>& result) const
    {
      getLbImpl(elInfo, result);
    }

  private:
    virtual void getLbImpl(ElInfo const* elInfo,
                           std::vector<DenseVector<double>>& result) const = 0;

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
                          DenseVector<double> const& uhAtQP,
                          DenseVector<WorldVector<double>> const& grdUhAtQP,
                          DenseVector<WorldMatrix<double>> const& D2UhAtQP,
                          DenseVector<double>& result,
                          double factor) const override;

  protected:
    /// Evaluation of \f$ \Lambda \cdot b\f$ if b contains the value 1.0
    /// in each component.
    void l1(DimVec<WorldVector<double>> const& Lambda,
            DenseVector<double>& Lb,
            double factor) const;

    /// Evaluation of \f$ \Lambda \cdot b\f$.
    void lb(DimVec<WorldVector<double>> const& Lambda,
            WorldVector<double> const& b,
            DenseVector<double>& Lb,
            double factor) const;

    /// Evaluation of \f$ \Lambda \cdot b\f$.
    void lb_one(DimVec<WorldVector<double>> const& Lambda,
                DenseVector<double>& Lb,
                double factor) const;
  };


  /* ----- BASIC OPERATOR-TERMS USED IN EXPRESONS --------------------------- */


  /// FirstOrder OperatorTerm for expressions: < 1 * expr() * grad(u), v >
  template <class Term>
  class GenericFirstOrderTerm_1 : public GenericOperatorTerm<Term, 1>
  {
    using Super = GenericOperatorTerm<Term, 1>;

  public:
    /// Constructor.
    template <class Term_, 
	      class = Requires_t< concepts::Compatible<Term, Term_> >>
    GenericFirstOrderTerm_1(Term_&& term_)
      : Super(std::forward<Term_>(term_))
    { }

  private:
    /// Implements FirstOrderTerm::getLb().
    virtual void getLbImpl(ElInfo const* elInfo,
                           std::vector<DenseVector<double>>& Lb) const override;
  };


  /// FirstOrder OperatorTerm for expressions: < e_i * expr() * grad(u), v >
  template <int I, class Term>
  class GenericFirstOrderTerm_i : public GenericOperatorTerm<Term, 1>
  {
    using Super = GenericOperatorTerm<Term, 1>;

  public:
    /// Constructor.
    template <class Term_, 
	      class = Requires_t< concepts::Compatible<Term, Term_> >>
    GenericFirstOrderTerm_i(Term_&& term_)
      : Super(std::forward<Term_>(term_))
    {
      this->FirstOrderTerm::bOne = I;
    }

    template <class Term_, class = 
      Requires_t< concepts::Compatible<Term, Term_> >>
    GenericFirstOrderTerm_i(Term_&& term_, int I0)
      : Super(std::forward<Term_>(term_))
    {
      this->FirstOrderTerm::bOne = I0;
      TEST_EXIT_DBG( I < 0 && I0 >= 0 )
      ("You should specify eather template<int I>, or constructor(expr, int I0)\n");
    }

  private:
    /// Implements FirstOrderTerm::getLb().
    virtual void getLbImpl(ElInfo const* elInfo,
                           std::vector<DenseVector<double>>& Lb) const override;
  };


  /// FirstOrder OperatorTerm for expressions: < Term() * grad(u), v >
  template <class Term>
  class GenericFirstOrderTerm_b : public GenericOperatorTerm<Term, 1>
  {
    using Super = GenericOperatorTerm<Term, 1>;

  public:
    /// Constructor.
    template <class Term_, 
	      class = Requires_t< concepts::Compatible<Term, Term_> >>
    GenericFirstOrderTerm_b(Term_&& term_)
      : Super(std::forward<Term_>(term_))
    { }

  private:
    /// Implements FirstOrderTerm::getLb().
    virtual void getLbImpl(ElInfo const* elInfo,
                           std::vector<DenseVector<double>>& Lb) const override;
  };

} // end namespace AMDiS

#include "FirstOrderTerm.hh"
