/** \file FirstOrderTerm.h */

#pragma once

#include <vector>

#include "AMDiS_fwd.h"
#include "OperatorTerm.h"

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

    /// Destructor.
    virtual ~FirstOrderTerm() {}

    /// Evaluation of \f$ \Lambda b \f$.
    void getLb(const ElInfo* elInfo,
               std::vector<DenseVector<double>>& result) const
    {
      getLbImpl(elInfo, result);
    }

  private:
    virtual void getLbImpl(const ElInfo* elInfo,
                           std::vector<DenseVector<double>>& result) const = 0;

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
                          const DenseVector<double>& uhAtQP,
                          const DenseVector<WorldVector<double>>& grdUhAtQP,
                          const DenseVector<WorldMatrix<double>>& D2UhAtQP,
                          DenseVector<double>& result,
                          double factor) const override;

  protected:
    /// Evaluation of \f$ \Lambda \cdot b\f$ if b contains the value 1.0
    /// in each component.
    void l1(const DimVec<WorldVector<double>>& Lambda,
            DenseVector<double>& Lb,
            double factor) const;

    /// Evaluation of \f$ \Lambda \cdot b\f$.
    void lb(const DimVec<WorldVector<double>>& Lambda,
            const WorldVector<double>& b,
            DenseVector<double>& Lb,
            double factor) const;

    /// Evaluation of \f$ \Lambda \cdot b\f$.
    void lb_one(const DimVec<WorldVector<double>>& Lambda,
                DenseVector<double>& Lb,
                double factor) const;
  };


  /* ----- BASIC OPERATOR-TERMS USED IN EXPRESONS --------------------------- */


  /// FirstOrder OperatorTerm for expressions: < 1 * expr() * grad(u), v >
  template <class Term>
  struct GenericFirstOrderTerm_1 : public GenericOperatorTerm<Term, 1>
  {
    template <class Term_>
    GenericFirstOrderTerm_1(Term_&& term_)
      : GenericOperatorTerm<Term, 1>(std::forward<Term_>(term_))
    { }

  private:
    /// Implements FirstOrderTerm::getLb().
    virtual void getLbImpl(const ElInfo* elInfo,
                           std::vector<DenseVector<double>>& Lb) const override;
  };


  /// FirstOrder OperatorTerm for expressions: < e_i * expr() * grad(u), v >
  template <int I, class Term>
  struct GenericFirstOrderTerm_i : public GenericOperatorTerm<Term, 1>
  {
    template <class Term_>
    GenericFirstOrderTerm_i(Term_&& term_)
      : GenericOperatorTerm<Term, 1>(std::forward<Term_>(term_))
    {
      this->FirstOrderTerm::bOne = I;
    }

    template <class Term_>
    GenericFirstOrderTerm_i(Term_&& term_, int I0)
      : GenericOperatorTerm<Term, 1>(std::forward<Term_>(term_))
    {
      this->FirstOrderTerm::bOne = I0;
      TEST_EXIT_DBG( I < 0 && I0 >= 0 )
      ("You should specify eather template<int I>, or constructor(expr, int I0)\n");
    }

  private:
    /// Implements FirstOrderTerm::getLb().
    virtual void getLbImpl(const ElInfo* elInfo,
                           std::vector<DenseVector<double>>& Lb) const override;
  };


  /// FirstOrder OperatorTerm for expressions: < Term() * grad(u), v >
  template <class Term>
  struct GenericFirstOrderTerm_b : public GenericOperatorTerm<Term, 1>
  {
    template <class Term_>
    GenericFirstOrderTerm_b(Term_&& term_)
      : GenericOperatorTerm<Term, 1>(std::forward<Term_>(term_))
    { }

  private:
    /// Implements FirstOrderTerm::getLb().
    virtual void getLbImpl(const ElInfo* elInfo,
                           std::vector<DenseVector<double>>& Lb) const override;
  };

  // ---------------------------------------------------------------------------
#if 0
  template <class Term>
  struct FOTWrapper
  {
    using OperatorTermType = GenericOperatorTerm<Term, 1>;
    FOTWrapper(OperatorTermType* ot_, FirstOrderType type_)
      : ot(ot_),
        type(type_) {}

    OperatorTermType* getOperatorTerm()
    {
      return ot;
    }
    auto getTerm() RETURNS( this->ot->term )
    FirstOrderType getType()
    {
      return type;
    }

  private:
    OperatorTermType* ot;
    FirstOrderType type;
  };


  template <class C, class C_ = Decay_t<C>>
  Requires_t<traits::is_scalar<Value_t<C>>, FOTWrapper<C_>>
                                        inline fot(C&& c, FirstOrderType type = GRD_PHI)
  {
    return {new GenericFirstOrderTerm_1<C_>(std::forward<C>(c)), type};
  }


  template <class C, class C_ = Decay_t<C>>
  Requires_t<traits::is_scalar<Value_t<C>>, FOTWrapper<C_>>
                                        inline fot(C&& c, int i, FirstOrderType type = GRD_PHI)
  {
    return {new GenericFirstOrderTerm_i<-1, C_>(std::forward<C>(c), i), type};
  }


  template <class C, class C_ = Decay_t<C>>
  Requires_t<traits::is_vector<Value_t<C>>, FOTWrapper<C_>>
                                        inline fot(C&& c, FirstOrderType type = GRD_PHI)
  {
    return {new GenericFirstOrderTerm_b<C_>(std::forward<C>(c)), type};
  }
#endif
} // end namespace AMDiS

#include "FirstOrderTerm.hh"
