/** \file SecondOrderTerm.h */

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
   * Describes the second order terms: \f$ \nabla \cdot A \nabla u(\vec{x}) \f$
   */
  class SecondOrderTerm : public OperatorTerm
  {
  public:
    /// Constructor.
    SecondOrderTerm(int deg)
      : OperatorTerm(deg)
    {}

    /// Evaluation of \f$ \Lambda A \Lambda^t \f$ at all quadrature points.
    void getLALt(const ElInfo* elInfo,
                 std::vector<mtl::dense2D<double>>& result) const
    {
      getLALtImpl(elInfo, result);
    }

    /// Evaluation of \f$ A \nabla u(\vec{x}) \f$ at all quadrature points.
    void weakEval(const std::vector<WorldVector<double>>& grdUhAtQP,
                  std::vector<WorldVector<double>>& result) const
    {
      weakEvalImpl(grdUhAtQP, result);
    }

  private:
    virtual void getLALtImpl(const ElInfo* elInfo,
                             std::vector<mtl::dense2D<double>>& result) const = 0;

    virtual void weakEvalImpl(const std::vector<WorldVector<double>>& grdUhAtQP,
                              std::vector<WorldVector<double>>& result) const = 0;

  protected:
    /// Evaluation of \f$ \Lambda \cdot A \cdot \Lambda^t\f$.
    void lalt(const DimVec<WorldVector<double>>& Lambda,
              const WorldMatrix<double>& matrix,
              mtl::dense2D<double>& LALt,
              bool symm,
              double factor) const;


    /// Evaluation of \f$ \Lambda \cdot A \cdot \Lambda^t\f$ for \f$ A \f$
    /// the matrix having a ONE in the position \f$ (K,L) \f$
    /// and ZEROS in all other positions.
    void lalt_kl(const DimVec<WorldVector<double>>& Lambda,
                 int k, int l,
                 mtl::dense2D<double>& LALt,
                 double factor) const;


    /// Evaluation of \f$ \Lambda \cdot A \cdot \Lambda^t\f$ for A equal to
    /// the identity.
    void l1lt(const DimVec<WorldVector<double>>& Lambda,
              mtl::dense2D<double>& LALt,
              double factor) const;

  };


  /* ----- BASIC OPERATOR-TERMS USED IN EXPRESONS --------------------------- */


  /// SecondOrder OperatorTerm for expressions: < expr() * grad(u), grad(v) >
  template <class Term>
  struct GenericSecondOrderTerm_1 : public GenericOperatorTerm<Term, 2>
  {
    using Super = GenericOperatorTerm<Term, 2>;

    template <class Term_>
    GenericSecondOrderTerm_1(Term_&& term_)
      : Super(std::forward<Term_>(term_))
    {
      this->setSymmetric(true);
    }

  private:
    /// Implements SecondOrderTerm::getLALt().
    virtual void getLALtImpl(const ElInfo* elInfo,
                             std::vector<mtl::dense2D<double>>& LALt) const override;

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
                          const DenseVector<double>& uhAtQP,
                          const DenseVector<WorldVector<double>>& grdUhAtQP,
                          const DenseVector<WorldMatrix<double>>& D2UhAtQP,
                          DenseVector<double>& result,
                          double f) const override;

    /// Implemetation of SecondOrderTerm::weakEval().
    virtual void weakEvalImpl(const std::vector<WorldVector<double>>& grdUhAtQP,
                              std::vector<WorldVector<double>>& result) const override;
  };


  /// SecondOrder OperatorTerm for expressions: < ExprMat() * grad(u), grad(v) >
  template <class Term, bool symmetric = false>
  struct GenericSecondOrderTerm_A : public GenericOperatorTerm<Term, 2>
  {
    using Super = GenericOperatorTerm<Term, 2>;

    template <class Term_>
    GenericSecondOrderTerm_A(Term_&& term_)
      : Super(std::forward<Term_>(term_))
    {
      this->setSymmetric(symmetric);
    }

  private:
    /// Implements SecondOrderTerm::getLALt().
    virtual void getLALtImpl(const ElInfo* elInfo,
                             std::vector<mtl::dense2D<double>>& LALt) const override;

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
                          const DenseVector<double>& uhAtQP,
                          const DenseVector<WorldVector<double>>& grdUhAtQP,
                          const DenseVector<WorldMatrix<double>>& D2UhAtQP,
                          DenseVector<double>& result,
                          double factor) const override;

    /// Implemetation of SecondOrderTerm::weakEval().
    virtual void weakEvalImpl(const std::vector<WorldVector<double>>& grdUhAtQP,
                              std::vector<WorldVector<double>>& result) const override;
  };



  /// SecondOrder OperatorTerm for expressions: < M * grad(u), grad(v) >, with M_ij = expr()
  template <int I, int J, class Term>
  struct GenericSecondOrderTerm_ij : public GenericOperatorTerm<Term, 2>
  {
    using Super = GenericOperatorTerm<Term, 2>;

    int row, col;

    template <class Term_>
    GenericSecondOrderTerm_ij(Term_&& term_)
      : Super(std::forward<Term_>(term_)), row(I), col(J)
    {
      this->setSymmetric(row == col);
    }

    template <class Term_>
    GenericSecondOrderTerm_ij(Term_&& term_, int I0, int J0)
      : Super(std::forward<Term_>(term_)), row(I0), col(J0)
    {
      TEST_EXIT_DBG( I < 0 && I0 >= 0 && J < 0 && J0 >= 0 )
      ("You yould specify eather template<int I, int J>, or constructor(int I0, int J0)\n");
      this->setSymmetric(row == col);
    }

  private:
    /// Implements SecondOrderTerm::getLALt().
    virtual void getLALtImpl(const ElInfo* elInfo,
                             std::vector<mtl::dense2D<double>>& LALt) const override;

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
                          const DenseVector<double>& uhAtQP,
                          const DenseVector<WorldVector<double>>& grdUhAtQP,
                          const DenseVector<WorldMatrix<double>>& D2UhAtQP,
                          DenseVector<double>& result,
                          double fac) const override;

    /// Implemetation of SecondOrderTerm::weakEval().
    virtual void weakEvalImpl(const std::vector<WorldVector<double>>& grdUhAtQP,
                              std::vector<WorldVector<double>>& result) const override;
  };

  // ---------------------------------------------------------------------------
#if 0
  template <class Term>
  struct SOTWrapper
  {
    using OperatorTermType = GenericOperatorTerm<Term, 2>;
    SOTWrapper(OperatorTermType* ot_)
      : ot(ot_) {}

    OperatorTermType* getOperatorTerm()
    {
      return ot;
    }
    auto& getTerm() RETURNS( this->ot->term )

  private:
    OperatorTermType* ot;
  };


  template <class C, class C_ = Decay_t<C>>
  Requires_t<traits::is_scalar<Value_t<C>>, SOTWrapper<C_>>
                                        inline sot(C&& c)
  {
    return {new GenericSecondOrderTerm_1<C_>(std::forward<C>(c))};
  }


  template <class C, class C_ = Decay_t<C>>
  Requires_t<traits::is_scalar<Value_t<C>>, SOTWrapper<C_>>
                                        inline sot(C&& c, int i, int j)
  {
    return {new GenericSecondOrderTerm_ij<-1, -1, C_>(std::forward<C>(c), i,j)};
  }


  template <class C, class C_ = Decay_t<C>>
  Requires_t<traits::is_matrix<Value_t<C>>, SOTWrapper<C_>>
                                        inline sot(C&& c)
  {
    return {new GenericSecondOrderTerm_A<C_>(std::forward<C>(c))};
  }
#endif
} // end namspace AMDiS

#include "SecondOrderTerm.hh"
