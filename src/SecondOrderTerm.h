/** \file SecondOrderTerm.h */

#pragma once

#include <vector>

#include "AMDiS_fwd.h"
#include "OperatorTerm.h"
#include "traits/traits.hpp"

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
    void getLALt(ElInfo const* elInfo,
                 std::vector<mtl::dense2D<double>>& result) const
    {
      getLALtImpl(elInfo, result);
    }

    /// Evaluation of \f$ A \nabla u(\vec{x}) \f$ at all quadrature points.
    void weakEval(std::vector<WorldVector<double>> const& grdUhAtQP,
                  std::vector<WorldVector<double>>& result) const
    {
      weakEvalImpl(grdUhAtQP, result);
    }

  private:
    virtual void getLALtImpl(ElInfo const*,
                             std::vector<mtl::dense2D<double>>&) const = 0;

    virtual void weakEvalImpl(std::vector<WorldVector<double>> const&,
                              std::vector<WorldVector<double>>&) const = 0;

  protected:
    /// Evaluation of \f$ \Lambda \cdot A \cdot \Lambda^t\f$.
    void lalt(DimVec<WorldVector<double>> const& Lambda, /* in */
              WorldMatrix<double> const& matrix,         /* in */
              mtl::dense2D<double>& LALt,                /* out */
              bool symm,
              double factor) const;


    /// Evaluation of \f$ \Lambda \cdot A \cdot \Lambda^t\f$ for \f$ A \f$
    /// the matrix having a ONE in the position \f$ (K,L) \f$
    /// and ZEROS in all other positions.
    void lalt_kl(DimVec<WorldVector<double>> const& Lambda, /* in */
                 int k, int l,
                 mtl::dense2D<double>& LALt,                /* out */
                 double factor) const;


    /// Evaluation of \f$ \Lambda \cdot A \cdot \Lambda^t\f$ for A equal to
    /// the identity.
    void l1lt(DimVec<WorldVector<double>> const& Lambda,   /* in */
              mtl::dense2D<double>& LALt,                  /* out */
              double factor) const;

  };


  /* ----- BASIC OPERATOR-TERMS USED IN EXPRESONS --------------------------- */


  /// SecondOrder OperatorTerm for expressions: < expr() * grad(u), grad(v) >
  template <class Term>
  class GenericSecondOrderTerm_1 : public GenericOperatorTerm<Term, 2>
  {
    using Super = GenericOperatorTerm<Term, 2>;

  public:
    /// Constructor.
    template <class Term_, 
	      class = Requires_t<concepts::Compatible<Term, Term_>>>
    GenericSecondOrderTerm_1(Term_&& term_)
      : Super(std::forward<Term_>(term_))
    {
      this->setSymmetric(true);
    }

  private:
    /// Implements SecondOrderTerm::getLALt().
    virtual void getLALtImpl(ElInfo const*,
                             std::vector<mtl::dense2D<double>>&) const override;

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
                          DenseVector<double> const& uhAtQP,
                          DenseVector<WorldVector<double>> const& grdUhAtQP,
                          DenseVector<WorldMatrix<double>> const& D2UhAtQP,
                          DenseVector<double>& result,
                          double f) const override;

    /// Implemetation of SecondOrderTerm::weakEval().
    virtual void weakEvalImpl(std::vector<WorldVector<double>> const&,
                              std::vector<WorldVector<double>>&) const override;
  };


  /// SecondOrder OperatorTerm for expressions: < ExprMat() * grad(u), grad(v) >
  template <class Term, bool symmetric = false>
  class GenericSecondOrderTerm_A : public GenericOperatorTerm<Term, 2>
  {
    using Super = GenericOperatorTerm<Term, 2>;

  public:
    template <class Term_, 
	      class = Requires_t<concepts::Compatible<Term, Term_>>>
    GenericSecondOrderTerm_A(Term_&& term_)
      : Super(std::forward<Term_>(term_))
    {
      this->setSymmetric(symmetric);
    }

  private:
    /// Implements SecondOrderTerm::getLALt().
    virtual void getLALtImpl(ElInfo const*,
                             std::vector<mtl::dense2D<double>>&) const override;

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
                          DenseVector<double> const& uhAtQP,
                          DenseVector<WorldVector<double>> const& grdUhAtQP,
                          DenseVector<WorldMatrix<double>> const& D2UhAtQP,
                          DenseVector<double>& result,
                          double factor) const override;

    /// Implemetation of SecondOrderTerm::weakEval().
    virtual void weakEvalImpl(std::vector<WorldVector<double>> const&,
                              std::vector<WorldVector<double>>&) const override;
  };



  /// SecondOrder OperatorTerm for expressions: < M * grad(u), grad(v) >, with M_ij = expr()
  template <int I, int J, class Term>
  class GenericSecondOrderTerm_ij : public GenericOperatorTerm<Term, 2>
  {
    using Super = GenericOperatorTerm<Term, 2>;

    int row, col;

  public:
    /// Constructor.
    template <class Term_, 
	      class = Requires_t<concepts::Compatible<Term, Term_>>>
    GenericSecondOrderTerm_ij(Term_&& term_)
      : Super(std::forward<Term_>(term_)), row(I), col(J)
    {
      this->setSymmetric(row == col);
    }

    /// Constructor.
    template <class Term_, 
	      class = Requires_t<concepts::Compatible<Term, Term_>>>
    GenericSecondOrderTerm_ij(Term_&& term_, int I0, int J0)
      : Super(std::forward<Term_>(term_)), row(I0), col(J0)
    {
      TEST_EXIT_DBG( I < 0 && I0 >= 0 && J < 0 && J0 >= 0 )
      ("You yould specify eather template<int I, int J>, or constructor(int I0, int J0)\n");
      this->setSymmetric(row == col);
    }

  private:
    /// Implements SecondOrderTerm::getLALt().
    virtual void getLALtImpl(ElInfo const*,
                             std::vector<mtl::dense2D<double>>&) const override;

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
                          DenseVector<double> const& uhAtQP,
                          DenseVector<WorldVector<double>> const& grdUhAtQP,
                          DenseVector<WorldMatrix<double>> const& D2UhAtQP,
                          DenseVector<double>& result,
                          double fac) const override;

    /// Implemetation of SecondOrderTerm::weakEval().
    virtual void weakEvalImpl(std::vector<WorldVector<double>> const&,
                              std::vector<WorldVector<double>>&) const override;
  };
  
} // end namspace AMDiS

#include "SecondOrderTerm.hh"
