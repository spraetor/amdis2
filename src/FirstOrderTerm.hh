/** \file FirstOrderTerm.hh */

namespace AMDiS
{
  template <class Term>
  void GenericFirstOrderTerm_1<Term>::getLbImpl(ElInfo const* elInfo,
      std::vector<DenseVector<double>>& Lb) const
  {
    auto const& grdLambda = elInfo->getGrdLambda();
    const int nPoints = static_cast<int>(Lb.size());

    for (int iq = 0; iq < nPoints; iq++)
      this->l1(grdLambda, Lb[iq], this->term.evalAtIdx(iq));
  }


  template <int I, class Term>
  void GenericFirstOrderTerm_i<I, Term>::getLbImpl(ElInfo const* elInfo,
      std::vector<DenseVector<double>>& Lb) const
  {
    auto const& grdLambda = elInfo->getGrdLambda();
    const int nPoints = static_cast<int>(Lb.size());

    for (int iq = 0; iq < nPoints; iq++)
      this->lb_one(grdLambda, Lb[iq], this->term.evalAtIdx(iq));
  }


  template <class Term>
  void GenericFirstOrderTerm_b<Term>::getLbImpl(ElInfo const* elInfo,
      std::vector<DenseVector<double>>& Lb) const
  {
    auto const& grdLambda = elInfo->getGrdLambda();
    const int nPoints = static_cast<int>(Lb.size());

    for (int iq = 0; iq < nPoints; iq++)
      this->lb(grdLambda, this->term.evalAtIdx(iq), Lb[iq], 1.0);
  }

} // end namespace AMDiS
