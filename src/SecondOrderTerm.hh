/** \file SecondOrderTerm.hh */

namespace AMDiS
{
  template <class Term>
  void GenericSecondOrderTerm_1<Term>::getLALtImpl(
    const ElInfo* elInfo,
    std::vector<mtl::dense2D<double>>& LALt) const
  {
    const DimVec<WorldVector<double>>& grdLambda = elInfo->getGrdLambda();
    const int nPoints = static_cast<int>(LALt.size());

    for (int iq = 0; iq < nPoints; iq++)
      this->l1lt(grdLambda, LALt[iq], this->term[iq]);
  }


  template <class Term>
  void GenericSecondOrderTerm_1<Term>::evalImpl(
    int nPoints,
    const DenseVector<double>& uhAtQP,
    const DenseVector<WorldVector<double>>& grdUhAtQP,
    const DenseVector<WorldMatrix<double>>& D2UhAtQP,
    DenseVector<double>& result,
    double f) const
  {
    int dow = Global::getGeo(WORLD);

    if (num_rows(D2UhAtQP) > 0)
    {
      for (int iq = 0; iq < nPoints; iq++)
      {
        double resultQP = 0.0;
        for (int i = 0; i < dow; i++)
          resultQP += D2UhAtQP[iq][i][i];

        result[iq] += resultQP * f * this->term[iq];
      }
    }
  }


  template <class Term>
  void GenericSecondOrderTerm_1<Term>::weakEvalImpl(
    const std::vector<WorldVector<double>>& grdUhAtQP,
    std::vector<WorldVector<double>>& result) const
  {
    int nPoints = grdUhAtQP.size();
    for (int iq = 0; iq < nPoints; iq++)
      result[iq] += this->term[iq] * grdUhAtQP[iq];
  }


  template <class Term, bool symmetric>
  void GenericSecondOrderTerm_A<Term, symmetric>::getLALtImpl(
    const ElInfo* elInfo,
    std::vector<mtl::dense2D<double>>& LALt) const
  {
    const DimVec<WorldVector<double>>& grdLambda = elInfo->getGrdLambda();
    const int nPoints = static_cast<int>(LALt.size());

    for (int iq = 0; iq < nPoints; iq++)
      this->lalt(grdLambda, this->term[iq], LALt[iq], symmetric, 1.0);
  }


  template <class Term, bool symmetric>
  void GenericSecondOrderTerm_A<Term, symmetric>::evalImpl(
    int nPoints,
    const DenseVector<double>& uhAtQP,
    const DenseVector<WorldVector<double>>& grdUhAtQP,
    const DenseVector<WorldMatrix<double>>& D2UhAtQP,
    DenseVector<double>& result,
    double factor) const
  {
    int dow = Global::getGeo(WORLD);

    for (int iq = 0; iq < nPoints; iq++)
    {
      double resultQP = 0.0;

      WorldMatrix<double> A = this->term[iq];

      if (num_rows(D2UhAtQP) > 0)
        for (int i = 0; i < dow; i++)
          for (int j = 0; j < dow; j++)
            resultQP += A[i][j] * D2UhAtQP[iq][j][i];
#if 0
      if (num_rows(grdUhAtQP) > 0)
        resultQP += (*divFct)(A) * grdUhAtQP[iq];
#endif
      result[iq] += resultQP * factor;
    }
  }


  template <class Term, bool symmetric>
  void GenericSecondOrderTerm_A<Term, symmetric>::weakEvalImpl(
    const std::vector<WorldVector<double>>& grdUhAtQP,
    std::vector<WorldVector<double>>& result) const
  {
    int nPoints = grdUhAtQP.size();
    WorldMatrix<double> A;
    for (int iq = 0; iq < nPoints; iq++)
      result[iq] += this->term[iq] * grdUhAtQP[iq];
  }


  template <int I, int J, class Term>
  void GenericSecondOrderTerm_ij<I, J, Term>::getLALtImpl(
    const ElInfo* elInfo,
    std::vector<mtl::dense2D<double>>& LALt) const
  {
    const DimVec<WorldVector<double>>& grdLambda = elInfo->getGrdLambda();
    const int nPoints = static_cast<int>(LALt.size());

    for (int iq = 0; iq < nPoints; iq++)
      this->lalt_kl(grdLambda, row, col, LALt[iq], this->term[iq]);
  }


  template <int I, int J, class Term>
  void GenericSecondOrderTerm_ij<I, J, Term>::evalImpl(
    int nPoints,
    const DenseVector<double>& uhAtQP,
    const DenseVector<WorldVector<double>>& grdUhAtQP,
    const DenseVector<WorldMatrix<double>>& D2UhAtQP,
    DenseVector<double>& result,
    double fac) const
  {
    if (num_rows(D2UhAtQP) > 0)
    {
      for (int iq = 0; iq < nPoints; iq++)
        result[iq] += D2UhAtQP[iq][row][col] * this->term[iq] * fac;
    }
  }


  template <int I, int J, class Term>
  void GenericSecondOrderTerm_ij<I, J, Term>::weakEvalImpl(
    const std::vector<WorldVector<double>>& grdUhAtQP,
    std::vector<WorldVector<double>>& result) const
  {
    int nPoints = grdUhAtQP.size();
    for (int iq = 0; iq < nPoints; iq++)
      result[iq][row] += grdUhAtQP[iq][col] * this->term[iq];
  }

} // end namspace AMDiS
