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



/** \file SecondOrderTerm.h */

#ifndef AMDIS_SECOND_ORDER_TERM_H
#define AMDIS_SECOND_ORDER_TERM_H

#include <vector>

#include "AMDiS_fwd.h"
#include "OperatorTerm.h"

namespace AMDiS {

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

    /// Destructor.
    virtual ~SecondOrderTerm() 
    {}

    /// Evaluation of \f$ \Lambda A \Lambda^t \f$ at all quadrature points.
    void getLALt(const ElInfo *elInfo, 
		 std::vector<mtl::dense2D<double> > &result) const
    {
      getLALtImpl(elInfo, result);
    }

    /// Evaluation of \f$ A \nabla u(\vec{x}) \f$ at all quadrature points.
    void weakEval(const std::vector<WorldVector<double> > &grdUhAtQP,
		  std::vector<WorldVector<double> > &result) const
    {
      weakEvalImpl(grdUhAtQP, result);
    }

  private:
    virtual void getLALtImpl(const ElInfo *elInfo, 
			     std::vector<mtl::dense2D<double> > &result) const = 0;

    virtual void weakEvalImpl(const std::vector<WorldVector<double> > &grdUhAtQP,
			      std::vector<WorldVector<double> > &result) const = 0;
			  
  protected:
    /// Evaluation of \f$ \Lambda \cdot A \cdot \Lambda^t\f$.
    void lalt(const DimVec<WorldVector<double> >& Lambda,
	      const WorldMatrix<double>& matrix,
	      mtl::dense2D<double>& LALt,
	      bool symm,
	      double factor) const;

	      
    /// Evaluation of \f$ \Lambda \cdot A \cdot \Lambda^t\f$ for \f$ A \f$
    /// the matrix having a ONE in the position \f$ (K,L) \f$
    /// and ZEROS in all other positions.
    void lalt_kl(const DimVec<WorldVector<double> >& Lambda,
		 int k, int l,
		 mtl::dense2D<double>& LALt,
		 double factor) const;

    
    /// Evaluation of \f$ \Lambda \cdot A \cdot \Lambda^t\f$ for A equal to 
    /// the identity.
    void l1lt(const DimVec<WorldVector<double> >& Lambda, 
		     mtl::dense2D<double>& LALt,
		     double factor) const;
  
  };
  
  
  /* ----- BASIC OPERATOR-TERMS USED IN EXPRESONS --------------------------- */
    
  
  /// SecondOrder OperatorTerm for expressions: < expr() * grad(u), grad(v) >
  template <class Term>
  struct GenericSecondOrderTerm_1 : public GenericOperatorTerm<Term, 2>
  {
    GenericSecondOrderTerm_1(const Term& term_)
      : GenericOperatorTerm<Term, 2>(term_) 
    {
      this->setSymmetric(true);
    }

  private:
    /// Implements SecondOrderTerm::getLALt().
    virtual void getLALtImpl(const ElInfo *elInfo, 
			     std::vector<mtl::dense2D<double> > &LALt) const override;

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
			  const mtl::dense_vector<double>& uhAtQP,
			  const mtl::dense_vector<WorldVector<double> >& grdUhAtQP,
			  const mtl::dense_vector<WorldMatrix<double> >& D2UhAtQP,
			  mtl::dense_vector<double>& result,
			  double f) override;

    /// Implemetation of SecondOrderTerm::weakEval().
    virtual void weakEvalImpl(const std::vector<WorldVector<double> > &grdUhAtQP,
			      std::vector<WorldVector<double> > &result) override;
  };


  /// SecondOrder OperatorTerm for expressions: < ExprMat() * grad(u), grad(v) >
  template <class Term, bool symmetric = false>
  struct GenericSecondOrderTerm_A : public GenericOperatorTerm<Term, 2>
  {  
    GenericSecondOrderTerm_A(const Term& term_)
      : GenericOperatorTerm<Term, 2>(term_) 
    {
      this->setSymmetric(symmetric);
    }
    
  private:
    /// Implements SecondOrderTerm::getLALt().
    virtual void getLALtImpl(const ElInfo *elInfo, 
			     std::vector<mtl::dense2D<double> > &LALt) const override;

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
			  const mtl::dense_vector<double>& uhAtQP,
			  const mtl::dense_vector<WorldVector<double> >& grdUhAtQP,
			  const mtl::dense_vector<WorldMatrix<double> >& D2UhAtQP,
			  mtl::dense_vector<double>& result,
			  double factor) override;

    /// Implemetation of SecondOrderTerm::weakEval().
    virtual void weakEvalImpl(const std::vector<WorldVector<double> > &grdUhAtQP,
			      std::vector<WorldVector<double> > &result) override;
  };



  /// SecondOrder OperatorTerm for expressions: < M * grad(u), grad(v) >, with M_ij = expr()
  template <int I, int J, class Term>
  struct GenericSecondOrderTerm_ij : public GenericOperatorTerm<Term, 2>
  {
    int row, col;
    
    GenericSecondOrderTerm_ij(const Term& term_)
      : GenericOperatorTerm<Term, 2>(term_), row(I), col(J)
    {
      this->setSymmetric(row == col);
    }
    
    GenericSecondOrderTerm_ij(const Term& term_, int I0, int J0)
      : GenericOperatorTerm<Term, 2>(term_), row(I0), col(J0)
    {
      this->setSymmetric(row == col);
      TEST_EXIT_DBG( I < 0 && I0 >= 0 && J < 0 && J0 >= 0 ) 
	("You yould specify eather template<int I, int J>, or constructor(int I0, int J0)\n");
    }

  private:
    /// Implements SecondOrderTerm::getLALt().
    virtual void getLALtImpl(const ElInfo *elInfo, 
			     std::vector<mtl::dense2D<double> > &LALt) const override;

    /// Implemetation of OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
			  const mtl::dense_vector<double>& uhAtQP,
			  const mtl::dense_vector<WorldVector<double> >& grdUhAtQP,
			  const mtl::dense_vector<WorldMatrix<double> >& D2UhAtQP,
			  mtl::dense_vector<double>& result,
			  double fac) override;

    /// Implemetation of SecondOrderTerm::weakEval().
    virtual void weakEvalImpl(const std::vector<WorldVector<double> > &grdUhAtQP,
			      std::vector<WorldVector<double> > &result) override;
  };
  
  
  /* ----- IMPLEMENTATION DETAILS ------------------------------------------- */
  
  
  template <class Term>
  void GenericSecondOrderTerm_1<Term>::getLALtImpl(
	  const ElInfo *elInfo, 
	  std::vector<mtl::dense2D<double> > &LALt) const
  {
    const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();
    const int nPoints = static_cast<int>(LALt.size());

    for (int iq = 0; iq < nPoints; iq++) 
      this->l1lt(grdLambda, LALt[iq], this->term(iq));
  }
  
  
  template <class Term>
  void GenericSecondOrderTerm_1<Term>::evalImpl(
	  int nPoints,
	  const mtl::dense_vector<double>& uhAtQP,
	  const mtl::dense_vector<WorldVector<double> >& grdUhAtQP,
	  const mtl::dense_vector<WorldMatrix<double> >& D2UhAtQP,
	  mtl::dense_vector<double>& result,
	  double f) 
  {
    int dow = Global::getGeo(WORLD);

    if (num_rows(D2UhAtQP) > 0) {
      for (int iq = 0; iq < nPoints; iq++) {
	double resultQP = 0.0;
	for (int i = 0; i < dow; i++) {
	  resultQP += D2UhAtQP[iq][i][i];
	}
	result[iq] += resultQP * f * this->term(iq);
      }
    }
  }

  
  template <class Term>
  void GenericSecondOrderTerm_1<Term>::weakEvalImpl(
	  const std::vector<WorldVector<double> > &grdUhAtQP,
	  std::vector<WorldVector<double> > &result) 
  {
    int nPoints = grdUhAtQP.size();
    for (int iq = 0; iq < nPoints; iq++)
      axpy(this->term(iq), grdUhAtQP[iq], result[iq]);
  }
  
  
  template <class Term, bool symmetric>
  void GenericSecondOrderTerm_A<Term, symmetric>::getLALtImpl(
	  const ElInfo *elInfo, 
	  std::vector<mtl::dense2D<double> > &LALt) const
  {
    const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();
    const int nPoints = static_cast<int>(LALt.size());

    for (int iq = 0; iq < nPoints; iq++)
      this->lalt(grdLambda, this->term(iq), LALt[iq], symmetric, 1.0);
  }

  
  template <class Term, bool symmetric>
  void GenericSecondOrderTerm_A<Term, symmetric>::evalImpl(
	  int nPoints,
	  const mtl::dense_vector<double>& uhAtQP,
	  const mtl::dense_vector<WorldVector<double> >& grdUhAtQP,
	  const mtl::dense_vector<WorldMatrix<double> >& D2UhAtQP,
	  mtl::dense_vector<double>& result,
	  double factor) 
  {
    int dow = Global::getGeo(WORLD);

    for (int iq = 0; iq < nPoints; iq++) {
      double resultQP = 0.0;

      WorldMatrix<double> A = this->term(iq);

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
	  const std::vector<WorldVector<double> > &grdUhAtQP,
	  std::vector<WorldVector<double> > &result)
  {
    int nPoints = grdUhAtQP.size();
    WorldMatrix<double> A;
    for (int iq = 0; iq < nPoints; iq++) {
      result[iq] += this->term(iq) * grdUhAtQP[iq];
    }    
  }
  
  
  template <int I, int J, class Term>
  void GenericSecondOrderTerm_ij<I, J Term>::getLALtImpl(
	  const ElInfo *elInfo, 
	  std::vector<mtl::dense2D<double> > &LALt) const
  {
    const DimVec<WorldVector<double> > &grdLambda = elInfo->getGrdLambda();
    const int nPoints = static_cast<int>(LALt.size());

    for (int iq = 0; iq < nPoints; iq++)
      this->lalt_kl(grdLambda, row, col, LALt[iq], this->term(iq));
  }  

  
  template <int I, int J, class Term>
  void GenericSecondOrderTerm_ij<I, J Term>::evalImpl(
	  int nPoints,
	  const mtl::dense_vector<double>& uhAtQP,
	  const mtl::dense_vector<WorldVector<double> >& grdUhAtQP,
	  const mtl::dense_vector<WorldMatrix<double> >& D2UhAtQP,
	  mtl::dense_vector<double>& result,
	  double fac) 
  {
    if (num_rows(D2UhAtQP) > 0) {
      for (int iq = 0; iq < nPoints; iq++)
	result[iq] += D2UhAtQP[iq][row][col] * this->term(iq) * fac;
    }
  }

  
  template <int I, int J, class Term>
  void GenericSecondOrderTerm_ij<I, J Term>::weakEvalImpl(
	  const std::vector<WorldVector<double> > &grdUhAtQP,
	  std::vector<WorldVector<double> > &result) 
  {
    int nPoints = grdUhAtQP.size();
    for (int iq = 0; iq < nPoints; iq++)
      result[iq][row] += grdUhAtQP[iq][col] * this->term(iq);
  }
  
} // end namspace AMDiS

#endif
