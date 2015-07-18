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



/** \file LazyOperatorTerm.h */

#ifndef AMDIS_LAZY_OPERATOR_TERM_H
#define AMDIS_LAZY_OPERATOR_TERM_H

#include "AMDiS_fwd.h"
#include "Traits.h"

namespace AMDiS 
{    
  struct LazyOperatorTermBase
  {
    template<typename List>
    void insertFeSpaces(List& feSpaces) const {}
    
    int getDegree() const { return 0; }
      
    template<typename OT>
    void initElement(OT* ot, const ElInfo* elInfo,
		      SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL) {}

    template<typename OT>
    void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
	    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL) {}
  protected: 
    virtual void helper() {}
  };


  template<typename Term>
  struct LazyOperatorTerm1 : public LazyOperatorTermBase
  {
    Term term;
    LazyOperatorTerm1(const Term& term_) : term(term_) {}
    
    template<typename List>
    inline void insertFeSpaces(List& feSpaces)
    {
      term.insertFeSpaces(feSpaces);
    }

    template<typename OT>
    inline void initElement(OT* ot, const ElInfo* elInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      term.initElement(ot, elInfo, subAssembler, quad, basisFct);
    }

    template<typename OT>
    inline void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      term.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
    }
    
    inline double operator()(const int& iq) const;
  };


  template<typename Term1, typename Term2>
  struct LazyOperatorTerm2 : public LazyOperatorTermBase
  {
    Term1 term1;
    Term2 term2;
    LazyOperatorTerm2(const Term1& term1_, const Term2& term2_) : term1(term1_), term2(term2_) {}
    
    template<typename List>
    inline void insertFeSpaces(List& feSpaces)
    {
      term1.insertFeSpaces(feSpaces);
      term2.insertFeSpaces(feSpaces);
    }

    template<typename OT>
    inline void initElement(OT* ot, const ElInfo* elInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      term1.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term2.initElement(ot, elInfo, subAssembler, quad, basisFct);
    }

    template<typename OT>
    inline void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      term1.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
      term2.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
    }
    
    double operator()(const int& iq) const;
  };


  template<typename Term1, typename Term2, typename Term3>
  struct LazyOperatorTerm3 : public LazyOperatorTermBase
  {
    Term1 term1;
    Term2 term2;
    Term3 term3;
    LazyOperatorTerm3(const Term1& term1_, const Term2& term2_, const Term3& term3_)
    : term1(term1_), term2(term2_), term3(term3_) {}
    
    template<typename List>
    inline void insertFeSpaces(List& feSpaces)
    {
      term1.insertFeSpaces(feSpaces);
      term2.insertFeSpaces(feSpaces);
      term3.insertFeSpaces(feSpaces);
    }

    template<typename OT>
    inline void initElement(OT* ot, const ElInfo* elInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      term1.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term2.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term3.initElement(ot, elInfo, subAssembler, quad, basisFct);
    }

    template<typename OT>
    inline void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      term1.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
      term2.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
      term3.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
    }
    
    inline double operator()(const int& iq) const;
  };

  /* 150203 added by Michael */
  template<typename Term1, typename Term2, typename Term3, typename Term4>
  struct LazyOperatorTerm4 : public LazyOperatorTermBase
  {
    Term1 term1;
    Term2 term2;
    Term3 term3;
    Term4 term4;
    LazyOperatorTerm4(const Term1& term1_, const Term2& term2_, const Term3& term3_, const Term4& term4_)
    : term1(term1_), term2(term2_), term3(term3_), term4(term4_) {}
    
    template<typename List>
    inline void insertFeSpaces(List& feSpaces)
    {
      term1.insertFeSpaces(feSpaces);
      term2.insertFeSpaces(feSpaces);
      term3.insertFeSpaces(feSpaces);
      term4.insertFeSpaces(feSpaces);
    }

    template<typename OT>
    inline void initElement(OT* ot, const ElInfo* elInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      term1.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term2.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term3.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term4.initElement(ot, elInfo, subAssembler, quad, basisFct);
    }

    template<typename OT>
    inline void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      term1.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
      term2.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
      term3.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
      term4.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
    }
    
    inline double operator()(const int& iq) const;
  };
  
  /* 150203 added by Michael */
  template<typename Term1, typename Term2, typename Term3, typename Term4, typename Term5>
  struct LazyOperatorTerm5 : public LazyOperatorTermBase
  {
    Term1 term1;
    Term2 term2;
    Term3 term3;
    Term4 term4;
    Term5 term5;
    LazyOperatorTerm5(const Term1& term1_, const Term2& term2_, const Term3& term3_, const Term4& term4_, const Term5& term5_)
    : term1(term1_), term2(term2_), term3(term3_), term4(term4_), term5(term5_) {}
    
    template<typename List>
    inline void insertFeSpaces(List& feSpaces)
    {
      term1.insertFeSpaces(feSpaces);
      term2.insertFeSpaces(feSpaces);
      term3.insertFeSpaces(feSpaces);
      term4.insertFeSpaces(feSpaces);
      term5.insertFeSpaces(feSpaces);
    }

    template<typename OT>
    inline void initElement(OT* ot, const ElInfo* elInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      term1.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term2.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term3.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term4.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term5.initElement(ot, elInfo, subAssembler, quad, basisFct);
    }

    template<typename OT>
    inline void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      term1.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
      term2.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
      term3.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
      term4.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
      term5.initElement(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct);
    }
    
    inline double operator()(const int& iq) const;
  };

} // end namespace AMDiS

#endif // AMDIS_LAZY_OPERATOR_TERM_H
