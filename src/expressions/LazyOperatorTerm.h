/** \file LazyOperatorTerm.h */

#pragma once

#include "AMDiS_fwd.h"
#include "Traits.h"

namespace AMDiS 
{  
  struct LazyOperatorTermBase
  {
    template <class List>
    void insertFeSpaces(List& feSpaces) const {}
    
    int getDegree() const { return 0; }
      
    // TODO: remove OT-template and implement class with variadic templates
    template <class OT>
    void initElement(OT* ot, const ElInfo* elInfo,
		      SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL) {}
		    
  protected: 
    virtual void helper() {}
  };


  template <class Term>
  struct LazyOperatorTerm1 : public LazyOperatorTermBase
  {
    Term term;
    LazyOperatorTerm1(const Term& term_) : term(term_) {}
    
    template <class List>
    void insertFeSpaces(List& feSpaces)
    {
      term.insertFeSpaces(feSpaces);
    }

    template <class OT>
    void initElement(OT* ot, const ElInfo* elInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      term.initElement(ot, elInfo, subAssembler, quad, basisFct);
    }
    
    double operator()(const int& iq) const;
  };


  template <class Term1, class Term2>
  struct LazyOperatorTerm2 : public LazyOperatorTermBase
  {
    Term1 term1;
    Term2 term2;
    LazyOperatorTerm2(const Term1& term1_, const Term2& term2_) : term1(term1_), term2(term2_) {}
    
    template <class List>
    void insertFeSpaces(List& feSpaces)
    {
      term1.insertFeSpaces(feSpaces);
      term2.insertFeSpaces(feSpaces);
    }

    template <class OT>
    void initElement(OT* ot, const ElInfo* elInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      term1.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term2.initElement(ot, elInfo, subAssembler, quad, basisFct);
    }
    
    double operator()(const int& iq) const;
  };


  template <class Term1, class Term2, class Term3>
  struct LazyOperatorTerm3 : public LazyOperatorTermBase
  {
    Term1 term1;
    Term2 term2;
    Term3 term3;
    LazyOperatorTerm3(const Term1& term1_, const Term2& term2_, const Term3& term3_)
    : term1(term1_), term2(term2_), term3(term3_) {}
    
    template <class List>
    void insertFeSpaces(List& feSpaces)
    {
      term1.insertFeSpaces(feSpaces);
      term2.insertFeSpaces(feSpaces);
      term3.insertFeSpaces(feSpaces);
    }

    template <class OT>
    void initElement(OT* ot, const ElInfo* elInfo,
			    SubAssembler* subAssembler, Quadrature *quad, 
			    const BasisFunction *basisFct = NULL)
    {
      term1.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term2.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term3.initElement(ot, elInfo, subAssembler, quad, basisFct);
    }
    
    double operator()(const int& iq) const;
  };

  /* 150203 added by Michael */
  template <class Term1, class Term2, class Term3, class Term4>
  struct LazyOperatorTerm4 : public LazyOperatorTermBase
  {
    Term1 term1;
    Term2 term2;
    Term3 term3;
    Term4 term4;
    LazyOperatorTerm4(const Term1& term1_, const Term2& term2_, const Term3& term3_, const Term4& term4_)
    : term1(term1_), term2(term2_), term3(term3_), term4(term4_) {}
    
    template <class List>
    void insertFeSpaces(List& feSpaces)
    {
      term1.insertFeSpaces(feSpaces);
      term2.insertFeSpaces(feSpaces);
      term3.insertFeSpaces(feSpaces);
      term4.insertFeSpaces(feSpaces);
    }

    template <class OT>
    void initElement(OT* ot, const ElInfo* elInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      term1.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term2.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term3.initElement(ot, elInfo, subAssembler, quad, basisFct);
      term4.initElement(ot, elInfo, subAssembler, quad, basisFct);
    }
    
    double operator()(const int& iq) const;
  };
  
  /* 150203 added by Michael */
  template <typename Term1, typename Term2, typename Term3, typename Term4, typename Term5>
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
  
  // ===========================================================================
  
  // TODO: replace boost::fusion::for_each by own implementation
  template <class... Terms>
  struct LazyOperatorTerms : public LazyOperatorTermBase
  {
    std::tuple<Terms...> terms;
    
    template <class... Terms_>
    LazyOperatorTerm5(Terms_&&... terms_)
      : terms(std::fordward<Terms_>(terms)...)
    { }
    
    template <class List>
    void insertFeSpaces(List& feSpaces)
    {
      boost::fusion::for_each(terms, 
	[feSpaces&](auto& term) { term.insertFeSpaces(feSpaces); });
    }

    void initElement(const ElInfo* elInfo,
		     SubAssembler* subAssembler, Quadrature *quad, 
		     const BasisFunction *basisFct = NULL)
    {
      boost::fusion::for_each(terms, 
	[elInfo, subAssembler, quad, basisFct](auto& term) { 
	    term.initElement(elInfo, subAssembler, quad, basisFct); 
	});
    }
    
    //TODO: add initElement for quad, and basisFct separately
    
    double operator()(const int& iq) const;
  };

} // end namespace AMDiS
