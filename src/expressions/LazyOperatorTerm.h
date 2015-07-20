/** \file LazyOperatorTerm.h */

#pragma once

#include "AMDiS_fwd.h"
#include "Traits.h"

namespace AMDiS 
{
  class SubAssembler;
  class Quadrature;
  class BasisFunction;
  class ElInfo;
  
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
    
    LazyOperatorTerm2(const Term1& term1_, const Term2& term2_) 
      : term1(term1_), term2(term2_) {}
    
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

  

  namespace detail 
  {
    /// Functor that initializes the feSpace list
    template <class List>
    struct InsertFeSpaces 
    {
      List& feSpaces;
      InsertFeSpaces(List& feSpaces_) : feSpaces(feSpaces_) {};
      
      template <class Term>
      void operator()(Term& term) {
	       term.insertFeSpaces(feSpaces);
      }
    };
    
    /// Functor that is called on each term to initialize it on an element
    template <class OT>
    struct InitElement 
    {
      OT* ot;
      const ElInfo *elInfo, *elInfo2;
      SubAssembler* subAssembler;
      Quadrature *quad;
      const BasisFunction *basisFct;
      
      InitElement(OT* ot_, const ElInfo* elInfo_, SubAssembler* subAssembler_, Quadrature *quad_, const BasisFunction *basisFct_)
      	: ot(ot_), elInfo(elInfo_), elInfo2(NULL), subAssembler(subAssembler_),
      	  quad(quad_), basisFct(basisFct_) {}
      
      InitElement(OT* ot_, const ElInfo* smallElInfo_, const ElInfo* largeElInfo_, SubAssembler* subAssembler_, Quadrature *quad_,  const BasisFunction *basisFct_)
      	: ot(ot_), elInfo(smallElInfo_), elInfo2(largeElInfo_), subAssembler(subAssembler_),
      	  quad(quad_), basisFct(basisFct_) {}
	  
      template<typename Term>
      void operator()(Term& term) {
      	if (elInfo2)
      	  term.initElement(ot, elInfo, elInfo2, subAssembler, quad, basisFct);
      	else
      	  term.initElement(ot, elInfo, subAssembler, quad, basisFct);
      }
    };
    
  } // end namespace detail
  

  /// Operator term with arbitrary number of sub-term (expressions)
  template <class... Terms>
  struct LazyOperatorTerms : public LazyOperatorTermBase
  {
    std::tuple<Terms...> term_tuple;
    
    template <class... Terms_>
    LazyOperatorTerms(Terms_&&... terms_)
      : term_tuple(std::forward<Terms_>(terms_)...) {}
    
    template <class List>
    void insertFeSpaces(List& feSpaces)
    {
      for_each(term_tuple, detail::InsertFeSpaces<List>(feSpaces));
    }

    template <class OT>
    void initElement(OT* ot, const ElInfo* elInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      for_each(term_tuple, detail::InitElement<OT>(ot, elInfo, subAssembler, quad, basisFct));
    }

    template <class OT>
    void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      for_each(term_tuple, detail::InitElement<OT>(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct));
    }
    
    double operator()(const int& iq) const;
  };
  
} // end namespace AMDiS
