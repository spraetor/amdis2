/** \file LazyOperatorTerm.h */

#pragma once

#include "AMDiS_fwd.h"
#include "Traits.h"
#include <traits/concepts.hpp>

namespace AMDiS 
{
  class SubAssembler;
  // class Quadrature;
  // class BasisFunction;
  // class ElInfo;
  
  struct LazyOperatorTermBase
  {
    template <class List>
    void insertFeSpaces(List& feSpaces) const {}
    
    constexpr int getDegree() const { return 0; }
      
    void initElement(ElInfo const* elInfo,
                     SubAssembler* subAssembler, 
                     Quadrature* quad, 
                     BasisFunction const* basisFct = NULL) {}
  };


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
    struct InitElement 
    {
      const ElInfo *elInfo;
      SubAssembler* subAssembler;
      Quadrature *quad;
      const BasisFunction *basisFct;
      
      InitElement(const ElInfo* elInfo_, SubAssembler* subAssembler_, 
                  Quadrature *quad_, const BasisFunction *basisFct_)
      	: elInfo(elInfo_), subAssembler(subAssembler_),
      	  quad(quad_), basisFct(basisFct_) {}
	  
      template <class Term>
      void operator()(Term& term) {
    	  term.initElement(elInfo, subAssembler, quad, basisFct);
      }
    };
    
  } // end namespace detail
  

  /// Operator term with arbitrary number of sub-term (expressions)
  template <class... Terms>
  struct LazyOperatorTerms : public LazyOperatorTermBase
  {
    using Self = LazyOperatorTerms;
    std::tuple<Terms...> terms;
    
    template <class... Terms_, class = typename enable_if< concepts::Term<Terms_...> >::type>
    constexpr LazyOperatorTerms(Terms_&&... terms_)
      : terms(std::forward<Terms_>(terms_)...) {}
    
    template <class List>
    void insertFeSpaces(List& feSpaces)
    {
      for_each(terms, detail::InsertFeSpaces<List>(feSpaces));
    }

    void initElement(const ElInfo* elInfo,
            		     SubAssembler* subAssembler, Quadrature *quad, 
            		     const BasisFunction *basisFct = NULL)
    {
      for_each(terms, detail::InitElement(elInfo, subAssembler, quad, basisFct));
    }
    
    template <int N>
    /*constexpr*/int getDegree(int_<N>) const
    {
      return std::get<N>(terms).getDegree();
    } 
    
    template <int N>
    typename std::tuple_element<N, decltype(Self::terms)>::type &
    getTerm(int_<N>)
    {
      return std::get<N>(Self::terms);
    } 
    
    template <int N>
    typename std::tuple_element<N, decltype(Self::terms)>::type const&
    getTerm(int_<N>) const
    {
      return std::get<N>(Self::terms);
    }
  };  
  
} // end namespace AMDiS
