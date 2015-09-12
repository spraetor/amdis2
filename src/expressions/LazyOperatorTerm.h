/** \file LazyOperatorTerm.h */

#pragma once

#include <AMDiS_fwd.h>
#include <Traits.h>
#include <traits/basic.hpp>
#include <utility/foreach.hpp>
#include "TermConcepts.hpp"

namespace AMDiS
{
  class SubAssembler;

  struct LazyOperatorTermBase
  {
    template <class List>
    void insertFeSpaces(List& feSpaces) const {}

    constexpr int getDegree() const
    {
      return 0;
    }

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
      InsertFeSpaces(List& feSpaces_)
        : feSpaces(feSpaces_)
      {}

      template <class Term>
      void operator()(Term& term)
      {
        term.insertFeSpaces(feSpaces);
      }

    private:
      List& feSpaces;
    };

    /// Functor that is called on each term to initialize it on an element
    struct InitElement
    {
      ElInfo const* elInfo;
      SubAssembler* subAssembler;
      Quadrature* quad;
      BasisFunction const* basisFct;

      InitElement(ElInfo const* elInfo_,
                  SubAssembler* subAssembler_,
                  Quadrature* quad_,
                  BasisFunction const* basisFct_)
        : elInfo(elInfo_),
          subAssembler(subAssembler_),
          quad(quad_),
          basisFct(basisFct_)
      {}

      template <class Term>
      void operator()(Term& term)
      {
        term.initElement(elInfo, subAssembler, quad, basisFct);
      }
    };

  } // end namespace detail


  /// Operator term with arbitrary number of sub-term (expressions)
  template <class... Terms>
  struct LazyOperatorTerms 
    : public LazyOperatorTermBase
  {
    using Self = LazyOperatorTerms;

    template <class... Terms_,
      class = Requires_t<concepts::Term<Terms_...>>>
    constexpr LazyOperatorTerms(Terms_&&... terms_)
      : terms(std::forward<Terms_>(terms_)...)
    {}

    template <class List>
    void insertFeSpaces(List& feSpaces)
    {
      for_each(terms, detail::InsertFeSpaces<List>(feSpaces));
    }

    void initElement(ElInfo const* elInfo,
                     SubAssembler* subAssembler,
                     Quadrature* quad,
                     BasisFunction const* basisFct = NULL)
    {
      for_each(terms, detail::InitElement(elInfo, subAssembler, quad, basisFct));
    }

    template <int N>
    int getDegree(int_<N>) const
    {
      return std::get<N>(terms).getDegree();
    }

    template <int N>
    typename std::tuple_element<N, decltype(Self::terms)>::type&
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
    
  private:
    std::tuple<Terms...> terms;
  };

} // end namespace AMDiS
