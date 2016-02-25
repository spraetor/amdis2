#pragma once

#include "AMDiS_fwd.hpp"
#include "Traits.hpp"

#include "expressions/TermConcepts.hpp"
#include "traits/basic.hpp"
#include "traits/traits.hpp"
#include "utility/foreach.hpp"

namespace AMDiS
{
  class SubAssembler;

  struct LazyOperatorTermBase
  {
    template <class List>
    void insertFeSpaces(List& /*feSpaces*/) const {}

    constexpr int getDegree() const
    {
      return 0;
    }

    void initElement(ElInfo const* /*elInfo*/,
                     SubAssembler* /*subAssembler*/,
                     Quadrature* /*quad*/,
                     BasisFunction const* /*basisFct*/ = NULL) {}
  };


  namespace detail
  {
    /// Functor that initializes the feSpace list
    template <class List>
    struct InsertFeSpaces
    {
      List& feSpaces;

      template <class Term>
      void operator()(Term& term)
      {
        term.insertFeSpaces(feSpaces);
      }

    private:
    };

    /// Functor that is called on each term to initialize it on an element
    struct InitElement
    {
      ElInfo const* elInfo;
      SubAssembler* subAssembler;
      Quadrature* quad;
      BasisFunction const* basisFct;

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

  private:
    std::tuple<Terms...> terms;

  public:
    template <class... Terms_,
      class = Requires_t< and_< concepts::Term<Terms_...>,
				concepts::Compatible<Types<Terms...>,
						     Types<Terms_...>> > >>
    constexpr LazyOperatorTerms(Terms_&&... terms_)
      : terms(std::forward<Terms_>(terms_)...)
    {}

    template <class List>
    void insertFeSpaces(List& feSpaces)
    {
#ifdef CXX14
      for_each([&feSpaces](auto& term)
      {
        term.insertFeSpaces(feSpaces);
      }, terms);
#else
      for_each(detail::InsertFeSpaces<List>{feSpaces}, terms);
#endif
    }

    void initElement(ElInfo const* elInfo,
                     SubAssembler* subAssembler,
                     Quadrature* quad,
                     BasisFunction const* basisFct = NULL)
    {
#ifdef CXX14
      for_each([=](auto& term)
      {
        term.insertFeSpaces(elInfo, subAssembler, quad, basisFct);
      }, terms);
#else
      for_each(detail::InitElement{elInfo, subAssembler, quad, basisFct}, terms);
#endif
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
  };

} // end namespace AMDiS
