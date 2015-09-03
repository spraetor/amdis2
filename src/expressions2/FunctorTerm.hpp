/** \file functorN_expr.h */

#pragma once

// std c++ headers
#include <tuple>
#include <utility>

// AMDiS headers
#include <expressions/LazyOperatorTerm.h>
#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>
#include <traits/meta_basic.hpp>

#include "TermConcepts.hpp"

namespace AMDiS
{
  namespace traits
  {
    /// get the degree of a functor by combining the degrees of the arguments
    template <class F, int N, class = void>
    struct functor_degree
    {
      template <class... Int>
      static constexpr int eval(F, Int...)
      {
        return 0;
      }
    };

    template <class F, int N>
    struct functor_degree<F, N, requires::TermFunctor<F,N>>
    {
      template <class... Int>
      static constexpr int eval(F const& f, Int... d)
      {
        return f.getDegree(d...);
      }
    };

  } // end namespace traits


  // the expressions
  // _____________________________________________________________________________

  /// Functor that takes arbitrary number of arguments
  template <class F, class Term1, class... Terms>
  struct FunctorTerm
    : public ShapedTerm_t<Term1, FunctorTerm<F, Term1, Terms...>>,
      public LazyOperatorTerms<Term1, Terms...>
  {
    using Self       = FunctorTerm;
    using Super      = LazyOperatorTerms<Term1, Terms...>;
    using value_type = typename std::result_of<F(Value_t<Term1>, Value_t<Terms>...)>::type;

    template <class... Terms_,
              class = Requires_t<traits::IsCompatible<Types<Term1,Terms...>, Types<Terms_...>>>>
                                 FunctorTerm(Terms_&&... terms_)
                                   : Super(std::forward<Terms_>(terms_)...),
                                     fct{}
    {}

    template <class F_, class... Terms_,
              class = Requires_t<traits::IsCompatible<Types<F,Term1,Terms...>, Types<F_,Terms_...>>>>
                                 FunctorTerm(F_&& f, Terms_&&... terms_)
                                   : Super(std::forward<Terms_>(terms_)...),
                                     fct(f)
    {}

    /// return the required quadrature degree to integrate term
    constexpr int getDegree() const
    {
      return getDegree(int_<N>());
    }

    /// eval at point with index iq
    value_type operator[](int iq) const
    {
      return eval(iq, int_<N>());
    }

    /// eval at point with coordinate x
    value_type operator()(WorldVector<double> const& x) const
    {
      return eval(x, int_<N>());
    }

  private:
    static constexpr int N = sizeof...(Terms) + 1;

    // call f.getDegree() function
    template <int I, class... Terms_>
    int getDegree(int_<I>, Terms_&& ... terms) const
    {
      return getDegree(int_<I-1>(), Super::getTerm(int_<I-1>()),
                       std::forward<Terms_>(terms)...);
    }

    template <class... Terms_>
    int getDegree(int_<0>, Terms_ const& ... terms) const
    {
      return traits::functor_degree<F,N>::eval(fct, terms.getDegree()...);
    }

    // call f.operator()(...)
    template <class Arg, int I, class... Terms_>
    value_type eval(Arg&& arg, int_<I>, Terms_&& ... terms) const
    {
      return eval(std::forward<Arg>(arg), int_<I-1>(), Super::getTerm(int_<I-1>()),
                  std::forward<Terms_>(terms)...);
    }

    template <class... Terms_>
    value_type eval(int iq, int_<0>, Terms_ const& ... terms) const
    {
      return fct(terms[iq]...);  // f(term1(iq), term2(iq), term3(iq),...)
    }

    template <class... Terms_>
    value_type eval(WorldVector<double> const& x, int_<0>, Terms_ const& ... terms) const
    {
      return fct(terms(x)...);  // f(term1(iq), term2(iq), term3(iq),...)
    }

  private:
    F fct; ///< the functor
  };


  namespace traits
  {
    /// \cond HIDDEN_SYMBOLS
    template <class F, class Term1, class... Terms>
    struct category<FunctorTerm<F, Term1, Terms...>>
    {
      using tag        = typename category<Term1>::tag;
      using value_type = typename std::result_of<F(Value_t<Term1>, Value_t<Terms>...)>::type;
      using size_type  = int;
                                                          };
    /// \endcond
                                                  }

} // end namespace AMDiS
