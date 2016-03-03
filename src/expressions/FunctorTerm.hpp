#pragma once

// std c++ headers
#include <tuple>
#include <utility>

// AMDiS headers
#include "expressions/LazyOperatorTerm.hpp"
#include "expressions/TermConcepts.hpp"
#include "expressions/ComponentView.hpp"
#include "traits/basic.hpp"
#include "traits/traits_fwd.hpp"
#include "traits/traits.hpp"
#include "traits/meta_basic.hpp"

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

  // forward declaration
  template <class F, class Term1, class... Terms>
  struct FunctorTerm;

  // extract the shape of the term from the shape of the functor result-type
  template <class F, class Term1, class... Terms>
  struct FunctorShape
  {
    using value_type = typename std::result_of<
        F(typename Term1::value_type, typename Terms::value_type...)
      >::type;
    using type = ShapedTerm_t<value_type, FunctorTerm<F, Term1, Terms...>>;
  };

  template <class F, class Term1, class... Terms>
  using FunctorShape_t = typename FunctorShape<F, Term1, Terms...>::type;


  // the expressions
  // ___________________________________________________________________________

  /// Functor that takes arbitrary number of arguments
  template <class F, class Term1, class... Terms>
  struct FunctorTerm
    : public FunctorShape_t<F, Term1, Terms...>,
      public LazyOperatorTerms<Term1, Terms...>/*,
      public ComponentView<Value_t<FunctorShape<F, Term1, Terms...>>,
                           FunctorTerm<F, Term1, Terms...>>*/
  {
    using Self       = FunctorTerm;
    using Super      = LazyOperatorTerms<Term1, Terms...>;
    using value_type = Value_t< FunctorShape<F, Term1, Terms...> >;

    FunctorTerm(Term1 const& term1_, Terms const&... terms_)
      : Super(term1_, terms_...),
        fct{}
    {}

    template <class F_,
      class = Requires_t<traits::IsCompatible<F, F_>> >
    FunctorTerm(F_&& f_, Term1 const& term1_, Terms const&... terms_)
      : Super(term1_, terms_...),
        fct(std::forward<F_>(f_))
    {}

    /// return the required quadrature degree to integrate term
    constexpr int getDegree() const
    {
      return getDegree(int_<N>());
    }

    /// eval at point with index iq
    value_type evalAtIdx(int iq) const
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
    int getDegree(int_<I>, Terms_&&... terms) const
    {
      return getDegree(int_<I-1>(),
                       Super::getTerm(int_<I-1>()),
                       std::forward<Terms_>(terms)...);
    }

    template <class... Terms_>
    int getDegree(int_<0>, Terms_ const&... terms) const
    {
      return traits::functor_degree<F,N>::eval(fct, terms.getDegree()...);
    }

    // call f.operator()(...)
    template <class Arg, int I, class... Terms_>
    value_type eval(Arg&& arg, int_<I>,
                    Terms_&& ... terms) const
    {
      return eval(std::forward<Arg>(arg),
                  int_<I-1>(),
                  Super::getTerm(int_<I-1>()),
                  std::forward<Terms_>(terms)...);
    }

    template <class... Terms_>
    value_type eval(int iq, int_<0>, Terms_ const&... terms) const
    {
      return fct(terms.evalAtIdx(iq)...);  // f(t1(iq), t2(iq), t3(iq),...)
    }

    template <class... Terms_>
    value_type eval(WorldVector<double> const& x,
                    int_<0>,
                    Terms_ const&... terms) const
    {
      return fct(terms(x)...);  // f(t1(iq), t2(iq), t3(iq),...)
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
      using value_type = Value_t< FunctorShape<F, Term1, Terms...> >;
      using size_type  = int;
    };
    /// \endcond

  } // end namespace traits

} // end namespace AMDiS
