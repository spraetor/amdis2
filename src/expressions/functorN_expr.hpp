/** \file functorN_expr.h */

#pragma once

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"
// #include "Functors.h"
// #include "expressions/functor_expr.hpp"
#include "operations/functors.hpp"
#include "traits/meta_basic.hpp"

#include <tuple>
#include <utility>
#include <functional>

namespace AMDiS
{
  namespace traits
  {
    /// get the degree of a functor by combining the degrees of the arguments
    template <class F, class Enable = void>
    struct functor_degree
    {
      template<typename... Int>
      static int eval(F f, Int... d)
      {
        return 0;
      }
    };

    template <class F>
    struct functor_degree<F, typename enable_if<std::is_base_of<FunctorBase, F>>::type>
    {
      template <class... Int>
      static int eval(F f, Int... d)
      {
        return f.getDegree(d...);
      }
    };

  } // end namespace traits

  namespace result_of
  {
    /// extract result type from function pointers
    template <class FPtr>
    struct Function;

    template <class R, class C, class... Args>
    struct Function<R (C::*)(Args...)>
    {
      typedef R type;
    };

    template <class R, class C, class... Args>
    struct Function<R (C::*)(Args...) const>
    {
      typedef R type;
    };

    template <class T>
    typename Function<T>::type function_helper(T);

    /// extract result type from functors
    template <class F>
    struct Functor
    {
      typedef decltype(function_helper(&F::operator())) type;
    };

    template <class R, class... Args>
    struct Functor<std::function<R(Args...)>>
                                           {
                                             typedef R type;
                                           };

  } // end namespace result_of


  // the expressions
  // _____________________________________________________________________________

  namespace expressions
  {
    /// Functor that takes arbitrary number of arguments
    template <class F, class... Terms>
    struct FunctionN : public LazyOperatorTerms<Terms...>
    {
      using Super = LazyOperatorTerms<Terms...>;
      static constexpr int N = sizeof...(Terms);

      using value_type = typename result_of::Functor<typename std::decay<F>::type>::type;
      static_assert( (!std::is_same<value_type, traits::no_valid_type>::value),
                     "********** ERROR: You have to define a result_type for your Functor **********" );

      F f; ///< the functor

      template <class... Terms_>
      FunctionN(const F& f_, Terms_&& ... terms_)
        : Super(std::forward<Terms_>(terms_)...), f(f_) {}

      // call f.getDegree() function
      template <int I, class... Terms_>
      int getDegree(int_<I>, const Terms_&... terms) const
      {
        return getDegree(int_<I-1>(), Super::getTerm(int_<I-1>()), terms...);
      }

      template <class... Terms_>
      int getDegree(int_<0>, const Terms_&... terms) const
      {
        return traits::functor_degree<F>::eval(f, terms.getDegree()...);
      }

      int getDegree() const
      {
        return getDegree(int_<N>());
      }

      // call f.operator()(...)
      template <int I, class... Terms_>
      value_type eval(const int& iq, int_<I>, const Terms_&... terms) const
      {
        return eval(iq, int_<I-1>(), Super::getTerm(int_<I-1>()), terms...);
      }

      template <class... Terms_>
      value_type eval(const int& iq, int_<0>, Terms_... terms) const
      {
        return f(terms(iq)...);  // f(term1(iq), term2(iq), term3(iq),...)
      }

      value_type operator()(const int& iq) const
      {
        return eval(iq, int_<N>());
      }
    };

    template <class F, class Term>
    using Function1 = FunctionN<F, Term>;

    template <class F, class Term1, class Term2>
    using Function2 = FunctionN<F, Term1, Term2>;

    template <class F, class Term1, class Term2, class Term3>
    using Function3 = FunctionN<F, Term1, Term2, Term3>;

    template <class F, class Term1, class Term2, class Term3, class Term4>
    using Function4 = FunctionN<F, Term1, Term2, Term3, Term4>;

  } // end namespace expressions


  namespace result_of
  {
    // result of generator-functions (used for enable_if constructs)
    template <class F, class... Terms>
    using FunctionN = std::enable_if
                      <
                      and_<traits::is_valid_arg<Terms>...>::value,
                      expressions::FunctionN<F, typename traits::to_expr<Terms>::type...>
                      >;

  } // end namespace result_of

  // generator-functions
  // _____________________________________________________________________________

  template <class F, class... Terms>
  inline typename result_of::FunctionN<F, Terms...>::type
  function_(F&& f, Terms... ts)
  {
    return {std::forward<F>(f), traits::to_expr<Terms>::get(ts)...};
  }

  template <class F, class... Terms>
  inline typename result_of::FunctionN<F, Terms...>::type
  func(F&& f, Terms... ts)
  {
    return {std::forward<F>(f), traits::to_expr<Terms>::get(ts)...};
  }

#if 0
  template <class F, class Term0, class... Terms>
  inline typename result_of::FunctionN<F, Term0, Terms...>::type
  eval(F&& f, Term0 t0, Terms... ts)
  {
    return expressions::FunctionN<F, typename traits::to_expr<Term0>::type,
           typename traits::to_expr<Terms>::type...>
           (std::forward<F>(f), traits::to_expr<Term0>::get(t0), traits::to_expr<Terms>::get(ts)...);
  }
#endif

} // end namespace AMDiS
