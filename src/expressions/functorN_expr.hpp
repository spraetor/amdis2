/** \file functorN_expr.h */

#pragma once

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"
// #include "Functors.h"
// #include "expressions/functor_expr.hpp"
#include "operations/functors.hpp"

#include <boost/static_assert.hpp>

#include <tuple>
#include <utility> 
#include <functional>

namespace AMDiS 
{

  template<int I>
  using int_ = std::integral_constant<int, I>;
  
  
  namespace traits
  {
    /// get the degree of a functor by combining the degrees of the arguments
    template<typename F, typename Enable = void>
    struct functor_degree
    {
      template<typename... Int>
      static int eval(F f, Int... d) { return 0; }
    };
    
    template<typename F>
    struct functor_degree<F, typename enable_if< boost::is_base_of<FunctorBase, F> >::type >
    {
      template<typename... Int>
      static int eval(F f, Int... d) { return f.getDegree(d...); }
    };
    
  } // end namespace traits
  
  namespace result_of
  {    
    /// extract result type from function pointers
    template <class FPtr>
    struct Function;

    template <class R, class C, class... As>
    struct Function<R (C::*)(As...)>
    {
	typedef R result_type;
    };

    template<class R, class C, class... As>
    struct Function<R (C::*)(As...) const>
    {
	typedef R type;
    };
    
    template<class T>
    typename Function<T>::type function_helper(T);

    /// extract result type from functors
    template <class F>
    struct Functor
    {
      typedef decltype(function_helper(&F::operator())) type;
    };
    
    template <class R, class... As>
    struct Functor<std::function<R(As...)> >
    {
      typedef R type;
    };
    
  } // end namespace result_of

    
  // the expressions
  // _____________________________________________________________________________

  namespace expressions 
  {
    /// Functor that takes arbitrary number of arguments
    template<typename F, typename... Terms>
    struct FunctionN : public LazyOperatorTerms<Terms...>
    {
      typedef LazyOperatorTerms<Terms...> super;
      static const int N = sizeof...(Terms);
      
      typedef typename result_of::Functor<F>::type value_type;
      BOOST_STATIC_ASSERT_MSG( (!boost::is_same<value_type, traits::no_valid_type>::value), "********** ERROR: You have to define a result_type for your Functor **********" );
      
      F f; ///< the functor
      
      template<typename... Terms_>
      FunctionN(const F& f_, Terms_... terms_)
	: super(terms_...), f(f_) {}
      
      // call f.getDegree() function    
      template<int I, typename... Terms_>
      int getDegree(int_<I>, const Terms_&... terms) const
      {
	return getDegree(int_<I-1>(), std::get<I-1>(super::term_tuple), terms...);
      }
      
      template<typename... Terms_>
      int getDegree(int_<0>, const Terms_&... terms) const
      {
	return traits::functor_degree<F>::eval(f, terms.getDegree()...);
      }
      
      int getDegree() const
      {
	return getDegree(int_<N>());
      }

      // call f.operator()(...)
      template<int I, typename... Terms_>
      value_type eval(const int& iq, int_<I>, const Terms_&... terms) const
      {
	       return eval(iq, int_<I-1>(), std::get<I-1>(super::term_tuple), terms...);
      }
      
      template<typename... Terms_>
      value_type eval(const int& iq, int_<0>, Terms_... terms) const
      {
	       return f(terms(iq)...);  // f(term1(iq), term2(iq), term3(iq),...)
      }
      
      value_type operator()(const int& iq) const { return eval(iq, int_<N>()); }
    };
    
    template<typename F, typename Term>
    using Function1 = FunctionN<F, Term>;
    
    template<typename F, typename Term1, typename Term2>
    using Function2 = FunctionN<F, Term1, Term2>;
    
    template<typename F, typename Term1, typename Term2, typename Term3>
    using Function3 = FunctionN<F, Term1, Term2, Term3>;
    
    template<typename F, typename Term1, typename Term2, typename Term3, typename Term4>
    using Function4 = FunctionN<F, Term1, Term2, Term3, Term4>;
    
#if 0
    /// A wrapper functor for AMDiS::AbstractFunctions
    template<typename TOut, typename TIn>
    struct Wrapper : public FunctorBase
    {
      typedef TOut result_type;
      Wrapper(AbstractFunction<TOut, TIn>* fct_) : fct(fct_) {}
      int getDegree(int degree) const { return fct->getDegree(); }
      
      TOut operator()(const TIn& x) const 
      {
	return (*fct)(x);
      }
      
    protected:
      AbstractFunction<TOut, TIn>* fct;
    };
#endif
  } // end namespace expressions


  namespace result_of
  {
    // result of generator-functions (used for enable_if constructs)
    template <typename F, typename... Terms>
    struct FunctionN : std::enable_if
      <
	and_< typename traits::is_valid_arg<Terms>::type... >::value,
	expressions::FunctionN< F, typename traits::to_expr<Terms>::type...> 
      > {};
      
    template <typename F, typename Term>
    struct FunctionN<F, Term> : std::enable_if
      <
	traits::is_valid_arg<Term>::value,
	expressions::FunctionN< F, typename traits::to_expr<Term>::type> 
      > {};
      
  } // end namespace result_of

  // generator-functions
  // _____________________________________________________________________________

  template<typename F, typename... Terms>
  inline typename result_of::FunctionN<F, Terms...>::type
  function_(F&& f, Terms... ts) 
  {
    return expressions::FunctionN<F, typename traits::to_expr<Terms>::to::type...>
	    (std::forward<F>(f), traits::to_expr<Terms>::to::get(ts)...); 
  }

  template<typename F, typename... Terms>
  inline typename result_of::FunctionN<F, Terms...>::type
  func(F&& f, Terms... ts) 
  {
    return expressions::FunctionN<F, typename traits::to_expr<Terms>::to::type...>
	    (std::forward<F>(f), traits::to_expr<Terms>::to::get(ts)...); 
  }

  template<typename F, typename Term0, typename... Terms>
  inline typename result_of::FunctionN<F, Term0, Terms...>::type
  eval(F&& f, Term0 t0, Terms... ts) 
  {
    return expressions::FunctionN<F, typename traits::to_expr<Term0>::to::type, 
				     typename traits::to_expr<Terms>::to::type...>
      (std::forward<F>(f), traits::to_expr<Term0>::to::get(t0), traits::to_expr<Terms>::to::get(ts)...); 
  }

#if 0
  // function wrapper for abstract functions
  // _____________________________________________________________________________
  template<typename TOut, typename TIn>
  inline expressions::Wrapper<TOut,TIn> wrap(AbstractFunction<TOut, TIn>* fct) 
  { return expressions::Wrapper<TOut,TIn>(fct); }
#endif

} // end namespace AMDiS
