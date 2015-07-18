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



/** \file functorN_expr.h */

#ifndef AMDIS_FUNCTOR_N_EXPRESSION_HPP
#define AMDIS_FUNCTOR_N_EXPRESSION_HPP

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"
#include "Functors.h"
// #include "expressions/functor_expr.hpp"
#include "operations/functors.hpp"

#include <boost/static_assert.hpp>

#include <tuple>
#include <utility> 
#include <functional>

namespace AMDiS 
{
  /// for_each for std::tuple
  template<std::size_t I = 0, typename FuncT, typename... Tp>
  inline typename std::enable_if<I == sizeof...(Tp), void>::type
  for_each(std::tuple<Tp...> &, FuncT) { }

  template<std::size_t I = 0, typename FuncT, typename... Tp>
  inline typename std::enable_if<I < sizeof...(Tp), void>::type
  for_each(std::tuple<Tp...>& t, FuncT f)
  {
    f(std::get<I>(t));
    for_each<I + 1, FuncT, Tp...>(t, f);
  }

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
    
    // specialization for abstract-functions
    template<typename R, typename... Args>
    struct functor_degree<AbstractFunction<R, Args...> >
    {
      template<typename... Int>
      static int eval(AbstractFunction<R, Args...> const& fct, Int... d) { return fct.getDegree(); }
    };
    
    template<typename R, typename Arg0, typename Arg1>
    struct functor_degree<BinaryAbstractFunction<R, Arg0, Arg1> >
    {
      template<typename... Int>
      static int eval(BinaryAbstractFunction<R, Arg0, Arg1> const& fct, Int... d) { return fct.getDegree(); }
    };
    
    template<typename R, typename Arg0, typename Arg1, typename Arg2>
    struct functor_degree<TertiaryAbstractFunction<R, Arg0, Arg1, Arg2> >
    {
      template<typename... Int>
      static int eval(TertiaryAbstractFunction<R, Arg0, Arg1, Arg2> const& fct, Int... d) { return fct.getDegree(); }
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
  

  namespace detail 
  {
    /// Functor that initializes the feSpace list
    template<typename List>
    struct InsertFeSpaces 
    {
      List& feSpaces;
      InsertFeSpaces(List& feSpaces_) : feSpaces(feSpaces_) {};
      
      template<typename Term>
      void operator()(Term& term) {
	term.insertFeSpaces(feSpaces);
      }
    };
    
    /// Functor that is called on each term to initialize it on an element
    template<typename OT>
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
  template<typename... Terms>
  struct LazyOperatorTerms : public LazyOperatorTermBase
  {
    std::tuple<Terms...> term_tuple;
    
    template<typename... Terms_>
    LazyOperatorTerms(Terms_... terms_)
    : term_tuple(terms_...) {}
    
    template<typename List>
    inline void insertFeSpaces(List& feSpaces)
    {
      for_each(term_tuple, detail::InsertFeSpaces<List>(feSpaces));
    }

    template<typename OT>
    inline void initElement(OT* ot, const ElInfo* elInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      for_each(term_tuple, detail::InitElement<OT>(ot, elInfo, subAssembler, quad, basisFct));
    }

    template<typename OT>
    inline void initElement(OT* ot, const ElInfo* smallElInfo, const ElInfo* largeElInfo,
		    SubAssembler* subAssembler, Quadrature *quad, 
		    const BasisFunction *basisFct = NULL)
    {
      for_each(term_tuple, detail::InitElement<OT>(ot, smallElInfo, largeElInfo, subAssembler, quad, basisFct));
    }
    
    inline double operator()(const int& iq) const;
  };

    
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
      inline value_type eval(const int& iq, int_<I>, const Terms_&... terms) const
      {
	return eval(iq, int_<I-1>(), std::get<I-1>(super::term_tuple), terms...);
      }
      
      template<typename... Terms_>
      inline value_type eval(const int& iq, int_<0>, Terms_... terms) const
      {
	return f(terms(iq)...);  // f(term1(iq), term2(iq), term3(iq),...)
      }
      
      inline value_type operator()(const int& iq) const { return eval(iq, int_<N>()); }
    };
    
    template<typename F, typename Term>
    using Function1 = FunctionN<F, Term>;
    
    template<typename F, typename Term1, typename Term2>
    using Function2 = FunctionN<F, Term1, Term2>;
    
    template<typename F, typename Term1, typename Term2, typename Term3>
    using Function3 = FunctionN<F, Term1, Term2, Term3>;
    
    template<typename F, typename Term1, typename Term2, typename Term3, typename Term4>
    using Function4 = FunctionN<F, Term1, Term2, Term3, Term4>;
    
    
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


  // function wrapper for abstract functions
  // _____________________________________________________________________________
  template<typename TOut, typename TIn>
  inline expressions::Wrapper<TOut,TIn> wrap(AbstractFunction<TOut, TIn>* fct) 
  { return expressions::Wrapper<TOut,TIn>(fct); }

} // end namespace AMDiS

#endif // AMDIS_FUNCTOR_N_EXPRESSION_HPP
