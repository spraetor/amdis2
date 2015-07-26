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

// TODO: remove this file

/** \file functor_expr.hpp */

#ifndef AMDIS_FUNCTOR_EXPRESSION_HPP
#define AMDIS_FUNCTOR_EXPRESSION_HPP


#if HAS_VARIADIC_TEMPLATES && HAS_ALIAS_TEMPLATES
// use c++11 variant of functors
#include "functorN_expr.hpp"
  
#else
// replacement for c++98 compilers

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"
#include "Functors.h"
#include "operations/functors.hpp"

#include <boost/static_assert.hpp>

namespace AMDiS 
{
  namespace traits
  {
    template<class T>
    struct void_{ typedef void type; };

    template<typename F, typename Enable = void>
    struct functor_result_type {};
    
    template<typename F>
    struct functor_result_type<F, typename void_<typename F::value_type>::type>
    {
      typedef typename F::value_type type;
    };
    
    template<typename F>
    struct functor_result_type<F, typename void_<typename F::result_type>::type>
    {
      typedef typename F::result_type type;
    };

  } // end namespace traits

  namespace expressions 
  {
    /// Expressions for a functor with one argument
    template<typename F, typename Term>
    struct Function1 : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      BOOST_STATIC_ASSERT_MSG( (boost::is_base_of<FunctorBase, F>::value), "********** ERROR: Only functors with base FunctorBase allowed **********" );
      
      typedef typename traits::functor_result_type<F>::type value_type;    
      BOOST_STATIC_ASSERT_MSG( (!boost::is_same<value_type, traits::no_valid_type>::value), "********** ERROR: You have to define a result_type for your Functor **********" );
      
      F f;    
      Function1(const F& f_, const Term& term_) : super(term_), f(f_) {}
      
      int getDegree() const
      {
	return f.getDegree(super::term.getDegree());
      }

      inline value_type operator()(const int& iq) const { return f(super::term(iq)); }
    };
    

    /// Expressions for a functor with two arguments
    template<typename F, typename Term1, typename Term2>
    struct Function2 : public LazyOperatorTerm2<Term1, Term2>
    {
      typedef LazyOperatorTerm2<Term1, Term2> super;
      BOOST_STATIC_ASSERT_MSG( (boost::is_base_of<FunctorBase, F>::value), "********** ERROR: Only functors with base FunctorBase allowed **********" );
      
      typedef typename traits::functor_result_type<F>::type value_type;    
      BOOST_STATIC_ASSERT_MSG( (!boost::is_same<value_type, traits::no_valid_type>::value), "********** ERROR: You have to define a result_type for your Functor **********" );
      
      F f;      
      Function2(const F& f_, const Term1& term1_, const Term2& term2_) : super(term1_, term2_), f(f_) {}
      
      int getDegree() const
      {
	return f.getDegree(super::term1.getDegree(), super::term2.getDegree());
      }

      inline value_type operator()(const int& iq) const { return f(super::term1(iq), super::term2(iq)); }
    };
    

    /// Expressions for a functor with 3 arguments
    template<typename F, typename Term1, typename Term2, typename Term3>
    struct Function3 : public LazyOperatorTerm3<Term1, Term2, Term3>
    {
      typedef LazyOperatorTerm3<Term1, Term2, Term3> super;
      BOOST_STATIC_ASSERT_MSG( (boost::is_base_of<FunctorBase, F>::value), "********** ERROR: Only functors with base FunctorBase allowed **********" );
      
      typedef typename traits::functor_result_type<F>::type value_type;    
      BOOST_STATIC_ASSERT_MSG( (!boost::is_same<value_type, traits::no_valid_type>::value), "********** ERROR: You have to define a result_type for your Functor **********" );
      
      F f;      
      Function3(const F& f_, const Term1& term1_, const Term2& term2_, const Term3& term3_)
	: super(term1_, term2_, term3_), f(f_) {}
      
      int getDegree() const
      {
	return f.getDegree(super::term1.getDegree(), super::term2.getDegree(), super::term3.getDegree());
      }

      inline value_type operator()(const int& iq) const { return f(super::term1(iq), super::term2(iq), super::term3(iq)); }
    };

    /* 150203 added by Michael */
    /// Expressions for a functor with 4 arguments
    template<typename F, typename Term1, typename Term2, typename Term3, typename Term4>
    struct Function4 : public LazyOperatorTerm4<Term1, Term2, Term3, Term4>
    {
      typedef LazyOperatorTerm4<Term1, Term2, Term3, Term4> super;
      BOOST_STATIC_ASSERT_MSG( (boost::is_base_of<FunctorBase, F>::value), "********** ERROR: Only functors with base FunctorBase allowed **********" );
      
      typedef typename traits::functor_result_type<F>::type value_type;    
      BOOST_STATIC_ASSERT_MSG( (!boost::is_same<value_type, traits::no_valid_type>::value), "********** ERROR: You have to define a result_type for your Functor **********" );
      
      F f;      
      Function4(const F& f_, const Term1& term1_, const Term2& term2_, const Term3& term3_, const Term4& term4_)
	: super(term1_, term2_, term3_, term4_), f(f_) {}
      
      int getDegree() const
      {
	return f.getDegree(super::term1.getDegree(), super::term2.getDegree(), super::term3.getDegree(), super::term2.getDegree());
      }

      inline value_type operator()(const int& iq) const { return f(super::term1(iq), super::term2(iq), super::term3(iq), super::term4(iq)); }
    };
    
    /* 150203 added by Michael */
    /// Expressions for a functor with 5 arguments
    template<typename F, typename Term1, typename Term2, typename Term3, typename Term4, typename Term5>
    struct Function5 : public LazyOperatorTerm5<Term1, Term2, Term3, Term4, Term5>
    {
      typedef LazyOperatorTerm5<Term1, Term2, Term3, Term4, Term5> super;
      BOOST_STATIC_ASSERT_MSG( (boost::is_base_of<FunctorBase, F>::value), "********** ERROR: Only functors with base FunctorBase allowed **********" );
      
      typedef typename traits::functor_result_type<F>::type value_type;    
      BOOST_STATIC_ASSERT_MSG( (!boost::is_same<value_type, traits::no_valid_type>::value), "********** ERROR: You have to define a result_type for your Functor **********" );
      
      F f;      
      Function5(const F& f_, const Term1& term1_, const Term2& term2_, const Term3& term3_, const Term4& term4_, const Term5& term5_)
	: super(term1_, term2_, term3_, term4_, term5_), f(f_) {}
      
      int getDegree() const
      {
	return f.getDegree(super::term1.getDegree(), super::term2.getDegree(), super::term3.getDegree(), super::term2.getDegree(), super::term5.getDegree());
      }

      inline value_type operator()(const int& iq) const { return f(super::term1(iq), super::term2(iq), super::term3(iq), super::term4(iq), super::term5(iq)); }
    };
    
    
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
    // function with one argument
    template<typename F, typename Term>
    struct Function1 : boost::enable_if
      <
	typename traits::is_valid_arg<Term>::type,
	expressions::Function1
	<
	  F, typename traits::to_expr<Term>::type
	> 
      > {};
      
    // function with two arguments
    template<typename F, typename Term1, typename Term2>
    struct Function2 : boost::enable_if
      <
	typename boost::mpl::and_
	<
	  typename traits::is_valid_arg<Term1>::type,
	  typename traits::is_valid_arg<Term2>::type
	>::type,
	expressions::Function2
	<
	  F, 
	  typename traits::to_expr<Term1>::type,
	  typename traits::to_expr<Term2>::type
	> 
      > {};
      
      
    // function with three arguments
    template<typename F, typename Term1, typename Term2, typename Term3>
    struct Function3 : boost::enable_if
      <
	typename boost::mpl::and_
	<
	  typename traits::is_valid_arg<Term1>::type,
	  typename traits::is_valid_arg<Term2>::type,
	  typename traits::is_valid_arg<Term3>::type
	>::type,
	expressions::Function3
	<
	  F, 
	  typename traits::to_expr<Term1>::type,
	  typename traits::to_expr<Term2>::type,
	  typename traits::to_expr<Term3>::type
	> 
      > {};

    /* 150203 added by Michael */
    // function with four arguments
    template<typename F, typename Term1, typename Term2, typename Term3, typename Term4>
    struct Function4 : boost::enable_if
      <
	typename boost::mpl::and_
	<
	  typename traits::is_valid_arg<Term1>::type,
	  typename traits::is_valid_arg<Term2>::type,
	  typename traits::is_valid_arg<Term3>::type,
	  typename traits::is_valid_arg<Term4>::type
	>::type,
	expressions::Function4
	<
	  F, 
	  typename traits::to_expr<Term1>::type,
	  typename traits::to_expr<Term2>::type,
	  typename traits::to_expr<Term3>::type,
	  typename traits::to_expr<Term4>::type
	> 
      > {};      
      
    /* 150203 added by Michael */
    // function with five arguments
    template<typename F, typename Term1, typename Term2, typename Term3, typename Term4, typename Term5>
    struct Function5 : boost::enable_if
      <
	typename boost::mpl::and_
	<
	  typename traits::is_valid_arg<Term1>::type,
	  typename traits::is_valid_arg<Term2>::type,
	  typename traits::is_valid_arg<Term3>::type,
	  typename traits::is_valid_arg<Term4>::type,
	  typename traits::is_valid_arg<Term5>::type
	>::type,
	expressions::Function5
	<
	  F, 
	  typename traits::to_expr<Term1>::type,
	  typename traits::to_expr<Term2>::type,
	  typename traits::to_expr<Term3>::type,
	  typename traits::to_expr<Term4>::type,
	  typename traits::to_expr<Term5>::type
	> 
      > {};
      
  } // end namespace result_of


  // call a function with 1 argument
  // _____________________________________________________________________________
  template<typename F, typename Term>
  inline typename result_of::Function1<F, Term>::type
  function_(const F& f, const Term& t)
  { 
    typedef typename traits::to_expr<Term>::to Expr;
    return expressions::Function1<F, typename Expr::type>(f, Expr::get(t)); 
  }


  template<typename F, typename Term>
  inline typename result_of::Function1<F, Term>::type
  function_(const Term& t) 
  { 
    typedef typename traits::to_expr<Term>::to Expr;
    return expressions::Function1<F, typename Expr::type>(F(), Expr::get(t)); 
  }


  // call a function with 2 arguments
  // _____________________________________________________________________________
  template<typename F, typename Term1, typename Term2>
  inline typename result_of::Function2<F, Term1, Term2>::type
  function_(const F& f, const Term1& t1, const Term2& t2) 
  { 
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    return expressions::Function2<F, typename Expr1::type, typename Expr2::type>
	    (f, Expr1::get(t1), Expr2::get(t2)); 
  }


  template<typename F, typename Term1, typename Term2>
  inline typename result_of::Function2<F, Term1, Term2>::type
  function_(const Term1& t1, const Term2& t2) 
  { 
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    return expressions::Function2<F, typename Expr1::type, typename Expr2::type>
	    (F(), Expr1::get(t1), Expr2::get(t2)); 
  }


  // call a function with 3 arguments
  // _____________________________________________________________________________
  template<typename F, typename Term1, typename Term2, typename Term3>
  inline typename result_of::Function3<F, Term1, Term2, Term3>::type
  function_(const F& f, const Term1& t1, const Term2& t2, const Term3& t3) 
  { 
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    typedef typename traits::to_expr<Term3>::to Expr3;
    return expressions::Function3<F, typename Expr1::type, typename Expr2::type, typename Expr3::type>
	    (f, Expr1::get(t1), Expr2::get(t2), Expr3::get(t3)); 
  }


  template<typename F, typename Term1, typename Term2, typename Term3>
  inline typename result_of::Function3<F, Term1, Term2, Term3>::type
  function_(const Term1& t1, const Term2& t2, const Term3& t3) 
  { 
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    typedef typename traits::to_expr<Term3>::to Expr3;
    return expressions::Function3<F, typename Expr1::type, typename Expr2::type, typename Expr3::type>
	    (F(), Expr1::get(t1), Expr2::get(t2), Expr3::get(t3)); 
  }

  /* 150203 added by Michael */
  // call a function with 4 arguments
  // _____________________________________________________________________________
  template<typename F, typename Term1, typename Term2, typename Term3, typename Term4>
  inline typename result_of::Function4<F, Term1, Term2, Term3, Term4>::type
  function_(const F& f, const Term1& t1, const Term2& t2, const Term3& t3, const Term4& t4) 
  { 
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    typedef typename traits::to_expr<Term3>::to Expr3;
    typedef typename traits::to_expr<Term4>::to Expr4;
    return expressions::Function4<F, typename Expr1::type, typename Expr2::type, typename Expr3::type, typename Expr4::type>
	    (f, Expr1::get(t1), Expr2::get(t2), Expr3::get(t3), Expr4::get(t4)); 
  }


  template<typename F, typename Term1, typename Term2, typename Term3, typename Term4>
  inline typename result_of::Function4<F, Term1, Term2, Term3, Term4>::type
  function_(const Term1& t1, const Term2& t2, const Term3& t3, const Term4& t4) 
  { 
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    typedef typename traits::to_expr<Term3>::to Expr3;
    typedef typename traits::to_expr<Term4>::to Expr4;
    return expressions::Function4<F, typename Expr1::type, typename Expr2::type, typename Expr3::type, typename Expr4::type>
	    (F(), Expr1::get(t1), Expr2::get(t2), Expr3::get(t3), Expr4::get(t4)); 
  }
  
  /* 150203 added by Michael */
  // call a function with 5 arguments
  // _____________________________________________________________________________
  template<typename F, typename Term1, typename Term2, typename Term3, typename Term4, typename Term5>
  inline typename result_of::Function5<F, Term1, Term2, Term3, Term4, Term5>::type
  function_(const F& f, const Term1& t1, const Term2& t2, const Term3& t3, const Term4& t4, const Term5& t5) 
  { 
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    typedef typename traits::to_expr<Term3>::to Expr3;
    typedef typename traits::to_expr<Term4>::to Expr4;
    typedef typename traits::to_expr<Term5>::to Expr5;
    return expressions::Function5<F, typename Expr1::type, typename Expr2::type, typename Expr3::type, typename Expr4::type, typename Expr5::type>
	    (f, Expr1::get(t1), Expr2::get(t2), Expr3::get(t3), Expr4::get(t4), Expr5::get(t5)); 
  }


  template<typename F, typename Term1, typename Term2, typename Term3, typename Term4, typename Term5>
  inline typename result_of::Function5<F, Term1, Term2, Term3, Term4, Term5>::type
  function_(const Term1& t1, const Term2& t2, const Term3& t3, const Term4& t4, const Term5& t5) 
  { 
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    typedef typename traits::to_expr<Term3>::to Expr3;
    typedef typename traits::to_expr<Term4>::to Expr4;
    typedef typename traits::to_expr<Term5>::to Expr5;
    return expressions::Function5<F, typename Expr1::type, typename Expr2::type, typename Expr3::type, typename Expr4::type, typename Expr5::type>
	    (F(), Expr1::get(t1), Expr2::get(t2), Expr3::get(t3), Expr4::get(t4), Expr5::get(t5)); 
  }


  // function wrapper for abstract functions
  // _____________________________________________________________________________
  template<typename TOut, typename TIn>
  inline expressions::Wrapper<TOut,TIn> wrap(AbstractFunction<TOut, TIn>* fct) 
  { return expressions::Wrapper<TOut,TIn>(fct); }

} // end namespace AMDiS

#endif

#endif // AMDIS_FUNCTOR_EXPRESSION_HPP
