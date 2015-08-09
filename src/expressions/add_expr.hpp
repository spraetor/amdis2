/** \file add_expr.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"

#include <string>
#include <boost/static_assert.hpp>

namespace AMDiS 
{
  namespace expressions 
  {
    /// Expression that represents "E1 + E2"
    template <class Term1, class Term2>
    struct Add : public LazyOperatorTerms<Term1, Term2>
    {
      using Super = LazyOperatorTerms<Term1, Term2>;
      using value_type = typename traits::add_type< Value_t<Term1>, Value_t<Term2> >::type;
      
      static_assert( !(std::is_same<value_type, traits::no_valid_type>::value), 
                     "********** ERROR: Can not add terms **********" );
      
      template <class Term0_, class Term1_>
      Add(Term0_&& term0_, Term1_&& term1_)
        : Super(std::forward<Term1>(term0_), std::forward<Term2>(term1_)) 
      {}
      
      int getDegree() const
      {
        return std::max(Super::getDegree(_0), Super::getDegree(_1));
      }

      value_type operator()(const int& iq) const 
      { 
        return Super::getTerm(_0)(iq) + Super::getTerm(_1)(iq); 
      }
      
      std::string str() const 
      { 
        return std::string("(") + Super::getTerm(_0).str() + " + " + Super::getTerm(_1).str() + ")"; 
      }
    };
    
    
    /// Expression that represents "E1 - E2"
    template <class Term1, class Term2>
    struct Subtract : public LazyOperatorTerms<Term1, Term2>
    {
      using Super = LazyOperatorTerms<Term1, Term2>;
      using value_type = typename traits::add_type<Value_t<Term1>, Value_t<Term2> >::type;
      
      static_assert( !(std::is_same<value_type, traits::no_valid_type>::value), 
                     "********** ERROR: Can not subtract terms **********" );
      
      template <class Term0_, class Term1_>
      Subtract(Term0_&& term0_, Term1_&& term1_)
        : Super(std::forward<Term1>(term0_), std::forward<Term2>(term1_))
      { }
      
      int getDegree() const
      {
        return std::max(Super::getDegree(_0), Super::getDegree(_1));
      }

      value_type operator()(const int& iq) const 
      {
        return Super::getTerm(_0)(iq) - Super::getTerm(_1)(iq); 
      }
      
      std::string str() const 
      { 
        return std::string("(") + Super::getTerm(_0).str() + " - " + Super::getTerm(_1).str() + ")"; 
      }
    };
    
    
    /// Expression that represents "-E"
    template <class Term>
    struct Negative : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      
      template <class Term_>
      Negative(Term_&& term_)
        : Super(std::forward<Term>(term_)) 
      { }
      
      int getDegree() const
      {
        return Super::getDegree(_0);
      }

      value_type operator()(const int& iq) const 
      { 
        return -Super::getTerm(_0)(iq); 
      }
      
      std::string str() const 
      { 
        return std::string("-(") + Super::getTerm(_0).str() + ")"; 
      }
    };
    
  } // end namespace expressions


  namespace result_of
  {
    template <class Term1, class Term2>
    using Add = enable_if
      < 
      	typename traits::is_valid_arg2<Term1, Term2>::type,
      	expressions::Add
      	<
      	  typename traits::to_expr<Term1>::type, 
      	  typename traits::to_expr<Term2>::type
      	>
      >;
      
      
    template <class Term>
    using Negative = enable_if
      < 
      	typename traits::is_valid_arg<Term>::type,
      	expressions::Negative
      	<
      	  typename traits::to_expr<Term>::type
      	>
      >;
      
      
    template <class Term1, class Term2>
    using Subtract = enable_if
      < 
      	typename traits::is_valid_arg2<Term1, Term2>::type,
      	expressions::Subtract
      	<
      	  typename traits::to_expr<Term1>::type, 
      	  typename traits::to_expr<Term2>::type
      	>
      >;
  }

#if 0
  // add two terms
  // _____________________________________________________________________________
  template<typename Term1, typename Term2>
  inline typename result_of::Add<Term1, Term2>::type
  operator+(Term1&& t1, Term2&& t2)
  {
    using Expr1 = traits::to_expr<Term1>;
    using Expr2 = traits::to_expr<Term2>;
    return {Expr1::get(std::forward<Term1>(t1)), Expr2::get(std::forward<Term2>(t2))};
  }


  // negative of a term
  // _____________________________________________________________________________
  template<typename Term>
  inline typename result_of::Negative<Term>::type
  operator-(Term&& t)
  {
    using Expr = traits::to_expr<Term>;
    return {Expr::get(std::forward<Term>(t))};
  }


  // substract two terms
  // _____________________________________________________________________________
  template<typename Term1, typename Term2>
  inline typename result_of::Subtract<Term1, Term2>::type
  operator-(Term1&& t1, Term2&& t2)
  {
    using Expr1 = traits::to_expr<Term1>;
    using Expr2 = traits::to_expr<Term2>;
    return {Expr1::get(std::forward<Term1>(t1)), Expr2::get(std::forward<Term2>(t2))};
  }
#endif
} // end namespace AMDiS
