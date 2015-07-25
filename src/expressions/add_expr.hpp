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
      
      Add(Term1&& term1, Term2&& term2)
        : Super(std::forward<Term1>(term1), std::forward<Term2>(term2)) 
      { }
      
      int getDegree() const
      {
        return std::max(degree<0>(*this), degree<1>(*this));
      }

      value_type operator()(const int& iq) const 
      { 
        return term<0>(*this)(iq) + term<1>(*this)(iq); 
      }
      
      std::string str() const 
      { 
        return std::string("(") + term<0>(*this).str() + " + " + term<1>(*this).str() + ")"; 
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
      
      Subtract(Term1&& term1_, Term2&& term2_)
        : Super(std::forward<Term1>(term1_), std::forward<Term2>(term2_))
      { }
      
      int getDegree() const
      {
        return std::max(degree<0>(*this), degree<1>(*this));
      }

      value_type operator()(const int& iq) const 
      {
        return term<0>(*this)(iq) - term<1>(*this)(iq); 
      }
      
      std::string str() const 
      { 
        return std::string("(") + term<0>(*this).str() + " - " + term<1>(*this).str() + ")"; 
      }
    };
    
    
    /// Expression that represents "-E"
    template <class Term>
    struct Negative : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      
      Negative(Term&& term_)
        : Super(std::forward<Term>(term_)) 
      { }
      
      int getDegree() const
      {
        return degree<0>(*this);
      }

      value_type operator()(const int& iq) const 
      { 
        return -term<0>(*this)(iq); 
      }
      
      std::string str() const 
      { 
        return std::string("-(") + term<0>(*this).str() + ")"; 
      }
    };
    
  } // end namespace expressions


  namespace result_of
  {
    template <class Term1, class Term2>
    using Add = boost::enable_if
      < 
      	typename traits::is_valid_arg2<Term1, Term2>::type,
      	expressions::Add
      	<
      	  typename traits::to_expr<Term1>::type, 
      	  typename traits::to_expr<Term2>::type
      	>
      >;
      
      
    template <class Term>
    using Negative = boost::enable_if
      < 
      	typename traits::is_valid_arg<Term>::type,
      	expressions::Negative
      	<
      	  typename traits::to_expr<Term>::type
      	>
      >;
      
      
    template <class Term1, class Term2>
    using Subtract = boost::enable_if
      < 
      	typename traits::is_valid_arg2<Term1, Term2>::type,
      	expressions::Subtract
      	<
      	  typename traits::to_expr<Term1>::type, 
      	  typename traits::to_expr<Term2>::type
      	>
      >;
  }


  // add two terms
  // _____________________________________________________________________________
  template<typename Term1, typename Term2>
  inline typename result_of::Add<Term1, Term2>::type
  operator+(const Term1& t1, const Term2& t2)
  {
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    return {Expr1::get(t1), Expr2::get(t2)};
  }


  // negative of a term
  // _____________________________________________________________________________
  template<typename Term>
  inline typename result_of::Negative<Term>::type
  operator-(const Term& t)
  {
    typedef typename traits::to_expr<Term>::to Expr;
    return {Expr::get(t)};
  }


  // substract two terms
  // _____________________________________________________________________________
  template<typename Term1, typename Term2>
  inline typename result_of::Subtract<Term1, Term2>::type
  operator-(const Term1& t1, const Term2& t2)
  {
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    return {Expr1::get(t1), Expr2::get(t2)};
  }

} // end namespace AMDiS
