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



/** \file add_expr.hpp */

#ifndef AMDIS_ADD_EXPRESSION_HPP
#define AMDIS_ADD_EXPRESSION_HPP

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"

#include <boost/static_assert.hpp>

namespace AMDiS 
{
  namespace expressions 
  {
    /// Expression that represents "E1 + E2"
    template<typename Term1, typename Term2>
    struct Add : public LazyOperatorTerm2<Term1, Term2>
    {
      typedef LazyOperatorTerm2<Term1, Term2> super;
      typedef typename traits::add_type
      <
	typename Term1::value_type,
	typename Term2::value_type
      >::type value_type;
      
      BOOST_STATIC_ASSERT_MSG( !(boost::is_same<value_type, traits::no_valid_type>::value), "********** ERROR: Can not add terms **********" );
      
      Add(const Term1& term1_, const Term2& term2_)
	: super(term1_, term2_) {}
      
      int getDegree() const
      {
	return std::max(super::term1.getDegree(), super::term2.getDegree());
      }

      inline value_type operator()(const int& iq) const { return super::term1(iq) + super::term2(iq); }
      
      std::string str() const { return std::string("(") + super::term1.str() + " + " + super::term2.str() + ")"; }
    };
    
    /// Expression that represents "E1 - E2"
    template<typename Term1, typename Term2>
    struct Subtract : public LazyOperatorTerm2<Term1, Term2>
    {
      typedef LazyOperatorTerm2<Term1, Term2> super;
      typedef typename traits::add_type
      <
	typename Term1::value_type,
	typename Term2::value_type
      >::type value_type;
      
      BOOST_STATIC_ASSERT_MSG( !(boost::is_same<value_type, traits::no_valid_type>::value), "********** ERROR: Can not subtract terms **********" );
      
      Subtract(const Term1& term1_, const Term2& term2_)
	: super(term1_, term2_) {}
      
      int getDegree() const
      {
	return std::max(super::term1.getDegree(), super::term2.getDegree());
      }

      inline value_type operator()(const int& iq) const { return super::term1(iq) - super::term2(iq); }
      
      std::string str() const { return std::string("(") + super::term1.str() + " - " + super::term2.str() + ")"; }
    };
    
    /// Expression that represents "-E"
    template<typename Term>
    struct Negative : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      
      Negative(const Term& term_)
	: super(term_) {}
      
      int getDegree() const
      {
	return super::term.getDegree();
      }

      inline value_type operator()(const int& iq) const { return -super::term(iq); }
      
      std::string str() const { return std::string("-(") + super::term.str() + ")"; }
    };
    
  } // end namespace expressions


  namespace result_of
  {
    template<typename Term1, typename Term2>
    struct Add : boost::enable_if
      < 
	typename traits::is_valid_arg2<Term1, Term2>::type,
	expressions::Add
	<
	  typename traits::to_expr<Term1>::type, 
	  typename traits::to_expr<Term2>::type
	>
      > {};
      
      
    template<typename Term>
    struct Negative : boost::enable_if
      < 
	typename traits::is_valid_arg<Term>::type,
	expressions::Negative
	<
	  typename traits::to_expr<Term>::type
	>
      > {};
      
    template<typename Term1, typename Term2>
    struct Subtract : boost::enable_if
      < 
	typename traits::is_valid_arg2<Term1, Term2>::type,
	expressions::Subtract
	<
	  typename traits::to_expr<Term1>::type, 
	  typename traits::to_expr<Term2>::type
	>
      > {};
  }


  // add two terms
  // _____________________________________________________________________________
  template<typename Term1, typename Term2>
  inline typename result_of::Add<Term1, Term2>::type
  operator+(const Term1& t1, const Term2& t2)
  {
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    return expressions::Add< typename Expr1::type, typename Expr2::type >
	    (Expr1::get(t1), Expr2::get(t2));
  }


  // negative of a term
  // _____________________________________________________________________________
  template<typename Term>
  inline typename result_of::Negative<Term>::type
  operator-(const Term& t)
  {
    typedef typename traits::to_expr<Term>::to Expr;
    return expressions::Negative< typename Expr::type >(Expr::get(t));
  }


  // substract two terms
  // _____________________________________________________________________________
  template<typename Term1, typename Term2>
  inline typename result_of::Subtract<Term1, Term2>::type
  operator-(const Term1& t1, const Term2& t2)
  {
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    return expressions::Subtract< typename Expr1::type, typename Expr2::type >
	    (Expr1::get(t1), Expr2::get(t2));
  }

} // end namespace AMDiS

#endif // AMDIS_ADD_EXPRESSION_HPP
