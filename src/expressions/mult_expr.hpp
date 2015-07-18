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



/** \file mult_expr.hpp */

#ifndef AMDIS_MULT_EXPRESSION_HPP
#define AMDIS_MULT_EXPRESSION_HPP

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"

namespace AMDiS 
{  
  namespace expressions 
  {
    
    /// Expressions that represents the multiplication of two expressions: E1 * E2
    template<typename Term1, typename Term2>
    struct Mult : public LazyOperatorTerm2<Term1, Term2>
    {
      typedef LazyOperatorTerm2<Term1, Term2> super;
      typedef typename traits::mult_type
      <
	typename Term1::value_type, 
	typename Term2::value_type
      >::type value_type;
      
      BOOST_STATIC_ASSERT_MSG( !(boost::is_same<value_type, traits::no_valid_type>::value), "********** ERROR: Can not multiply terms **********" );
      
      Mult(const Term1& term1_, const Term2& term2_)
	: super(term1_, term2_) {}
      
      int getDegree() const
      {
	return super::term1.getDegree() + super::term2.getDegree();
      }

      inline value_type operator()(const int& iq) const { return super::term1(iq) * super::term2(iq); }
      
      std::string str() const { return std::string("(") + super::term1.str() + " * " + super::term2.str() + ")"; }
    };
    
    
    /// Expressions that represents the inverse of an expressions: 1/E
    template<typename Term>
    struct MultInverse : public LazyOperatorTerm1<Term>
    {
      typedef LazyOperatorTerm1<Term> super;
      typedef typename Term::value_type value_type;
      
      MultInverse(const Term& term_)
	: super(term_) {}
      
      int getDegree() const
      {
	return super::term.getDegree();
      }

      // works only for scalar types
      // TODO: extend implementation to inverse of matrices
      inline value_type operator()(const int& iq) const { return 1.0 / super::term(iq); } 
      
      std::string str() const { return std::string("(") + super::term.str() + ")^(-1)"; }
    };
    
  } // end namespace expressions
  

  namespace result_of
  {
    template<typename Term1, typename Term2>
    struct Mult : boost::enable_if
      < 
	typename traits::is_valid_arg2<Term1, Term2>::type,
	expressions::Mult
	<
	  typename traits::to_expr<Term1>::type, 
	  typename traits::to_expr<Term2>::type
	>
      > {};
      
      
    template<typename Term1, typename Term2>
    struct Divide : boost::enable_if
      < 
	typename traits::is_valid_arg2<Term1, Term2>::type,
	expressions::Mult
	<
	  typename traits::to_expr<Term1>::type, 
	  expressions::MultInverse< typename traits::to_expr<Term2>::type >
	>
      > {};
      
  } // end namespace result_of
  

  // multiply two terms
  // _____________________________________________________________________________
  template<typename Term1, typename Term2>
  inline typename result_of::Mult<Term1, Term2>::type
  operator*(const Term1& t1, const Term2& t2)
  {
    typedef typename traits::to_expr<Term1>::to Expr1;
    typedef typename traits::to_expr<Term2>::to Expr2;
    return expressions::Mult< typename Expr1::type, typename Expr2::type >
	    (Expr1::get(t1), Expr2::get(t2));
  }


  // divide two terms
  // _____________________________________________________________________________
  template<typename Term1, typename Term2>
  inline typename result_of::Divide<Term1, Term2>::type
  operator/(const Term1& t1, const Term2& t2)
  {
    typedef typename traits::to_expr<Term2>::to Expr2;
    return t1 * expressions::MultInverse< typename Expr2::type >(Expr2::get(t2));
  }

} // end namespace AMDiS

#endif // AMDIS_MULT_EXPRESSION_HPP
