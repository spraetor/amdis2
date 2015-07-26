/** \file mult_expr.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"

namespace AMDiS 
{  
  namespace expressions 
  {
    
    /// Expressions that represents the multiplication of two expressions: E1 * E2
    template <class Term0, class Term1>
    struct Mult : public LazyOperatorTerms<Term0, Term1>
    {
      using Super = LazyOperatorTerms<Term0, Term1>;
      using value_type = typename traits::mult_type< Value_t<Term0>, Value_t<Term1> >::type;
      
      static_assert( !(std::is_same<value_type, traits::no_valid_type>::value), 
                     "********** ERROR: Can not multiply terms **********" );
      
      template <class Term0_, class Term1_>
      Mult(Term0_&& term0_, Term1_&& term1_)
        : Super(std::forward<Term0_>(term0_), std::forward<Term1_>(term1_)) {}
      
      int getDegree() const
      {
        return Super::getDegree(_0) + Super::getDegree(_1);
      }

      inline value_type operator()(const int& iq) const { return Super::getTerm(_0)(iq) * Super::getTerm(_1)(iq); }
      
      std::string str() const { return std::string("(") + Super::getTerm(_0).str() + " * " + Super::getTerm(_1).str() + ")"; }
    };
    
    
    /// Expressions that represents the inverse of an expressions: 1/E
    template <class Term>
    struct MultInverse : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      
      template <class Term_>
      MultInverse(Term_&& term_)
        : Super(std::forward<Term_>(term_)) {}
      
      int getDegree() const
      {
        return Super::getDegree(_0);
      }

      // works only for scalar types
      // TODO: extend implementation to inverse of matrices
      inline value_type operator()(const int& iq) const { return 1.0 / Super::getTerm()(iq); } 
      
      std::string str() const { return std::string("(") + Super::getTerm().str() + ")^(-1)"; }
    };
    
  } // end namespace expressions
  

  namespace result_of
  {
    template <class Term0, class Term1>
    using Mult = enable_if
      < 
      	traits::is_valid_arg2<Term0, Term1>,
      	expressions::Mult
      	<
      	  typename traits::to_expr<Term0>::type, 
      	  typename traits::to_expr<Term1>::type
      	>
      >;
      
      
    template <class Term0, class Term1>
    using Divide = enable_if
      < 
      	traits::is_valid_arg2<Term0, Term1>,
      	expressions::Mult
      	<
      	  typename traits::to_expr<Term0>::type, 
      	  expressions::MultInverse< typename traits::to_expr<Term1>::type >
      	>
      >;
      
  } // end namespace result_of
  

  // multiply two terms
  // _____________________________________________________________________________
  template <class Term0, class Term1>
  inline typename result_of::Mult<Term0, Term1>::type
  operator*(const Term0& t0, const Term1& t1)
  {
    using Expr0 = traits::to_expr<Term0>;
    using Expr1 = traits::to_expr<Term1>;
    return {Expr0::get(t0), Expr1::get(t1)};
  }


  // divide two terms
  // _____________________________________________________________________________
  template <class Term0, class Term1>
  inline typename result_of::Divide<Term0, Term1>::type
  operator/(const Term0& t0, const Term1& t1)
  {
    using Expr1 = traits::to_expr<Term1>;
    return t0 * expressions::MultInverse< typename Expr1::type >(Expr1::get(t1));
  }

} // end namespace AMDiS
