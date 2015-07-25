/** \file mult_expr.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "LazyOperatorTerm.h"

namespace AMDiS 
{  
  namespace expressions 
  {
    
    /// Expressions that represents the multiplication of two expressions: E1 * E2
    template <class Term1, class Term2>
    struct Mult : public LazyOperatorTerms<Term1, Term2>
    {
      using Super = LazyOperatorTerms<Term1, Term2>;
      using value_type = typename traits::mult_type< Value_t<Term1>, Value_t<Term2> >::type;
      
      static_assert( !(std::is_same<value_type, traits::no_valid_type>::value), 
                     "********** ERROR: Can not multiply terms **********" );
      
      Mult(Term1&& term1_, Term2&& term2_)
        : Super(std::forward<Term1>(term1_), std::forward<Term2>(term2_)) {}
      
      int getDegree() const
      {
        return degree<0>(*this) + degree<1>(*this);
      }

      inline value_type operator()(const int& iq) const { return term<0>(*this)(iq) * term<1>(*this)(iq); }
      
      std::string str() const { return std::string("(") + term<0>(*this).str() + " * " + term<1>(*this).str() + ")"; }
    };
    
    
    /// Expressions that represents the inverse of an expressions: 1/E
    template <class Term>
    struct MultInverse : public LazyOperatorTerms<Term>
    {
      using Super = LazyOperatorTerms<Term>;
      using value_type = Value_t<Term>;
      
      MultInverse(Term&& term_)
        : Super(std::forward<Term>(term_)) {}
      
      int getDegree() const
      {
        return degree<0>(*this);
      }

      // works only for scalar types
      // TODO: extend implementation to inverse of matrices
      inline value_type operator()(const int& iq) const { return 1.0 / term<0>(*this)(iq); } 
      
      std::string str() const { return std::string("(") + term<0>(*this).str() + ")^(-1)"; }
    };
    
  } // end namespace expressions
  

  namespace result_of
  {
    template <class Term1, class Term2>
    using Mult = enable_if
      < 
      	traits::is_valid_arg2<Term1, Term2>,
      	expressions::Mult
      	<
      	  typename traits::to_expr<Term1>::type, 
      	  typename traits::to_expr<Term2>::type
      	>
      >;
      
      
    template <class Term1, class Term2>
    using Divide = enable_if
      < 
      	traits::is_valid_arg2<Term1, Term2>,
      	expressions::Mult
      	<
      	  typename traits::to_expr<Term1>::type, 
      	  expressions::MultInverse< typename traits::to_expr<Term2>::type >
      	>
      >;
      
  } // end namespace result_of
  

  // multiply two terms
  // _____________________________________________________________________________
  template <class Term1, class Term2>
  inline typename result_of::Mult<Term1, Term2>::type
  operator*(const Term1& t1, const Term2& t2)
  {
    using Expr1 = traits::to_expr<Term1>;
    using Expr2 = traits::to_expr<Term2>;
    return {Expr1::get(t1), Expr2::get(t2)};
  }


  // divide two terms
  // _____________________________________________________________________________
  template <class Term1, class Term2>
  inline typename result_of::Divide<Term1, Term2>::type
  operator/(const Term1& t1, const Term2& t2)
  {
    using Expr2 = traits::to_expr<Term2>;
    return t1 * expressions::MultInverse< typename Expr2::type >(Expr2::get(t2));
  }

} // end namespace AMDiS
