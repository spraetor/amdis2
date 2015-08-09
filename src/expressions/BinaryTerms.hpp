/** \file BinaryTerms.hpp */

#pragma once

#include <string>

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>
#include <operations/functors.hpp>

#include "BaseTerms.hpp" // for ShapedExpr

namespace AMDiS 
{
  /// Expression with two arguments
  template <class Term1, class Term2, class Functor>
  struct BinaryTerm
      : public ShapedTerm_t<Term1, BinaryTerm<Term1, Term2, Functor> >,
        public LazyOperatorTerms<Term1, Term2>
  {
    using Self       = BinaryTerm;
    using Super      = LazyOperatorTerms<Term1, Term2>;    
    using value_type = typename std::result_of<Functor(Value_t<Term1>, Value_t<Term2>)>::type;
    
    /// constructor takes two terms
    template <class Term1_, class Term2_>
    constexpr BinaryTerm(Term1_&& A, Term2_&& B) 
      : Super(std::forward<Term1_>(A), std::forward<Term2_>(B)),
        term1(A), term2(A) 
    {}
    
    /// constructor takes two terms and the functor
    template <class Term1_, class Term2_, class Functor_>
    constexpr BinaryTerm(Term1_&& A, Term2_&& B, Functor_&& f) 
      : Super(std::forward<Term1_>(A), std::forward<Term2_>(B)),
        term1(A), term2(B),
        fct{f}
    {}
     
    int getDegree() const
    {
      return fct.getDegree( term1.getDegree(), term2.getDegree() );
    }
    
    /// access the elements of an term.
    value_type operator()(int i) const
    { 
      return fct( term1(i), term2(i) );
    }
    
    std::string str() const 
    { 
      return to_string(Functor()) + "(" + term1.str() + ", " + term2.str() + ")"; 
    }
    
    Term1 const& getTerm1() const { return term1; }
    Term2 const& getTerm2() const { return term2; }
    
  private:
    Term1 const& term1;
    Term2 const& term2;
    
    Functor fct;
  };
  
  
  /// Size of ElementwiseBinaryExpr
  template <class Term1, class Term2, class F>
  size_t size(BinaryTerm<Term1, Term2, F> const& term)
  {
    return size(term.getTerm1());
  }
  
  /// number of rows of ElementwiseBinaryExpr
  template <class Term1, class Term2, class F>
  size_t num_rows(BinaryTerm<Term1, Term2, F> const& term)
  {
    return num_rows(term.getTerm1());
  }
  
  /// number of columns of ElementwiseBinaryExpr
  template <class Term1, class Term2, class F>
  size_t num_cols(BinaryTerm<Term1, Term2, F> const& term)
  {
    return num_cols(term.getTerm1());
  }
  
  namespace traits 
  {
    /// \cond HIDDEN_SYMBOLS
    template <class Term1, class Term2, class F>
    struct category< BinaryTerm<Term1,Term2,F> > 
    {
      using tag        = typename category<Term1>::tag;
      using value_type = Result_t<F>;
      using size_type  = int; //max_size_type<Term1,Term2>;
    };
    /// \endcond
  }

  // Term1 + Term2
  template <class Term1, class Term2>
  using PlusTerm =
    BinaryTerm<Term1, Term2, 
      functors::plus<Value_t<Term1>, Value_t<Term2> > >;
      
  // Term1 - Term2
  template <class Term1, class Term2>
  using MinusTerm =
    BinaryTerm<Term1, Term2, 
      functors::minus<Value_t<Term1>, Value_t<Term2> > >;

  // Term1 * Term2
  template <class Term1, class Term2>
  using MultipliesTerm =
    BinaryTerm<Term1, Term2, 
      functors::multiplies<Value_t<Term1>, Value_t<Term2> > >;
      
  // Term1 / Term2
  template <class Term1, class Term2>
  using DividesTerm =
    BinaryTerm<Term1, Term2, 
      functors::divides<Value_t<Term1>, Value_t<Term2> > >;
  
} // end namespace AMDiS
