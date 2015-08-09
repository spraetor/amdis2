/** \file UnaryTerms.hpp */

#pragma once

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>

#include "BaseTerms.hpp" // for ShapedExpr

namespace AMDiS 
{
  /// Expression with one argument
  template <class Term, class Functor>
  struct UnaryTerm
      : public ShapedTerm_t<Term, UnaryTerm<Term, Functor> >,
        public LazyOperatorTerms<Term>
  {
    using Self       = UnaryTerm;
    using Super      = LazyOperatorTerms<Term>;    
    using value_type = typename std::result_of<Functor(Value_t<Term>)>::type;
    
    /// constructor takes two termessions
    template <class Term_>
    constexpr UnaryTerm(Term_&& A) 
      : Super(std::forward<Term_>(A)),
        term(A)
    {}
    
    /// constructor takes two termessions
    template <class Term_, class F_>
    constexpr UnaryTerm(Term_&& A, F_&& f) 
      : Super(std::forward<Term_>(A)),
        term(A),
        fct{f}
    {}
     
    constexpr int getDegree() const
    {
      return fct.getDegree( term.getDegree() );
    }
    
    /// access the elements of an term.
    value_type operator()(int i) const
    { 
      return fct( term(i) );
    }
    
    std::string str() const 
    { 
      return to_string(Functor()) + "(" + term.str() + ")"; 
    }
    
    Term const& getTerm() const { return term; }
    
  private:
    Term const& term;
    
    Functor fct;
  };
  
  
  /// Size of ElementwiseBinaryExpr
  template <class Term, class F>
  size_t size(UnaryTerm<Term, F> const& term)
  {
    return size(term.getTerm());
  }
  
  /// number of rows of ElementwiseBinaryExpr
  template <class Term, class F>
  size_t num_rows(UnaryTerm<Term, F> const& term)
  {
    return num_rows(term.getTerm());
  }
  
  /// number of columns of ElementwiseBinaryExpr
  template <class Term, class F>
  size_t num_cols(UnaryTerm<Term, F> const& term)
  {
    return num_cols(term.getTerm());
  }
  
  namespace traits 
  {
    /// \cond HIDDEN_SYMBOLS
    template <class Term, class F>
    struct category< UnaryTerm<Term,F> > 
    {
      using tag        = typename category<Term>::tag;
      using value_type = Result_t<F>;
      using size_type  = int; //Size_t<Term>;
    };
    /// \endcond
    
  } // end namespace traits
  
  template <class Term>
  using NegateTerm = UnaryTerm<Term, functors::negate<Value_t<Term> > >;

} // end namespace AMDiS
