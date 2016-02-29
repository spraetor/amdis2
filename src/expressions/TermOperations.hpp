#pragma once

// AMDiS includes
#include "expressions/FunctorTerm.hpp"
#include "expressions/TermConcepts.hpp"
#include "expressions/TermGenerator.hpp"

#include "operations/functors.hpp"


// Macro to generate unary operators on term-expressions
#ifndef AMDIS_UNARY_TERM_OPERATION
#define AMDIS_UNARY_TERM_OPERATION( OP, FUNCTOR )                         \
  template <class T>                                                      \
  FunctorTerm<functors::FUNCTOR<Value_t<T>>, T>                           \
  inline OP (BaseTerm<T> const& term)                                     \
  {                                                                       \
    return {term.sub()};                                                  \
  }
#endif

// Macro to generate binary operators on term-expressions
#ifndef AMDIS_BINARY_TERM_OPERATION
#define AMDIS_BINARY_TERM_OPERATION( OP, FUNCTOR )                        \
  template <class Term1, class Term2>                                     \
  using FUNCTOR ## _term =                                                \
    FunctorTerm<functors::FUNCTOR<Value_t<Term1>, Value_t<Term2>>,        \
                Term1, Term2>;                                            \
  template <class T1, class T2>                                           \
    requires::Term<FUNCTOR ## _term<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>  \
  inline OP (T1 const& t1, T2 const& t2)                                  \
  {                                                                       \
    return {toTerm(t1), toTerm(t2)};                                      \
  }
#endif

namespace AMDiS
{
  // ______ generator functions for pointwise basic arithmetic operations _____

  AMDIS_BINARY_TERM_OPERATION( operator+ ,  plus )
  AMDIS_BINARY_TERM_OPERATION( operator- ,  minus )
  AMDIS_BINARY_TERM_OPERATION( operator* ,  multiplies )
  AMDIS_BINARY_TERM_OPERATION( operator/ ,  divides )

  AMDIS_UNARY_TERM_OPERATION( operator- ,  negate )
  

  // ______ generator functions for pointwise basic logical operations ________
  
  AMDIS_BINARY_TERM_OPERATION( operator== ,  equal )
  AMDIS_BINARY_TERM_OPERATION( operator!= ,  unequal )
  AMDIS_BINARY_TERM_OPERATION( operator<  ,  less )
  AMDIS_BINARY_TERM_OPERATION( operator>  ,  greater )

  AMDIS_BINARY_TERM_OPERATION( operator&& ,  logical_and )
  AMDIS_BINARY_TERM_OPERATION( operator|| ,  logical_or )

  
  // ______ generator functions for pointwise tertiary logical operations ________
  
  template <class Term0, class Term1, class Term2>
  using ConditionalTerm =
    FunctorTerm< functors::conditional<Value_t<Term1>, Value_t<Term2>>,
                 Term0, Term1, Term2>;
                 
  // if_(t0, t1, t2) = t0 ? t1 : t2
  template <class T0, class T1, class T2>
    requires::Term<ConditionalTerm<ToTerm_t<T0>, ToTerm_t<T1>, ToTerm_t<T2>>, T0, T1, T2>
  inline if_(T0 const& t0, T1 const& t1, T2 const& t2)
  {
    return {toTerm(t0), toTerm(t1), toTerm(t2)};
  }
  

} // end namespace AMDiS
