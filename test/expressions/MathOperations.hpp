#pragma once

// AMDiS includes
#include "expressions/FunctorTerm.hpp"
#include "expressions/ConstantTerms.hpp"
#include "expressions/CoordsTerm.hpp"
#include "expressions/TermConcepts.hpp"
#include "expressions/TermGenerator.hpp"

#include "expressions/TermOperations.hpp"

#include "expressions/cmath_functors.hpp"


namespace AMDiS
{
  // ______ generator functions for constants and references __________________

  /// Returns an expression the evaluates to a constant value given as
  /// argument \param value.
  template <class T>
  inline RTConstant<T> constant(T const& value)
  {
    return {value};
  }

  /// Returns an expression the evaluates to a compile-time constant value,
  /// given as template argument.
  template <int I>
  inline CTConstant<I> constant()
  {
    return {};
  }


  /// Returns an expression the evaluates to a value given as argument
  /// \param value. The value is stored as reference.
  template <class T>
  inline Reference<T> ref(T const& value)
  {
    return {value};
  }

  /// return an expression the evaluates to a value given as argument \param value.
  /// The value is stored as reference.
  template <class T>
  inline Reference<T> var(T const& value)
  {
    return {value};
  }


  /// Returns an expression that evaluates to a coordinate vector.
  inline CoordsTerm<-1> X()
  {
    return {};
  }

  /// Returns an expression that evaluates to a component of coordinate vector,
  /// where the component index is given as function argument.
  inline CoordsTerm<-2> X(int comp)
  {
    return {comp};
  }

  /// Returns an expression that evaluates to a component of coordinate vector,
  /// where the component index is given as template argument.
  template <int Comp>
  inline CoordsTerm<Comp> X()
  {
    return {};
  }


  /// Returns an expression that evaluates a functor \param f at the evaluated terms \param ts.
  template <class F, class... Terms>
  requires::Term<FunctorTerm<Decay_t<F>, ToTerm_t<Decay_t<Terms>>...>, Decay_t<Terms>...>
  inline func(F&& f, Terms&& ... terms)
  {
    return {std::forward<F>(f), toTerm(std::forward<Terms>(terms))...};
  }

  // ______ generator functions for pointwise cmath operations __________________
#if 0

  /// expression for pos<p>(V)
  template <int p, class T>
  FunctorTerm<functors::pow<p, Value_t<T>>, T>
  inline pow(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for sqr(V)
  template <class T>
  FunctorTerm<functors::pow<2, Value_t<T>>, T>
  inline sqr(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for pos<p>(V)
  template <int p, class T>
  FunctorTerm<functors::root<p, Value_t<T>>, T>
  inline root(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for sqrt(V)
  template <class T>
  FunctorTerm<functors::root<2, Value_t<T>>, T>
  inline sqrt(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }


  // --------------

  template <class Term1, class Term2>
  using MinTerm =
    FunctorTerm<functors::min<Common_t<Value_t<Term1>, Value_t<Term2>>>, Term1, Term2>;

  /// expression for min(V, W)
  template <class T1, class T2>
  requires::Term<MinTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>
  inline min(T1 const& t1, T2 const& t2)
  {
    return {toTerm(t1), toTerm(t2)};
  }

  // ------------

  template <class Term1, class Term2>
  using MaxTerm =
    FunctorTerm<functors::max<Common_t<Value_t<Term1>, Value_t<Term2>>>, Term1, Term2>;

  /// expression for max(V, W)
  template <class T1, class T2>
  requires::Term<MaxTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>
  inline max(T1 const& t1, T2 const& t2)
  {
    return {toTerm(t1), toTerm(t2)};
  }


  // ------------


  /// expression for abs(V)
  template <class T>
  FunctorTerm<functors::abs<Value_t<T>>, T>
  inline abs(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for signum(V)
  template <class T>
  FunctorTerm<functors::Signum<Value_t<T>>, T>
  inline signum(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for floor(V)
  template <class T>
  FunctorTerm<functors::Floor<Value_t<T>>, T>
  inline floor(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for ceil(V)
  template <class T>
  FunctorTerm<functors::Ceil<Value_t<T>>, T>
  inline ceil(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }


  // ------------


  /// expression for abs(V)
  template <class T>
  FunctorTerm<functors::Exp<Value_t<T>>, T>
  inline exp(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for abs(V)
  template <class T>
  FunctorTerm<functors::Log<Value_t<T>>, T>
  inline log(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }


  // ------------


  /// expression for sin(V)
  template <class T>
  FunctorTerm<functors::Sin<Value_t<T>>, T>
  inline sin(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for cos(V)
  template <class T>
  FunctorTerm<functors::Cos<Value_t<T>>, T>
  inline cos(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for tan(V)
  template <class T>
  FunctorTerm<functors::Tan<Value_t<T>>, T>
  inline tan(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }


  /// expression for asin(V)
  template <class T>
  FunctorTerm<functors::Asin<Value_t<T>>, T>
  inline asin(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for acos(V)
  template <class T>
  FunctorTerm<functors::Acos<Value_t<T>>, T>
  inline acos(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for atan(V)
  template <class T>
  FunctorTerm<functors::Atan<Value_t<T>>, T>
  inline atan(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  // ------------

  template <class Term1, class Term2>
  using ATan2Term =
    FunctorTerm<functors::Atan2<Value_t<Term1>, Value_t<Term2>>, Term1, Term2>;

  /// expression for atan2(V, W)
  template <class T1, class T2>
  requires::Term<ATan2Term<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>
      inline atan2(T1 const& t1, T2 const& t2)
  {
    return {toTerm(t1), toTerm(t2)};
  }

  // ------------

  /// expression for sin(V)
  template <class T>
  FunctorTerm<functors::Sinh<Value_t<T>>, T>
  inline sinh(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for cos(V)
  template <class T>
  FunctorTerm<functors::Cosh<Value_t<T>>, T>
  inline cosh(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for tan(V)
  template <class T>
  FunctorTerm<functors::Tanh<Value_t<T>>, T>
  inline tanh(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }


  /// expression for asin(V)
  template <class T>
  FunctorTerm<functors::Asinh<Value_t<T>>, T>
  inline asinh(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for acos(V)
  template <class T>
  FunctorTerm<functors::Acosh<Value_t<T>>, T>
  inline acosh(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  /// expression for atan(V)
  template <class T>
  FunctorTerm<functors::Atanh<Value_t<T>>, T>
  inline atanh(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
#endif
} // end namespace AMDiS
