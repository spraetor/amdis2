/** \file MatrixVectorOperations.hpp */

#pragma once

#include "UnaryTerms.hpp"
#include "BinaryTerms.hpp"
#include "FunctorTerm.hpp"
#include "ConstantTerms.hpp"
#include "CoordsTerm.hpp"
#include "DOFVectorTerms.hpp"
#include "TermConcepts.hpp"
#include "TermGenerator.hpp"

#include "cmath_functors.hpp"

namespace AMDiS 
{
  // ___________________________________________________________________________
  // generator functions for constants and references

  template <class T>
  constexpr RTConstant<T> constant(T&& value) { return {std::forward<T>(value)}; }
  
  template <int I>
  constexpr CTConstant<I> constant() { return {}; }
  
  template <class T>
  constexpr Reference<T> ref(T const& value) { return {value}; }
  
  
  /// return coords vector
  inline CoordsTerm<-1> X() { return {}; }
  
  /// return component of coords vector, given as function argument
  inline CoordsTerm<-2> X(int comp) { return {comp}; }
  
  /// return component of coords vector, given as template argument
  template <int Comp>
  inline CoordsTerm<Comp> X() { return {}; }
  
  // ---------------------------------------------------------------------------
  
  template <class Name = _unknown, class T>
  constexpr ValueOf<DOFVector<T>, Name> 
  valueOf(DOFVector<T> const& vector) { return {vector}; }

  template <class Name = _unknown, class T>
  constexpr ValueOf<DOFVector<T>, Name> 
  valueOf(DOFVector<T> const* vector) { return {vector}; }
  
  template <class Name = _unknown, class T>
  constexpr GradientOf<DOFVector<T>, Name> 
  gradientOf(DOFVector<T> const& vector) { return {vector}; }

  template <class Name = _unknown, class T>
  constexpr GradientOf<DOFVector<T>, Name> 
  gradientOf(DOFVector<T> const* vector) { return {vector}; }
  
  // ---------------------------------------------------------------------------
  
  template <class F, class... Terms>
  constexpr requires::Term< FunctorTerm< F, ToTerm_t<Terms>... >, Terms... >
  func(F&& f, Terms&&... ts)
  {
    return {std::forward<F>(f), toTerm(std::forward<Terms>(ts))...}; 
  }
  
  // ---------------------------------------------------------------------------
  // Operations (elementwise)
  
  /// expression for V + W
  template <class T1, class T2>
  constexpr requires::Term< PlusTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2 > 
  operator+(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
  
  /// expression for V - W
  template <class T1, class T2>
  constexpr requires::Term< MinusTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2 > 
  operator-(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
  
  /// expression for -V
  template <class T>
  constexpr NegateTerm<T>
  operator-(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  
  /// expression for V .* W
  template <class T1, class T2>
  constexpr requires::Term< MultipliesTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2 >  
  operator*(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
  
  /// expression for V ./ W
  template <class T1, class T2>
  constexpr requires::Term< DividesTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2 > 
  operator/(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
  
  // ---------------------------------------------------------------------------
  
  
  /// expression for pos<p>(V)
  template <int p, class T>
  UnaryTerm<T, functors::pow<p, Value_t<T>> >
  constexpr pow(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for sqr(V)
  template <class T>
  UnaryTerm<T, functors::pow<2, Value_t<T>> >
  constexpr sqr(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for pos<p>(V)
  template <int p, class T>
  UnaryTerm<T, functors::root<p, Value_t<T>> >
  constexpr root(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for sqrt(V)
  template <class T>
  UnaryTerm<T, functors::root<2, Value_t<T>> >
  constexpr sqrt(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  
  // ---------------------------------------------------------------------------
  
  
  /// expression for min(V, W)
  template <class T1, class T2, 
            class T1_ = ToTerm_t<T1>, class T2_ = ToTerm_t<T2>, 
            class V_  = Common_t<Value_t<T1_>, Value_t<T2_>> >
  requires::Term< BinaryTerm<T1_, T2_, functors::min<V_>>, T1, T2> 
  constexpr min(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
  /// expression for max(V, W)
  template <class T1, class T2, 
            class T1_ = ToTerm_t<T1>, class T2_ = ToTerm_t<T2>, 
            class V_  = Common_t<Value_t<T1_>, Value_t<T2_>> >
  requires::Term< BinaryTerm<T1_, T2_, functors::max<V_>>, T1, T2> 
  constexpr max(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
  
  // ---------------------------------------------------------------------------
  
  
  /// expression for abs(V)
  template <class T>
  UnaryTerm<T, functors::abs<Value_t<T>> >
  constexpr abs(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for signum(V)
  template <class T>
  UnaryTerm<T, functors::Signum<Value_t<T>> >
  constexpr signum(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for floor(V)
  template <class T>
  UnaryTerm<T, functors::Floor<Value_t<T>> >
  constexpr floor(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for ceil(V)
  template <class T>
  UnaryTerm<T, functors::Ceil<Value_t<T>> >
  constexpr ceil(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  
  // ---------------------------------------------------------------------------
  
  
  /// expression for abs(V)
  template <class T>
  UnaryTerm<T, functors::Exp<Value_t<T>> >
  constexpr exp(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for abs(V)
  template <class T>
  UnaryTerm<T, functors::Log<Value_t<T>> >
  constexpr log(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  
  // ---------------------------------------------------------------------------
  
  
  /// expression for sin(V)
  template <class T>
  UnaryTerm<T, functors::Sin<Value_t<T>> >
  constexpr sin(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for cos(V)
  template <class T>
  UnaryTerm<T, functors::Cos<Value_t<T>> >
  constexpr cos(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for tan(V)
  template <class T>
  UnaryTerm<T, functors::Tan<Value_t<T>> >
  constexpr tan(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  
  /// expression for asin(V)
  template <class T>
  UnaryTerm<T, functors::Asin<Value_t<T>> >
  constexpr asin(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for acos(V)
  template <class T>
  UnaryTerm<T, functors::Acos<Value_t<T>> >
  constexpr acos(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for atan(V)
  template <class T>
  UnaryTerm<T, functors::Atan<Value_t<T>> >
  constexpr atan(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for V ./ W
  template <class T1, class T2, 
            class T1_ = ToTerm_t<T1>, class T2_ = ToTerm_t<T2> >
  requires::Term< BinaryTerm<T1_, T2_, functors::Atan2<Value_t<T1_>, Value_t<T2_>>>, T1, T2> 
  constexpr atan2(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }
  
  
  /// expression for sin(V)
  template <class T>
  UnaryTerm<T, functors::Sinh<Value_t<T>> >
  constexpr sinh(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for cos(V)
  template <class T>
  UnaryTerm<T, functors::Cosh<Value_t<T>> >
  constexpr cosh(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for tan(V)
  template <class T>
  UnaryTerm<T, functors::Tanh<Value_t<T>> >
  constexpr tanh(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  
  /// expression for asin(V)
  template <class T>
  UnaryTerm<T, functors::Asinh<Value_t<T>> >
  constexpr asinh(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for acos(V)
  template <class T>
  UnaryTerm<T, functors::Acosh<Value_t<T>> >
  constexpr acosh(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
  /// expression for atan(V)
  template <class T>
  UnaryTerm<T, functors::Atanh<Value_t<T>> >
  constexpr atanh(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }
  
} // end namespace AMDiS
