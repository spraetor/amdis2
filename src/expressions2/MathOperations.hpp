/** \file MatrixVectorOperations.hpp */

#pragma once

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
  RTConstant<T>
  inline constant(T&& value)
  {
    return {std::forward<T>(value)};
  }

  template <int I>
  CTConstant<I>
  inline constant()
  {
    return {};
  }

  template <class T>
  Reference<T>
  inline ref(T const& value)
  {
    return {value};
  }


  /// return coords vector
  inline CoordsTerm<-1> X()
  {
    return {};
  }

  /// return component of coords vector, given as function argument
  inline CoordsTerm<-2> X(int comp)
  {
    return {comp};
  }

  /// return component of coords vector, given as template argument
  template <int Comp>
  inline CoordsTerm<Comp> X()
  {
    return {};
  }

  // ---------------------------------------------------------------------------

  template <class Name = _unknown, class T>
  ValueOf<DOFVector<T>, Name>
  inline valueOf(DOFVector<T> const& vector)
  {
    return {vector};
  }

  template <class Name = _unknown, class T>
  ValueOf<DOFVector<T>, Name>
  inline valueOf(DOFVector<T> const* vector)
  {
    return {vector};
  }

  template <class Name = _unknown, class T>
  GradientOf<DOFVector<T>, Name>
  inline gradientOf(DOFVector<T> const& vector)
  {
    return {vector};
  }

  template <class Name = _unknown, class T>
  GradientOf<DOFVector<T>, Name>
  inline gradientOf(DOFVector<T> const* vector)
  {
    return {vector};
  }

  // ---------------------------------------------------------------------------

  template <class F, class... Terms>
  requires::Term<FunctorTerm<F, ToTerm_t<Terms>...>, Terms...>
  inline func(F&& f, Terms&& ... ts)
  {
    return {std::forward<F>(f), toTerm(std::forward<Terms>(ts))...};
  }

  // ---------------------------------------------------------------------------
  // Operations (elementwise)


  template <class Term1, class Term2>
  using PlusTerm =
    FunctorTerm<functors::plus<Value_t<Term1>, Value_t<Term2>>, Term1, Term2>;

  /// expression for V + W
  template <class T1, class T2>
  requires::Term<PlusTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>
      inline operator+(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }

  // ------------

  template <class Term1, class Term2>
  using MinusTerm =
    FunctorTerm<functors::minus<Value_t<Term1>, Value_t<Term2>>, Term1, Term2>;

  /// expression for V - W
  template <class T1, class T2>
  requires::Term<MinusTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>
      inline operator-(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }

  // ------------

  /// expression for -V
  template <class T>
  FunctorTerm<functors::negate<Value_t<T>>, T>
                                        inline operator-(BaseTerm<T> const& term)
  {
    return {term.sub()};
  }

  // ------------

  template <class Term1, class Term2>
  using MultipliesTerm =
    FunctorTerm<functors::multiplies<Value_t<Term1>, Value_t<Term2>>, Term1, Term2>;

  /// expression for V .* W
  template <class T1, class T2>
  requires::Term<MultipliesTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>
      inline operator*(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }

  // ------------

  template <class Term1, class Term2>
  using DividesTerm =
    FunctorTerm<functors::divides<Value_t<Term1>, Value_t<Term2>>, Term1, Term2>;

  /// expression for V ./ W
  template <class T1, class T2>
  requires::Term<DividesTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>
      inline operator/(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }


  // ---------------------------------------------------------------------------


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


  // ---------------------------------------------------------------------------


  template <class Term1, class Term2>
  using MinTerm =
    FunctorTerm<functors::min<Common_t<Value_t<Term1>, Value_t<Term2>>>, Term1, Term2>;

  /// expression for min(V, W)
  template <class T1, class T2>
  requires::Term<MinTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>
      inline min(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }

  // ------------

  template <class Term1, class Term2>
  using MaxTerm =
    FunctorTerm<functors::max<Common_t<Value_t<Term1>, Value_t<Term2>>>, Term1, Term2>;

  /// expression for max(V, W)
  template <class T1, class T2>
  requires::Term<MaxTerm<ToTerm_t<T1>, ToTerm_t<T2>>, T1, T2>
      inline max(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
  }


  // ---------------------------------------------------------------------------


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


  // ---------------------------------------------------------------------------


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


  // ---------------------------------------------------------------------------


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
      inline atan2(T1&& t1, T2&& t2)
  {
    return {toTerm(std::forward<T1>(t1)), toTerm(std::forward<T2>(t2))};
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

} // end namespace AMDiS
