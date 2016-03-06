/** \file TermGenerator.hpp */

#pragma once

// AMDiS includes
#include "expressions/BaseTerms.hpp"
#include "expressions/ConstantTerms.hpp"
#include "expressions/CoordsTerm.hpp"
#include "expressions/FunctorTerm.hpp"
#include "expressions/TermConcepts.hpp"

#include "traits/basic.hpp"
#include "traits/meta_basic.hpp"

namespace AMDiS
{
  namespace detail
  {
    // named categories
    struct _term {};
    struct _functor {};
    struct _default {};

    template <class T>
    using TermCategory =
        if_then_else< concepts::Term<T>::value,            _term,
        if_then_else< concepts::CoordsFunctor<T>::value,   _functor,
                                                           _default > >;


    /// Generator class that takes a type as argument
    /// and returns the corresponding term-type.
    /// By default all arguments are interpreted as constants.
    template <class T, class Category> struct ToTerm;

    template <class T>
    struct ToTerm<T, _default>
    {
      using type = RTConstant<T>;

      static constexpr type get(T const& t)
      {
        return {t};
      }
    };

    // if the argument is already a term return it directly.
    template <class T>
    struct ToTerm<T, _term>
    {
      using type = T const&;

      static constexpr type get(T const& t)
      {
        return t;
      }
    };

    // if the argument is a functor, that takes a WorldVector
    // and returns a double, see concept \ref CoordsFunctor.
    template <class F>
    struct ToTerm<F, _functor>
    {
      using type = FunctorTerm<F, CoordsTerm<>>;

      static type get(F const& f)
      {
        return {f, CoordsTerm<>()};
      }
    };

  } // end namespace detail

  template <class T>
  using ToTerm = detail::ToTerm<Decay_t<T>, detail::TermCategory<Decay_t<T>> >;

  template <class T>
  using ToTerm_t = Decay_t<typename ToTerm<T>::type>;

  template <class T>
  inline typename ToTerm<T>::type toTerm(T const& t) { return ToTerm<T>::get(t); }

} // end namespace AMDiS
