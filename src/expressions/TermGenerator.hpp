/** \file TermGenerator.hpp */

#pragma once

// AMDiS includes
#include "expressions/BaseTerms.hpp"
#include "expressions/ConstantTerms.hpp"
#include "expressions/CoordsTerm.hpp"
#include "expressions/FunctorTerm.hpp"
#include "expressions/TermConcepts.hpp"

#include "traits/basic.hpp"

namespace AMDiS
{
  namespace detail
  {
    /// Generator class that takes a type as argument
    /// and returns the corresponding term-type.
    /// By default all arguments are interpreted as constants.
    template <class T, class = void>
    struct ToTerm
    {
      using type = RTConstant<T>;

      template <class T_>
      static constexpr type get(T_&& t)
      {
        return {std::forward<T_>(t)};
      }
    };

    // if the argument is already a term return it directly.
    template <class T>
    struct ToTerm<T, Requires_t<concepts::Term<T>>>
    {
      using type = T;

      template <class T_>
      static constexpr T_&& get(T_&& t)
      {
        return std::forward<T_>(t);
      }
    };

    // if the argument is a functor, that takes a WorldVector
    // and returns a double, see concept \ref CoordsFunctor.
    template <class F>
    struct ToTerm<F, Requires_t<concepts::CoordsFunctor<F>>>
    {
      using type = FunctorTerm<F, CoordsTerm<>>;

      template <class F_>
      static type get(F_&& f)
      {
        return {std::forward<F_>(f), CoordsTerm<>()};
      }
    };

  } // end namespace detail

  template <class T>
  using ToTerm = detail::ToTerm<Decay_t<T>>;

  template <class T>
  using ToTerm_t = typename ToTerm<T>::type;

  template <class T>
  inline auto toTerm(T&& t) RETURNS( ToTerm<T>::get(std::forward<T>(t)) )

} // end namespace AMDiS
