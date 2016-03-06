/** \file TermConcepts.hpp */

#pragma once

// std c++ headers
#include <utility>
#include <type_traits>

// AMDiS headers
#include "MatrixVector_fwd.hpp" // WorldVector<T>
#include "expressions/LazyOperatorTermBase.hpp"
#include "traits/basic.hpp"
#include "traits/concepts_base.hpp"
#include "traits/meta_basic.hpp"

namespace AMDiS
{
  namespace concepts
  {
    HAS_MEMBER_GENERATE(getDegree)

    // expressions concepts
    template <class... Ts>
    struct Term : and_<Term<Ts>...> {};

    template <class T>
    struct Term<T>
    {
    private:
      using Type = Decay_t<T>;

      static constexpr bool value0 = std::is_base_of<LazyOperatorTermBase, Type>::value;
      static constexpr bool value2 = traits::HasValueType<Type>::value;

      template <bool, class T1> struct check1 : traits::IsFunctor<T1, Value_t<T1>(WorldVector<double>)> {};
      template <      class T1> struct check1<false, T1> : false_ {};

    public:
      static constexpr bool value = value0 && value2 && check1<value2, Type>::value;
    };


    // check wether a functor has a member getDegree, with N arguments of type int
    template <class F, int N>
    struct TermFunctor
    {
      template <int I, class... Ints>
      static constexpr bool has_member(int_<I>, Ints...)
      {
        return has_member(int_<I-1>(), int{0}, Ints(0)...);
      }

      template <class... Ints>
      static constexpr bool has_member(int_<0>, Ints...)
      {
        return HAS_MEMBER(getDegree)<F, int(Ints...)>::value;
      }

      static constexpr bool value = has_member(int_<N>());
    };


    template <class F>
    using CoordsFunctor = and_<traits::IsFunctor<F, double(WorldVector<double>)>, not_<Term<F>>>;

  } // end namespace concepts


  namespace requires
  {
    // test whether one of the arguments is term
    template <class Result, class... Args>
    using Term = Requires_t<or_<concepts::Term<Args>...>, Result>;

    template <class F, int N>
    using TermFunctor = Requires_t<concepts::TermFunctor<F,N>>;

    template <class Result, class F>
    using CoordsFunctor = Requires_t<concepts::CoordsFunctor<F>, Result>;

  } // end namespace require
} // end namespace AMDiS
