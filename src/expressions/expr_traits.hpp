/** \file expr_traits.hpp */

#pragma once

#include "Traits.h"
#include "LazyOperatorTerm.h"
#include "value_expr.hpp"

namespace AMDiS
{
  namespace traits
  {

    // type-trait for short-cuts for constants
    // ___________________________________________________________________________

    template <class T>
    struct is_constant : std::is_arithmetic<T>::type {};

    template <class T>
    struct is_constant<WorldVector<T>> : bool_<true> {};

    template <class T>
    struct is_constant<WorldMatrix<T>> : bool_<true> {};


    // type-traits for terms
    // ___________________________________________________________________________
    template <class T>
    using is_expr = typename std::is_base_of<LazyOperatorTermBase, T>::type;


    // type-traits for arguments, filter terms and constants
    // ___________________________________________________________________________
    template <class T, class S = typename std::decay<T>::type>
    using is_valid_arg = or_<is_expr<S>,
          is_constant<S>>;

    template <class... Args>
    using is_valid_args = and_<and_<is_valid_arg<Args>...>,
          or_<is_expr<Args>...>>; // one of the args must be an expression

    template <class T1, class T2>
    using is_valid_arg2 = is_valid_args<T1, T2>;

    template <class T1, class T2, class T3>
    using is_valid_arg3 = is_valid_args<T1, T2, T3>;

    // expressions
    //     template <class T>
    //     struct category<T, typename boost::enable_if<is_expr<T> >::type >
    //     {
    //       using tag = tag::expression;
    //       using value_type = Value_t<T>;
    //       // typedef size_t               size_type;
    //     };


    // type-conversion
    // ___________________________________________________________________________

    template <class T, bool IsValid, bool IsConstant>
    struct to_expr_aux
    {
      using to = to_expr_aux;
      using type = void;
    };

    template <class T>
    struct to_expr_aux<T, /*IsValid: */true, /*IsConstant: */false>
    {
      using to = to_expr_aux;
      using type = T;

      static constexpr const type& get(const T& t)
      {
        return t;
      }
    };

    template <class T> // T is numeric constant
    struct to_expr_aux<T, /*IsValid: */true, /*IsConstant: */true>
    {
      using to = to_expr_aux;
      using type = ::AMDiS::expressions::RValue<T>;

      template <class T_>
      static constexpr type get(T_&& t)
      {
        return {t};
      }
    };

    template <class T>
    using to_expr = to_expr_aux<T, is_valid_arg<T>::value, is_constant<T>::value>;

  } // end namespace traits

} // end namespace AMDiS
