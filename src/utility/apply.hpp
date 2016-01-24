/** \file int_seq.hpp */

#pragma once

#include <tuple>
#include <array>

#include <traits/basic.hpp>
#include <traits/size.hpp>
#include "int_seq.hpp"

namespace AMDiS
{
  namespace detail
  {
    // return f(t[0], t[1], t[2], ...)
    template <class F, class Tuple, int... I>
    constexpr auto apply_impl(F&& f, Tuple&& t, AMDiS::Seq<I...>) RETURNS
    (
      (std::forward<F>(f))(std::get<I>(std::forward<Tuple>(t))...)
    )
    
  } // end namespace detail
  

  template <class F, class Tuple, int N = Size<Decay_t<Tuple>>::value>
  constexpr auto apply(F&& f, Tuple&& t) RETURNS
  (
    detail::apply_impl(std::forward<F>(f), std::forward<Tuple>(t), MakeSeq_t<N>{})
  )
  
} // end namespace AMDiS
