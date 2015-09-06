/** \file foreach.hpp */

#pragma once

#include <tuple>
#include <utility>
#include <traits/meta_basic.hpp>

namespace AMDiS
{
  namespace detail
  {
    template <class Tuple>
    constexpr int size(Tuple const&) { return std::tuple_size<Tuple>::value; }
    
    template <template<class...> class Tuple, class... Ts>
    constexpr int size(Tuple<Ts...> const&) { return sizeof...(Ts); }
    
    
    template <class Tuple, class F, int N>
    inline void for_each(Tuple, F, int_<N>, int_<N>) {}
    
    template <class Tuple, class F, int I, int N>
    inline void for_each(Tuple&& t, F&& f, int_<I>, int_<N>)
    {
      f(std::get<I>(std::forward<Tuple>(t)));
      for_each(std::forward<Tuple>(t), std::forward<F>(f), int_<I+1>(), int_<N>());
    }
    
  } // end namespace detail
  
  template <class Tuple, class F>
  inline void for_each(Tuple&& t, F&& f)
  {
    static constexpr int N = detail::size(t);
    detail::for_each(std::forward<Tuple>(t), std::forward<F>(f), int_<0>(), int_<N>());
  }

} // end namespace AMDiS
