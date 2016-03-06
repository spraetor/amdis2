#pragma once

// std c++ headers
#include <tuple>
#include <utility>

// AMDiS includes
#include "traits/meta_basic.hpp"
#include "traits/size.hpp"

namespace AMDiS
{
  namespace detail
  {    
    template <class F, class Tuple, int N>
    inline void for_each(F&&, Tuple&&, empty_range_<N>) {}
    
    template <class F, class Tuple, int I, int N>
    inline void for_each(F&& f, Tuple&& t, range_<I,N>)
    {
      (std::forward<F>(f))(std::get<I>(std::forward<Tuple>(t)));
      for_each(std::forward<F>(f), std::forward<Tuple>(t), range_<I+1,N>{});
    }
    
  } // end namespace detail
  
  template <class F, class Tuple, int N = Size<Decay_t<Tuple>>::value>
  inline void for_each(F&& f, Tuple&& t)
  {
    detail::for_each(std::forward<F>(f), std::forward<Tuple>(t), range_<0,N>{});
  }

} // end namespace AMDiS
