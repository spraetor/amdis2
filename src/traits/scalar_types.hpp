/** \file scalar_types.hpp */

#pragma once

// std c++ headers
#include <type_traits>

namespace AMDiS
{
  namespace traits
  {    
    template <class T>
    using IsIntegral = std::is_integral<typename std::decay<T>::type>;

    template <class T>
    using IsArithmetic = std::is_arithmetic<typename std::decay<T>::type>;

  } // end namespace traits

  namespace concepts
  {
    template <class T>
    using Integral = traits::IsIntegral<T>;
    
    template <class T>
    using Arithmetic = traits::IsArithmetic<T>;
    
  } // end namespace concepts
  
} // end namespace AMDiS
