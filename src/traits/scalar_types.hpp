/** \file scalar_types.hpp */

#pragma once

#include <type_traits>

namespace AMDiS 
{
  namespace traits 
  {
    // type-traits
    // _________________________________________________________________________
    template <class T>
    using is_integral = std::is_integral<typename std::decay<T>::type>;
    
    template <class T>
    using is_arithmetic = std::is_arithmetic<typename std::decay<T>::type>;
    
  } // end namespace traits
  
} // end namespace AMDiS
