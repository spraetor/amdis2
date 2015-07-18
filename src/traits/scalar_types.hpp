/******************************************************************************
 *
 * AMDiS - Adaptive multidimensional simulations
 *
 * Copyright (C) 2013 Dresden University of Technology. All Rights Reserved.
 * Web: https://fusionforge.zih.tu-dresden.de/projects/amdis
 *
 * Authors: 
 * Simon Vey, Thomas Witkowski, Andreas Naumann, Simon Praetorius, et al.
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 *
 * This file is part of AMDiS
 *
 * See also license.opensource.txt in the distribution.
 * 
 ******************************************************************************/



/** \file scalar_types.hpp */

#ifndef AMDIS_TYPE_TRAITS_SCALAR_TYPES_HPP
#define AMDIS_TYPE_TRAITS_SCALAR_TYPES_HPP

#ifdef HAS_CPP11
#include <type_traits>
#endif

namespace AMDiS 
{
  namespace traits 
  {
	
    // type-traits
    // _________________________________________________________________________
    template<typename T>
    struct is_integer : boost::mpl::or_
      <
	typename boost::is_signed<T>::type,
	typename boost::is_unsigned<T>::type
      >::type {};
    
    template<typename T>
    struct is_numeric : boost::mpl::or_
      <
	typename boost::is_floating_point<T>::type,
	typename is_integer<T>::type
      >::type {};
    
  } // end namespace traits
  
} // end namespace AMDiS

#endif // AMDIS_TYPE_TRAITS_TYPES_HPP
