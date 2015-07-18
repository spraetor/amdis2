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



/** \file to_string.hpp */

#ifndef AMDIS_UTILITY_TO_STRING_HPP
#define AMDIS_UTILITY_TO_STRING_HPP

#include <string>
#include <boost/lexical_cast.hpp>

namespace AMDiS {

  template <class T>
  inline std::string to_string(const T& value)
  {
    return boost::lexical_cast<std::string>(value);
  }

} // end namspace AMDiS

#endif // AMDIS_UTILITY_TO_STRING_HPP
