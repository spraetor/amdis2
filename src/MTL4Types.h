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


/** \file Mapper.h */

#ifndef MTL4TYPES_H
#define MTL4TYPES_H
#include <boost/numeric/mtl/mtl.hpp>

namespace AMDiS {

  namespace MTLTypes {
    typedef double value_type;
    typedef unsigned int size_type;
    typedef mtl::matrix::parameters<mtl::row_major, mtl::index::c_index, mtl::non_fixed::dimensions, false, size_type> para;
    typedef mtl::matrix::compressed2D< value_type, para > MTLMatrix;
    typedef mtl::vector::dense_vector< value_type > MTLVector;
    
#if defined(HAVE_PARALLEL_MTL4)
    typedef mtl::matrix::distributed< MTLMatrix > PMTLMatrix;
    typedef mtl::vector::distributed< MTLVector > PMTLVector;
#endif
  }

}

#endif

