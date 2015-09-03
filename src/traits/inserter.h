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



/** \file inserter.hpp */

#ifndef AMDIS_COLLECTION_INSERTER_HPP
#define AMDIS_COLLECTION_INSERTER_HPP

#include "MTL4Types.h"

namespace AMDiS
{
  template<typename T>
  struct Collection {};

  template<>
  struct Collection<MTLTypes::MTLMatrix>
  {
    typedef mtl::matrix::inserter<MTLTypes::MTLMatrix> Inserter;
  };

  template<>
  struct Collection<MTLTypes::MTLVector>
  {
    typedef mtl::vector::inserter<MTLTypes::MTLVector> Inserter;
    typedef MTLTypes::MTLMatrix PreconditionMatrix;
  };
}

#endif // AMDIS_COLLECTION_INSERTER_HPP

