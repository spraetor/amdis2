/** \file inserter.hpp */

#pragma once

// AMDiS headers
#include "MTL4Types.hpp"

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
  
} // end namespace AMDiS
