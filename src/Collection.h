/** \file Collection.h */

#pragma once

#include "MTL4Types.h"

namespace AMDiS {
  template<typename T>
  struct Collection {};

  template<>
    struct Collection< MTLTypes::MTLMatrix > {
      typedef mtl::matrix::inserter< MTLTypes::MTLMatrix > Inserter;
    };

  template<>
    struct Collection< MTLTypes::MTLVector > {
      typedef mtl::vector::inserter< MTLTypes::MTLVector > Inserter;
      typedef MTLTypes::MTLMatrix PreconditionMatrix;
    };
    
} // end namespace AMDiS
