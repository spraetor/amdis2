#pragma once

#include <boost/numeric/mtl/mtl.hpp>

// TODO: find replacement for this explicit typedefs
namespace AMDiS
{
  namespace MTLTypes
  {
    typedef double value_type;
    typedef unsigned long size_type;
    typedef mtl::matrix::parameters<mtl::row_major, mtl::index::c_index, mtl::non_fixed::dimensions, false, size_type> para;
    typedef mtl::matrix::compressed2D<value_type, para> MTLMatrix;
    typedef mtl::vector::dense_vector<value_type> MTLVector;

#if defined(HAVE_PARALLEL_MTL4)
    typedef mtl::matrix::distributed<MTLMatrix> PMTLMatrix;
    typedef mtl::vector::distributed<MTLVector> PMTLVector;
#endif
  }

} // endnamespace AMDiS
