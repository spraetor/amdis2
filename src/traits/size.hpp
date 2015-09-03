/** \file size.hpp */

#pragma once

#include <AMDiS_fwd.h>
#include <traits/traits_fwd.hpp>

namespace AMDiS
{
  namespace traits
  {
    template <typename T>
    struct size<T, typename enable_if<is_scalar<T>>::type>
    {
      typedef size_t type;
      type operator()(const T& v) const
      {
        return 1;
      }
    };

    // == vectors ===

    template <typename T>
    struct size<T, typename enable_if<and_<is_vector<T>, is_mtl<T>>>::type>
    {
      typedef typename T::size_type type;
      type operator()(const T& v) const
      {
        return ::mtl::vector::size(v);
      }
    };

    /// size implementation for AMDiS::DOFVector
    template <typename Value>
    struct size<AMDiS::DOFVector<Value>>
    {
      typedef typename category<AMDiS::DOFVector<Value>>::size_type   type;
      type operator()(const AMDiS::DOFVector<Value>& v) const
      {
        return v.getUsedSize();
      }
    };

    // === matrices ===

    template <typename T>
    struct size<T, typename enable_if<and_<is_matrix<T>, is_mtl<T>>>::type>
    {
      typedef typename T::size_type type;
      type operator()(const T& v) const
      {
        return ::mtl::matrix::size(v);
      }
    };

    // === multi-vector or multi-matrix ===

    /// size implementation for AMDiS::VectorOfFixVecs
    template <typename FixVecType>
    struct size<AMDiS::VectorOfFixVecs<FixVecType>>
    {
      typedef typename category<AMDiS::VectorOfFixVecs<FixVecType>>::size_type   type;
      type operator()(const AMDiS::VectorOfFixVecs<FixVecType>& v) const
      {
        return v.getSize() * v.getSizeOfFixVec();
      }
    };

    /// AMDiS::MatrixOfFixVecs
    template <typename FixVecType>
    struct size<AMDiS::MatrixOfFixVecs<FixVecType>>
    {
      typedef typename category<AMDiS::MatrixOfFixVecs<FixVecType>>::size_type   type;
      type operator()(const AMDiS::MatrixOfFixVecs<FixVecType>& v) const
      {
        return v.getNumberOfRows() == 0 ? 0 : v.getNumberOfRows() * v.getNumberOfColumns() * v[0].getSizeOfFixVec();
      }
    };

  } // end namespace traits

} // end namespace AMDiS
