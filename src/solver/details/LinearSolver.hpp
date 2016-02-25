#pragma once

#include <iostream>

#include <boost/mpl/bool.hpp>
#include <boost/numeric/mtl/utility/is_distributed.hpp>

#include "Timer.hpp"
#include "solver/MatrixStreams.hpp"

#ifdef HAVE_PARALLEL_MTL4
#include <boost/numeric/mtl/par/distribution.hpp>
#endif

namespace AMDiS
{

  namespace tag
  {
    struct distributed {};
    struct non_distributed {};
  }

  namespace traits
  {
    template<class T>
    struct distributed_tag : boost::mpl::if_<
      mtl::traits::is_distributed<T>,
      tag::distributed,
      tag::non_distributed
      > {};
  }

  namespace dispatch
  {
    /// init matrix
    template<typename MatrixOut, typename M>
    void initMatrix(MatrixOut& /*m*/, MapperBase<M>& /*mapper*/) {}

    /// init vector
    template<typename VectorOut, typename MatrixT>
    void initVector(VectorOut& /*v*/, const MatrixT& /*source*/) {}


    /// init systemmatrix depending on Mapper parameters.
    template<typename M>
    void initMatrix(MTLTypes::MTLMatrix& m, MapperBase<M>& mapper)
    {
      m.change_dim(mapper.getNumRows(), mapper.getNumCols());
      set_to_zero(m);
    }

    template<typename MatrixOut, typename MatrixIn, typename M>
    void fillMatrix(MatrixOut& m, const MatrixIn& source, MapperBase<M>& mapper)
    {
      MatMap<const MatrixIn, M> matMap{source, mapper.self()};
      m << matMap;
    }

#ifdef HAVE_PARALLEL_MTL4
    /// init systemmatrix depending on Mapper parameters,
    /// specialized for distributed matrices
    template<typename M>
    void initMatrix(MTLTypes::PMTLMatrix& m, MapperBase<M>& mapper)
    {
      mtl::par::block_distribution dist(mapper.getNumRows());
      dist.setup_from_local_size(mapper.getMap().getLocalDofs());
      m.change_dim(0, 0);
      m.init_distribution(dist, dist, mapper.getNumRows(), mapper.getNumRows());
      set_to_zero(m);
    }

    /// init MTL-vector depending on Mapper parameters,
    /// specialized for distributed matrices
    template<typename MatrixT>
    void initVector(MTLTypes::PMTLVector& v, const MatrixT& matrix)
    {
      v.init_distribution(row_distribution(matrix), num_rows(matrix));
      set_to_zero(v);
    }
#endif

    /// init MTL-vector depending on Mapper parameters.
    template<typename MatrixT>
    void initVector(MTLTypes::MTLVector& v, const MatrixT& matrix)
    {
      v.change_dim(num_rows(matrix));
      set_to_zero(v);
    }

    /// fill MTL-vector
    template<typename VectorOut, typename VectorIn, typename M>
    void fillVector(VectorOut& v, const VectorIn& source, MapperBase<M>& mapper)
    {
      VecMap<const VectorIn, M> srcVecMap{source, mapper.self()};
      v << srcVecMap;
    }



  } // end namespace dispatch
} // end namespace AMDiS
