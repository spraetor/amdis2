/** \file MatrixStreams.h */

#pragma once

// AMDiS includes
#include "DOFIterator.hpp"
#include "DOFMatrix.hpp"
#include "DOFVector.hpp"
#include "SystemVector.hpp"
#include "solver/SolverMatrix.hpp"
#include "traits/category.hpp"
#include "traits/inserter.hpp"

namespace AMDiS
{
  template <class VectorType, class Mapper,
            class Updater = mtl::operations::update_store<typename traits::category<VectorType>::value_type>>
  struct VecMap
  {
    VectorType& vec;
    Mapper& mapper;
  };

  template <class MatrixType, class Mapper>
  struct MatMap
  {
    MatrixType& mat;
    Mapper& mapper;
  };
  
  // copy System-matrix to MTL-matrix
  // ---------------------------------------------------------------------------

  // requires MatrixType to have an inserter
  template <class MatrixType, class Mapper>
  void operator<<(MatrixType& matrix, MatMap<const Matrix<DOFMatrix*>, Mapper>& Asolver)
  {
    using Ins = typename Collection<MatrixType>::Inserter;
    using MappedInserter = typename Mapper::template Inserter<Ins>::type;
    
    const Matrix<DOFMatrix*>& A = Asolver.mat;
    int ns = A.getNumRows();

    Mapper& mapper(Asolver.mapper);
    set_to_zero(matrix);

    size_t nnz = 0;
    for (int rb = 0; rb < ns; ++rb)
      for (int cb = 0; cb < ns; ++cb)
        if (A[rb][cb])
          nnz += A[rb][cb]->getBaseMatrix().nnz();

    {
      MappedInserter ins(matrix, mapper, size_t(1.2 * nnz / matrix.dim1()));
      for (int rb = 0; rb < ns; ++rb)
      {
        mapper.setRow(rb);
        for (int cb = 0; cb < ns; ++cb)
        {
          mapper.setCol(cb);
          if (A[rb][cb])
            ins << A[rb][cb]->getBaseMatrix();
        }
      }
    }
  }

  // requires MatrixType to have an inserter
  template <class MatrixType, class Mapper>
  void operator<<(MatrixType& matrix, MatMap<const SolverMatrix<Matrix<DOFMatrix*>>, Mapper>& Asolver)
  {
    MatMap<const Matrix<DOFMatrix*>, Mapper> mapMat{*Asolver.mat.getOriginalMat(), Asolver.mapper};
    matrix << mapMat;
  }


  // requires MatrixType to have an inserter
  template <class MatrixType, class Mapper>
  void operator<<(MatrixType& matrix, MatMap<const Matrix<MatrixType*>, Mapper>& Asolver)
  {
    using Ins = typename Collection<MatrixType>::Inserter;
    using MappedInserter = typename Mapper::template Inserter<Ins>::type;
      
    Matrix<MatrixType*>& A = *(Asolver.mat);
    int ns = A.getNumRows();

    Mapper& mapper(Asolver.mapper);
    set_to_zero(matrix);

    size_t nnz = 0;
    for (int rb = 0; rb < ns; ++rb)
      for (int cb = 0; cb < ns; ++cb)
        if (A[rb][cb])
          nnz += A[rb][cb]->nnz();

    {
      MappedInserter ins(matrix, mapper, size_t(1.2 * nnz / matrix.dim1()));
      for (int rb = 0; rb < ns; ++rb)
      {
        mapper.setRow(rb);
        for (int cb = 0; cb < ns; ++cb)
        {
          mapper.setCol(cb);
          if (A[rb][cb])
            ins << (*A[rb][cb]);
        }
      }
    }
  }
  
  
  // copy System-vector and DOFVector to MTL-vector and vice versa
  // ---------------------------------------------------------------------------

  template <class Vector, class Mapper, class Updater>
  void operator<<(Vector& dest, VecMap<const DOFVector<double>, Mapper, Updater>& rhs)
  {
    using Inserter = typename mtl::vector::inserter<Vector, Updater>;
    using MappedInserter = mtl::vector::mapped_inserter<Inserter, Mapper>;
    
    DOFConstIterator<double> it_x(&rhs.vec, USED_DOFS);
    size_t counter = 0;

    Inserter swapIns(dest);
    MappedInserter ins(swapIns, rhs.mapper);

    for (it_x.reset(); !it_x.end(); ++it_x)
    {
      ins[counter] << *it_x;
      counter++;
    }
  }

  template <class Vector, class Mapper, class Updater>
  void operator>>(Vector const& source, VecMap<DOFVector<double>, Mapper, Updater>& rhs)
  {
    DOFVector<double>::Iterator it_x(&rhs.vec, USED_DOFS);
    {
      mtl::vector::extracter<Vector> extracter(source);
      size_t counter = 0;
      for (it_x.reset(); !it_x.end(); ++it_x)
      {
        extracter[rhs.mapper.row(counter)] >> *it_x ;
        counter++;
      }
    }
  }

  template <class Vector, class Mapper, class Updater>
  void operator<<(Vector& dest, VecMap<const SystemVector, Mapper, Updater>& rhs)
  {
    using Inserter = typename mtl::vector::inserter<Vector, Updater>;
    using MappedInserter = mtl::vector::mapped_inserter<Inserter, Mapper>;
    
    int ns = rhs.vec.getSize();  // Number of systems.

    // Copy vectors
    Inserter swapInserter(dest);
    MappedInserter ins(swapInserter, rhs.mapper);
    for (int i = 0; i < ns; i++)
    {
      DOFConstIterator<double> it_source(rhs.vec.getDOFVector(i), USED_DOFS);
      size_t counter(0);
      rhs.mapper.setRow(i);
      for (it_source.reset(); !it_source.end(); ++it_source)
      {
        ins[counter] << *it_source;
        counter++;
      }
    }
  }

  template <class Vector, class Mapper, class Updater>
  void operator>>(Vector const& source, VecMap<SystemVector, Mapper, Updater>& rhs)
  {
    int ns = rhs.vec.getSize();  // Number of systems.

    // Copy vectors
    for (int i = 0; i < ns; i++)
    {
      DOFVector<double>& dofvec(*(rhs.vec.getDOFVector(i)));
      rhs.mapper.setRow(i);
      VecMap<DOFVector<double>, Mapper> swap{dofvec, rhs.mapper};
      source >> swap;
    }
  }
}
