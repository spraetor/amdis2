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


/** \file MatrixStreams.h */

#ifndef AMDIS_MATRIXSTREAMS_H
#define AMDIS_MATRIXSTREAMS_H

#include "DOFMatrix.h"
#include "DOFVector.h"
#include "DOFIterator.h"
#include "solver/SolverMatrix.h"
#include "SystemVector.h"
#include "traits/inserter.h"
#include "traits/category.hpp"

namespace AMDiS
{

  // TODO reduce dependencies and remove using namespace
  using namespace mtl::operations;

  template<typename VecT, typename CurMap,
           typename Updater = update_store<typename traits::category<VecT>::value_type>>
  struct VecMap
  {
    VecT& vec;
    CurMap& mapper;
    VecMap(VecT& vec, CurMap& mapper):
      vec(vec),mapper(mapper) {}
  };

  template<typename MatT, typename CurMap>
  struct MatMap
  {
    MatT& mat;
    CurMap& mapper;
    MatMap(MatT& mat, CurMap& m):
      mat(mat), mapper(m) {}
  };

  // requires MatrixType to have an inserter
  template<typename MatrixType, typename Mapper>
  void operator<<(MatrixType& matrix, MatMap<const Matrix<DOFMatrix*>, Mapper>& Asolver)
  {
    const Matrix<DOFMatrix*>& A = Asolver.mat;
    int ns = A.getNumRows();

    Mapper& mapper(Asolver.mapper);
    set_to_zero(matrix);

    int nnz = 0;
    for (int rb = 0; rb < ns; ++rb)
      for (int cb = 0; cb < ns; ++cb)
        if (A[rb][cb])
          nnz += A[rb][cb]->getBaseMatrix().nnz();

    {
      typedef typename Mapper::template Inserter<typename Collection<MatrixType>::Inserter>::type MappedInserter;
      MappedInserter ins(matrix, mapper, int(1.2 * nnz / matrix.dim1()));
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
  template<typename MatrixType, typename Mapper>
  void operator<<(MatrixType& matrix, MatMap<const SolverMatrix<Matrix<DOFMatrix*>>, Mapper>& Asolver)
  {
    MatMap<const Matrix<DOFMatrix*>, Mapper> mapMat(*Asolver.mat.getOriginalMat(), Asolver.mapper);
    matrix << mapMat;
  }


  // requires MatrixType to have an inserter
  template<typename MatrixType, typename Mapper>
  void operator<<(MatrixType& matrix, MatMap<const Matrix<MatrixType*>, Mapper>& Asolver)
  {
    Matrix<MatrixType*>& A = *(Asolver.mat);
    int ns = A.getNumRows();

    Mapper& mapper(Asolver.mapper);
    set_to_zero(matrix);

    int nnz = 0;
    for (int rb = 0; rb < ns; ++rb)
      for (int cb = 0; cb < ns; ++cb)
        if (A[rb][cb])
          nnz += A[rb][cb]->nnz();

    {
      typedef typename Mapper::template Inserter<typename Collection<MatrixType>::Inserter>::type MappedInserter;
      MappedInserter ins(matrix, mapper, int(1.2 * nnz / matrix.dim1()));
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

  template<typename Vector, typename CurMap, typename Updater>
  void operator<<(Vector& dest, VecMap<const DOFVector<double>, CurMap, Updater>& rhs)
  {
    DOFConstIterator<double> it_x(&rhs.vec, USED_DOFS);
    size_t counter(0);
    typedef typename mtl::vector::inserter<Vector, Updater> Inserter;

    Inserter swapIns(dest);
    typedef mtl::vector::mapped_inserter<Inserter, CurMap> MappedInserter;
    MappedInserter ins(swapIns, rhs.mapper);

    for (it_x.reset(); !it_x.end(); ++it_x)
    {
      ins[counter] << *it_x;
      counter++;
    }
  }

  template<typename Vector, typename CurMap, typename Updater>
  inline void operator>>(const Vector& dest, VecMap<DOFVector<double>, CurMap, Updater>& rhs)
  {
    DOFVector<double>::Iterator it_x(&rhs.vec, USED_DOFS);
    size_t counter(0);
    {
      mtl::vector::extracter<Vector> extracter(dest);
      for (it_x.reset(); !it_x.end(); ++it_x)
      {
        extracter[rhs.mapper.row(counter)] >> *it_x ;
        counter++;
      }
    }
  }

  template<typename Vector, typename CurMap, typename Updater>
  void operator<<(Vector& dest, VecMap<const SystemVector, CurMap, Updater>& rhs)
  {
    int ns = rhs.vec.getSize();  // Number of systems.

    // Copy vectors
    typedef typename mtl::vector::inserter<Vector, Updater> Inserter;
    Inserter swapInserter(dest);
    typedef mtl::vector::mapped_inserter<Inserter, CurMap> MappedInserter;
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

  template<typename Vector, typename CurMap, typename Updater>
  void operator>>(const Vector& dest, VecMap<SystemVector, CurMap, Updater>& rhs)
  {
    int ns = rhs.vec.getSize();  // Number of systems.

    // Copy vectors
    for (int i = 0; i < ns; i++)
    {
      DOFVector<double>& dofvec(*(rhs.vec.getDOFVector(i)));
      rhs.mapper.setRow(i);
      VecMap<DOFVector<double>, CurMap> swap(dofvec, rhs.mapper);
      dest >> swap;
    }
  }
}
#endif // AMDIS_MATRIXSTREAMS_H
