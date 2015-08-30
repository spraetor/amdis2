#ifdef HAVE_SEQ_PETSC

namespace AMDiS
{

  inline void operator<<(Mat& mat, const DOFMatrix& rhs)
  {
    bool initMatrix = false;
    if (mat == PETSC_NULL)
    {
      std::vector<PetscInt> nnz(rhs.getSize());
      for (size_t k = 0; k < static_cast<size_t>(rhs.getSize()); k++)
        nnz[k] = rhs.getBaseMatrix().nnz_local(k);

      MatCreateSeqAIJ(PETSC_COMM_SELF, num_rows(rhs.getBaseMatrix()), num_cols(rhs.getBaseMatrix()), 0, &(nnz[0]), &mat);
      initMatrix = true;
    }

    std::vector<PetscInt> indices;
    for (size_t i = 0; i < rhs.getBaseMatrix().ref_minor().size(); i++)
      indices.push_back(rhs.getBaseMatrix().ref_minor()[i]);

    for (PetscInt i = 0; i < static_cast<PetscInt>(num_rows(rhs.getBaseMatrix())); i++)
    {
      if  (rhs.getBaseMatrix().nnz_local(i) > 0)
        MatSetValues(mat, 1, &i, rhs.getBaseMatrix().nnz_local(i),
                     &(indices[rhs.getBaseMatrix().ref_major()[i]]),
                     &(rhs.getBaseMatrix().value_from_offset(rhs.getBaseMatrix().ref_major()[i])), ADD_VALUES);
    }

    if (initMatrix)
    {
      MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);
    }
  }


  template<typename Mapper>
  inline void operator<<(PetscMatrix& mat, MatMap<const SolverMatrix<Matrix<DOFMatrix*>>, Mapper>& rhs)
  {
    if (mat.assembled)
      mat.destroy();

    const Matrix<DOFMatrix*>& A = *(rhs.mat.getOriginalMat());

    bool initMatrix = false;
    unsigned nRows(0);
    //     unsigned nEntries(0);
    std::vector<int> nRowsBlock(A.getNumRows());
    for(int i(0); i < A.getNumRows(); ++i)
    {
      int j(0);
      for( ; j<A.getNumCols() && A[i][j] == NULL; ++j) ;
      if( j == A.getNumCols() )
      {
        std::stringstream ss;
        ss << "ERROR: solver matrix has empty row " << i << "\n";
        throw std::runtime_error(ss.str());
      }
      nRowsBlock[i]=num_rows(A[i][j]->getBaseMatrix());
      nRows+=nRowsBlock[i];
    }
    std::vector<PetscInt> nnz(nRows);

    //initialize the nnz (could be avoided, but gets complicated)
    for (int i(0); i < nRows; ++i)
      nnz[i] = 0;

    //build a list of pairs (col, value) for each row
    std::vector<std::list<std::pair<int, double>>> rowData(nRows);

    PetscInt maxNnz(0);
    Mapper& mapper = rhs.mapper;
    for (int i(0); i < A.getNumRows(); ++i)
    {
      mapper.setRow(i);
      //compute the number of nonzeros per global row
      for (int j(0); j < A.getNumCols(); ++j)
      {
        if (A[i][j] != NULL)
        {
          mapper.setCol(j);
          const DOFMatrix::base_matrix_type& mtlMat(A[i][j]->getBaseMatrix());
          for (unsigned k(0);  k < nRowsBlock[i]; ++k)
          {
            unsigned curRow(mapper.row(k));
            nnz[curRow] += mtlMat.nnz_local(k);
            maxNnz = std::max(maxNnz, nnz[curRow]);
            unsigned curPos(mtlMat.ref_major()[k]);
            for(unsigned l(0); l < mtlMat.nnz_local(k); ++l, ++curPos)
            {
              rowData[curRow].push_back(std::make_pair(mapper.col(mtlMat.ref_minor()[curPos]), mtlMat.data[curPos]));
            }
          }
        }
      }

    }
    if (mat.matrix == PETSC_NULL)
    {
      MatCreateSeqAIJ(PETSC_COMM_SELF, nRows, nRows, 0, &(nnz[0]), &mat.matrix);
      initMatrix = true;
    }
    //reduce mallocs..
    std::vector<PetscInt> colIndices(maxNnz);
    std::vector<PetscScalar> values(maxNnz);
    std::list<std::pair<int, double>>::iterator it;
    unsigned j;
    for (PetscInt i(0); i < nRows; ++i)
    {
      j=0;
      it = rowData[i].begin();
      for (; it != rowData[i].end(); ++it, ++j)
      {
        colIndices[j] = it->first;
        values[j] = it->second;
      }
      MatSetValues(mat.matrix, 1, &i, nnz[i], &colIndices[0], &values[0], INSERT_VALUES);
    }
    if (initMatrix)
    {
      MatAssemblyBegin(mat.matrix, MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(mat.matrix, MAT_FINAL_ASSEMBLY);
    }

    mat.assembled = true;
  }


  template<typename Mapper>
  inline void operator<<(PetscMatrixNested& mat, MatMap<const SolverMatrix<Matrix<DOFMatrix*>>, Mapper>& mapper)
  {
    if (mat.assembled)
      mat.destroy();

    const Matrix<DOFMatrix*>& A = *(mapper.mat.getOriginalMat());

    std::vector<PetscInt> nnz;
    mat.nestMat.resize(A.getNumRows() * A.getNumCols());

    for (unsigned int i = 0; i < A.getNumRows(); i++)
    {
      for (unsigned int j = 0; j < A.getNumCols(); j++)
      {
        size_t idx = i * A.getNumCols() + j;

        if (!(A[i][j])
            || num_rows(A[i][j]->getBaseMatrix()) == 0
            || num_cols(A[i][j]->getBaseMatrix()) == 0
            || A[i][j]->getBaseMatrix().nnz() == 0)
        {
          mat.nestMat[idx] = PETSC_NULL;
          continue;
        }

        mat.nestMat[idx] << *(A[i][j]);
      }
    }

    // create nested matrix from a vector of block matrices
    MatCreateNest(PETSC_COMM_SELF, A.getNumRows(), PETSC_NULL, A.getNumCols(), PETSC_NULL, &(mat.nestMat[0]), &mat.matrix);
    mat.assembled = true;
  }


  inline void operator<<(Vec& petscVec, const DOFVector<double>& vec)
  {
    // Traverse all used DOFs in the dof vector.
    DOFConstIterator<double> dofIt(&vec, USED_DOFS);
    PetscInt index = 0;
    for (dofIt.reset(); !dofIt.end(); ++dofIt, ++index)
    {
      double value = *dofIt;
      VecSetValue(petscVec, index, value, ADD_VALUES);
    }
  }


  inline void operator>>(const Vec& petscVec, DOFVector<double>& vec)
  {
    // Traverse all used DOFs in the dof vector.
    DOFVector<double>::Iterator dofIt(&vec, USED_DOFS);
    PetscInt index = 0;
    for (dofIt.reset(); !dofIt.end(); ++dofIt, ++index)
    {
      double value = 0.0;
      VecGetValues(petscVec, 1, &index, &value);
      *dofIt = value;
    }
  }


  inline void operator<<(Vec& petscVec, const SystemVector& vec)
  {
#ifdef VecType // PETSc uses MACROS instead of typedefs in Versions 3.3x
    const VecType vecType;
#else
    VecType vecType;
#endif
    VecGetType(petscVec, &vecType);

    if (strcmp(vecType, VECNEST) == 0)
    {
      for (size_t i = 0; i < static_cast<size_t>(vec.getSize()); i++)
      {
        Vec v;
        VecNestGetSubVec(petscVec, i, &v);
        v << *(vec.getDOFVector(i));
      }
    }
    else
    {
      PetscInt index = 0;
      for (size_t i = 0; i < static_cast<size_t>(vec.getSize()); i++)
      {
        DOFConstIterator<double> dofIt(vec.getDOFVector(i), USED_DOFS);
        for (dofIt.reset(); !dofIt.end(); ++dofIt, ++index)
        {
          double value = *dofIt;
          VecSetValue(petscVec, index, value, ADD_VALUES);
        }
      }
    }
  }


  inline void operator>>(const Vec& petscVec, SystemVector& vec)
  {
#ifdef VecType // PETSc uses MACROS instead of typedefs in Versions 3.3x
    const VecType vecType;
#else
    VecType vecType;
#endif
    VecGetType(petscVec, &vecType);

    if (strcmp(vecType, VECNEST) == 0)
    {
      for (size_t i = 0; i < static_cast<size_t>(vec.getSize()); i++)
      {
        Vec v;
        VecNestGetSubVec(petscVec, i, &v);
        v >> *(vec.getDOFVector(i));
      }
    }
    else
    {
      PetscInt n = 0;
      for (size_t i = 0; i < static_cast<size_t>(vec.getSize()); i++)
        n += vec.getDOFVector(i)->getUsedSize();

      PetscInt N = 0;
      VecGetSize(petscVec, &N);
      assert(n == N);

      PetscInt index = 0;
      for (size_t i = 0; i < static_cast<size_t>(vec.getSize()); i++)
      {
        DOFVector<double>::Iterator dofIt(vec.getDOFVector(i), USED_DOFS);
        for (dofIt.reset(); !dofIt.end(); ++dofIt, ++index)
        {
          double value = 0.0;
          VecGetValues(petscVec, 1, &index, &value);
          *dofIt = value;
        }
      }
    }
  }


  template<typename Mapper>
  void operator>>(const PetscVector& dest, VecMap<SystemVector, Mapper>& rhs)
  {
    dest.vector >> rhs.vec;
  }


  template<typename Mapper>
  void operator>>(const PetscVectorNested& dest, VecMap<SystemVector, Mapper>& rhs)
  {
    dest.vector >> rhs.vec;
  }


  template<typename Mapper>
  inline void operator<<(PetscVector& petscVec, VecMap<const SystemVector, Mapper>& rhs)
  {
    size_t nComponents = rhs.vec.getSize();
    if (petscVec.assembled)
      petscVec.destroy();

    PetscInt n = 0;
    for (size_t i = 0; i < nComponents; i++)
      n += rhs.vec.getDOFVector(i)->getUsedSize();

    VecCreateSeq(PETSC_COMM_SELF, n, &(petscVec.vector));
    petscVec.vector << rhs.vec;

    petscVec.assembled = true;
  }


  template<typename Mapper>
  inline void operator<<(PetscVectorNested& petscVec, VecMap<const SystemVector, Mapper>& rhs)
  {
    size_t nComponents = rhs.mapper.getNumComponents();
    if (petscVec.assembled)
      petscVec.destroy();

    petscVec.nestVec.resize(nComponents);

    for (size_t i = 0; i < nComponents; i++)
    {
      VecCreateSeq(PETSC_COMM_SELF,rhs.mapper.getNumRows(i), &(petscVec.nestVec[i]));

      petscVec.nestVec[i] << *(rhs.vec.getDOFVector(i));

      VecAssemblyBegin(petscVec.nestVec[i]);
      VecAssemblyEnd(petscVec.nestVec[i]);
    }

    // create nested vector from vector of block vectors
    VecCreateNest(PETSC_COMM_SELF, nComponents, PETSC_NULL, &(petscVec.nestVec[0]), &(petscVec.vector));

    petscVec.assembled = true;
  }

} // end namespace AMDiS

#endif
