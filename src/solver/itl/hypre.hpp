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

// Written by Andreas Naumann


#ifndef AMDIS_HYPRE_SOLVER_HPP
#define AMDIS_HYPRE_SOLVER_HPP

#ifdef MTL_HAS_HYPRE

#include "mpi.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

namespace mtl
{

  class HypreMatrix
  {
  public:
    HypreMatrix() : initialized(false) {}

    template<typename Matrix>
    HypreMatrix(const Matrix& matrix)
    {
      init(matrix);
    }

    template<typename Matrix>
    void init(const Matrix& matrix)
    {
      MPI_Comm comm = MPI_COMM_WORLD;
      int ilower(0);
      int jlower(ilower);
      int nRows(num_rows(matrix));
      int nCols(num_cols(matrix));
      int iupper(ilower + nRows-1);
      int jupper(jlower + nCols-1);
      HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &ij_matrix);
      HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &ij_matrix);
      HYPRE_IJMatrixSetObjectType(ij_matrix, HYPRE_PARCSR);
      HYPRE_IJMatrixInitialize(ij_matrix);
      int* ncols = new int[nRows];
      int* rows = new int[nRows];
      int nnz(0);
      for(int i(0); i < nRows; ++i)
      {
        ncols[i] = matrix.nnz_local(i);
        rows[i] = i;
        nnz += ncols[i];
      }
      int* cols = new int[nnz];
      const typename Matrix::size_type* minors = matrix.address_minor();
      for(int i(0); i < nnz; ++i)
      {
        cols[i] = minors[i];
      }
      const double* values = matrix.address_data();
      HYPRE_IJMatrixSetValues(ij_matrix, nRows, ncols, rows, cols, values);
      HYPRE_IJMatrixAssemble(ij_matrix);
      delete [] ncols;
      delete [] rows;
      delete [] cols;

      initialized = true;
    }

    ~HypreMatrix()
    {
      if(initialized)
        HYPRE_IJMatrixDestroy(ij_matrix);
      initialized = false;
    }

    inline operator HYPRE_IJMatrix()
    {
      return ij_matrix;
    }

  private:
    HYPRE_IJMatrix ij_matrix;
    bool initialized;
  };

  class HypreVector
  {
  public:
    template<typename Vector>
    HypreVector(const Vector& mtlVec, int jlower = 0)
      : jlower(jlower)
    {
      int jupper(num_rows(mtlVec)-1 + jlower);
      HYPRE_IJVectorCreate(MPI_COMM_WORLD, jlower, jupper, &vec);
      HYPRE_IJVectorSetObjectType(vec, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize(vec);
      int* vecIndices = new int[num_rows(mtlVec)];
      for (size_t i(0); i < num_rows(mtlVec); ++i)
        vecIndices[i] = i+jlower;
      HYPRE_IJVectorSetValues(vec, num_rows(mtlVec), vecIndices, mtlVec.address_data());
      HYPRE_IJVectorAssemble(vec);
      delete [] vecIndices;
      HYPRE_ParVector x;
      HYPRE_IJVectorGetObject(vec, (void**) &x);
    }

    ~HypreVector()
    {
      HYPRE_IJVectorDestroy(vec);
    }

    inline operator HYPRE_IJVector()
    {
      return vec;
    }

  public:
    int jlower;

  private:
    HYPRE_IJVector vec;
  };

  template<typename Vector>
  void convert(HypreVector& src, Vector& dest)
  {
    int* indices = new int[num_rows(dest)];
    for (size_t i(0); i < num_rows(dest); ++i)
      indices[i] = src.jlower + i;
    HYPRE_IJVectorGetValues(src, num_rows(dest), indices, dest.address_data());
    delete [] indices;
  }

  class HypreParVector
  {
  public:
    template<typename MtlVector>
    HypreParVector(const MtlVector& mtlVec)
      : vec2(mtlVec)
    {
      HYPRE_IJVectorGetObject(vec2, (void**) &vec);
    }

    inline operator HYPRE_ParVector()
    {
      return vec;
    }
    inline HypreVector& getHypreVector()
    {
      return vec2;
    }
    inline const HypreVector& getHypreVector() const
    {
      return vec2;
    }

  private:
    HYPRE_ParVector vec;
    HypreVector vec2;
  };

  class AMGConfigurator
  {
  public:
    enum CycleDirection { DOWN=1, UP=2, COARSEST=3 };

    AMGConfigurator()
      : _debugFlag(0), _printLevel(0), _logLevel(0),
        _tolerance(1.0e-7), _cycle(-1), _finestSweep(1),
        _interpolation(-1), _relaxation(-1), _maxIter(20)
    {
      for(int i(2); i>= 0; --i)
        _sweeps[i] = 1;
    }

    AMGConfigurator& debugFlag(int df)
    {
      _debugFlag = df;
      return *this;
    }
    AMGConfigurator& printLevel(int pl)
    {
      _printLevel = pl;
      return *this;
    }
    AMGConfigurator& logLevel(int ll)
    {
      _logLevel = ll;
      if (ll > 1) throw std::runtime_error("logging > 1 not supported\n2");
      return *this;
    }
    AMGConfigurator& tolerance(double tol)
    {
      _tolerance = tol;
      return *this;
    }
    AMGConfigurator& cycle(int c)
    {
      _cycle = c;
      return *this;
    }
    AMGConfigurator& finestSweep(int nr)
    {
      _finestSweep = nr;
      return *this;
    }

    AMGConfigurator& sweep(CycleDirection dir, int nr)
    {
      _sweeps[dir-1] = nr;
      return *this;
    }

    AMGConfigurator& interpolationType(int type)
    {
      _interpolation = type;
      return *this;
    }
    AMGConfigurator& relaxType(int type)
    {
      _relaxation = type;
      return *this;
    }
    AMGConfigurator& maxIter(int iter)
    {
      _maxIter = iter;
      return *this;
    }

    inline void operator()(HYPRE_Solver& solver) const
    {
      HYPRE_BoomerAMGSetDebugFlag(solver, _debugFlag);
      HYPRE_BoomerAMGSetPrintLevel(solver, _printLevel);
      HYPRE_BoomerAMGSetLogging(solver, _logLevel);
      HYPRE_BoomerAMGSetTol(solver, _tolerance);
      if (_cycle >= 0)
        HYPRE_BoomerAMGSetCycleType(solver, _cycle);
      HYPRE_BoomerAMGSetNumSweeps(solver, _finestSweep);
      for(int i(3); i >= 1; --i)
      {
        HYPRE_BoomerAMGSetCycleNumSweeps(solver, _sweeps[i-1], i);
      }
      if (_interpolation >= 0)
        HYPRE_BoomerAMGSetInterpType(solver, _interpolation);
      if (_relaxation >= 0)
        HYPRE_BoomerAMGSetRelaxType(solver, _relaxation);

      HYPRE_BoomerAMGSetMaxIter(solver, _maxIter);
    }

  protected:
    int _debugFlag;
    int _printLevel;
    int _logLevel;
    double _tolerance;
    int _cycle;
    int _finestSweep;
    int _sweeps[3];
    int _interpolation;
    int _relaxation;
    int _maxIter;
  };

} //namespace mtl

#endif // MTL_HAS_HYPRE

#endif // HYPRE_SOLVER_HPP
