/** \file BlockMTLMatrix.h */

#pragma once

#include "MTL4Types.hpp"
#include "solver/SolverMatrix.hpp"
#include "solver/Mapper.hpp"

namespace AMDiS
{
  /// A wrapper for AMDiS::SolverMatrix to be used in MTL/ITL solvers
  struct BlockMTLMatrix
  {
    using size_type = mtl::Collection<MTLTypes::MTLMatrix>::size_type;

    size_type n; // overall number of rows
    size_type m; // overall number of columns

    size_t n_rows; // number of row blocks
    size_t n_cols; // number of column blocks

    BlockMTLMatrix()
      : initialized(false), 
	A(NULL)
    { }

    BlockMTLMatrix(SolverMatrix<Matrix<DOFMatrix*>> const& A_)
      : initialized(false), 
	A(&A_)
    {
      init();
    }

    void setMatrix(SolverMatrix<Matrix<DOFMatrix*>> const& A_)
    {
      A = &A_;
      init();
    }

    void init()
    {
      RectangularMapper mapper(*A);
      n = mapper.getNumRows();
      m = mapper.getNumCols();
      n_rows = mapper.getNumRowComponents();
      n_cols = mapper.getNumColComponents();

      r_rows.resize(n_rows);
      r_cols.resize(n_cols);
      int start = 0;
      mapper.setCol(0);
      for (size_t i = 0; i < n_rows; i++)
      {
        mapper.setRow(i+1);
        int finish = mapper.row(0);
        r_rows[i].set(start, finish);
        start = finish;
      }
      start = 0;
      mapper.setRow(0);
      for (size_t i = 0; i < n_cols; i++)
      {
        mapper.setCol(i+1);
        int finish = mapper.col(0);
        r_cols[i].set(start, finish);
        start = finish;
      }
      initialized = true;
    }

    mtl::irange const& getRowRange(size_t i) const
    {
      return r_rows[i];
    }

    mtl::irange const& getColRange(size_t i) const
    {
      return r_cols[i];
    }

    DOFMatrix::base_matrix_type const& getSubMatrix(size_t i, size_t j) const
    {
      return A->getSubMatrix(i, j);
    }

    /// perform blockwise multiplication A*b -> x
    template <class VectorIn, class VectorOut, class Assign>
    void mult(VectorIn const& b, VectorOut& x, Assign) const
    {
      TEST_EXIT_DBG(initialized)
	("Block matrix not initialized. Assign a SolverMatrix first!\n");

      for (size_t i = 0; i < n_rows; i++)
      {
        VectorOut x_i(x[r_rows[i]]);
        bool first = true;
        for (size_t j = 0; j < n_cols; j++)
        {
          if ((*A->getOriginalMat())[i][j])
          {
            const VectorIn b_j(b[r_cols[j]]);
            if (first)
            {
              Assign::first_update(x_i, A->getSubMatrix(i, j) * b_j);
              first = false;
            }
            else
            {
              Assign::update(x_i, A->getSubMatrix(i, j) * b_j);
            }
          }
        }
      }
    }

    template <class VectorIn>
    mtl::vector::mat_cvec_multiplier<BlockMTLMatrix, VectorIn> 
    operator*(VectorIn const& v) const
    {
      return {*this, v};
    }

  protected:
    bool initialized;

  private:
    SolverMatrix<Matrix<DOFMatrix*>> const* A;
    std::vector<mtl::irange> r_rows, r_cols;
  };


  namespace dispatch
  {
    template<typename M>
    void initMatrix(BlockMTLMatrix& m, MapperBase<M>& mapper)
    { }

    template<typename M>
    void fillMatrix(BlockMTLMatrix& m, 
		    SolverMatrix<Matrix<DOFMatrix*>> const& source, 
		    MapperBase<M>& mapper)
    {
      m.setMatrix(source);
    }

  } // end namespace dispatch


  namespace traits
  {
    template <>
    struct category<BlockMTLMatrix>
    {
      using tag        = tag::matrix;
      using value_type = MTLTypes::value_type;
      using size_type  = MTLTypes::size_type;
    };
  }


  inline size_t size(BlockMTLMatrix const& A)
  {
    return (A.n) * (A.m);
  }

  inline size_t num_rows(BlockMTLMatrix const& A)
  {
    return A.n;
  }

  inline size_t num_cols(BlockMTLMatrix const& A)
  {
    return A.m;
  }

  inline void set_to_zero(BlockMTLMatrix& A) {}

} // end namespace AMDiS

namespace mtl
{
  template <>
  struct Collection<AMDiS::BlockMTLMatrix>
  {
    using value_type = AMDiS::MTLTypes::value_type;
    using size_type  = AMDiS::MTLTypes::size_type;
  };

  namespace ashape
  {
    template <> 
    struct ashape_aux<AMDiS::BlockMTLMatrix>
    {
      using type = nonscal;
    };

  } // end namespace ashape

} // end namespace mtl
