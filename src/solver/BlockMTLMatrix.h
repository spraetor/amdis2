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
 

/** \file BlockMTLMatrix.h */

#ifndef AMDIS_BLOCK_MTL_MATRIX_H
#define AMDIS_BLOCK_MTL_MATRIX_H

#include "solver/SolverMatrix.h"
#include "solver/Mapper.h"
#include "MTL4Types.h"

namespace AMDiS {
  
  /// A wrapper for AMDiS::SolverMatrix to be used in MTL/ITL solvers
  struct BlockMTLMatrix
  {
    typedef mtl::Collection<MTLTypes::MTLMatrix>::size_type size_type;
    
    size_type n; // overall number of rows
    size_type m; // overall number of columns
    
    size_t n_rows; // number of row blocks
    size_t n_cols; // number of column blocks
    
    BlockMTLMatrix()
      : initialized(false), A(NULL) 
    { }
    
    BlockMTLMatrix(const SolverMatrix<Matrix<DOFMatrix*> >& A_)
      : initialized(false), A(&A_)
    {
      init();
    }
    
    void setMatrix(const SolverMatrix<Matrix<DOFMatrix*> >& A_) 
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
      for (size_t i = 0; i < n_rows; i++) {
	mapper.setRow(i+1);
	int finish = mapper.row(0);
	r_rows[i].set(start, finish);
	start = finish;
      }
      start = 0;
      mapper.setRow(0);
      for (size_t i = 0; i < n_cols; i++) {
	mapper.setCol(i+1);
	int finish = mapper.col(0);
	r_cols[i].set(start, finish);
	start = finish;
      }
      initialized = true;
    }
    
    inline const mtl::irange& getRowRange(size_t i) const
    {
      return r_rows[i];
    } 
    
    inline const mtl::irange& getColRange(size_t i) const
    {
      return r_cols[i];
    } 
    
    const DOFMatrix::base_matrix_type& getSubMatrix(size_t i, size_t j) const 
    {
      return A->getSubMatrix(i, j);
    }

    /// perform blockwise multiplication A*b -> x
    template <typename VectorIn, typename VectorOut, typename Assign>
    void mult(const VectorIn& b, VectorOut& x, Assign) const
    {
      TEST_EXIT_DBG(initialized)("Block matrix not initialized. Assign a SolverMatrix first!\n");
      
      for (size_t i = 0; i < n_rows; i++) {
	VectorOut x_i(x[r_rows[i]]);
	bool first = true;
	for (size_t j = 0; j < n_cols; j++) {
	  if ((*A->getOriginalMat())[i][j]) {
	    const VectorIn b_j(b[r_cols[j]]);
	    if (first) {
	      Assign::first_update(x_i, A->getSubMatrix(i, j) * b_j);
	      first = false;
	    } else {
	      Assign::update(x_i, A->getSubMatrix(i, j) * b_j);
	    }
	  }
	}
      }
    }

    template <typename VectorIn>
    mtl::vector::mat_cvec_multiplier<BlockMTLMatrix, VectorIn> operator*(const VectorIn& v) const
    {   return mtl::vector::mat_cvec_multiplier<BlockMTLMatrix, VectorIn>(*this, v);    }

  protected:
    bool initialized;
    
  private:
    const SolverMatrix<Matrix<DOFMatrix*> >* A;    
    std::vector<mtl::irange> r_rows, r_cols;
  };
  
  namespace dispatch
  {
    template< typename M >
    void initMatrix(BlockMTLMatrix& m, MapperBase<M>& mapper)
    { }
    
    template< typename M >
    void fillMatrix(BlockMTLMatrix& m, const SolverMatrix<Matrix<DOFMatrix*> >& source, MapperBase<M>& mapper)
    {
      m.setMatrix(source);
    }
    
  } // end namespace dispatch

} // end namespace AMDiS


inline AMDiS::BlockMTLMatrix::size_type size(const AMDiS::BlockMTLMatrix& A) { return (A.n) * (A.m); }
inline AMDiS::BlockMTLMatrix::size_type num_rows(const AMDiS::BlockMTLMatrix& A) { return A.n; }
inline AMDiS::BlockMTLMatrix::size_type num_cols(const AMDiS::BlockMTLMatrix& A) { return A.m; }

namespace mtl 
{
  template <>
  struct Collection<AMDiS::BlockMTLMatrix>
  {
    typedef double value_type;
    typedef AMDiS::BlockMTLMatrix::size_type size_type;
  };

  namespace ashape 
  {
    template <> struct ashape_aux<AMDiS::BlockMTLMatrix> 
    { 
      typedef nonscal type;
    };
  } // end namespace ashape
  
} // end namespace mtl

#endif // AMDIS_BLOCK_MTL_MATRIX_H