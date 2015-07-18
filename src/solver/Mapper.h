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


/** \file Mapper.h */

#ifndef AMDIS_MAPPER_H
#define AMDIS_MAPPER_H

#include "solver/SolverMatrix.h"
#include "MTL4Types.h"
#include <vector>

namespace AMDiS {
  
  /**
   * \brief BaseClass for all mapper types
   * 
   * A mapper assigned a local row/column index of a cell of
   * a block matrix to a global matrix index.
   * Call \ref setRow and \ref setCol first, to select a block
   * in the block matrix. Then get with \ref row, respective \ref col,
   * the global matrix index to the assigned local index.
   **/  
  template< typename Derived >
  struct MapperBase
  {
    typedef MTLTypes::size_type size_type;
    
    template <class BaseInserter>
    struct Inserter {
      typedef mtl::matrix::mapped_inserter<BaseInserter, Derived> type;
    };
    
    /// set the current block row
    inline void setRow( unsigned int r) { self().setRow(r); }
    
    /// set the current block columns
    inline void setCol( unsigned int c) { self().setCol(c); }
    
    /// return global matrix row, for local row in the current matrix block
    inline size_type row(size_type r) const { return self().row(r); }
    
    /// return global matrix column, for local column in the current matrix block
    inline size_type col(size_type c) const { return self().col(c); }
    
    /// return overall number of rows
    inline size_type getNumRows() const { return self().getNumRows(); }
    
    /// return overall number of columns
    inline size_type getNumCols() const { return self().getNumCols(); }

    inline size_type getNumRows(unsigned int comp) const { return self().getNumRows(comp); }
    inline size_type getNumCols(unsigned int comp) const { return self().getNumCols(comp); }
    
    /// return number of components/blocks
    inline unsigned int getNumComponents() const { return self().getNumComponents(); }
    
    inline size_type row(unsigned int block_r, unsigned int block_c, size_type r)
    {
      setRow(block_r);
      setCol(block_c);
      return row(r);
    }
        
    inline size_type col(unsigned int block_r, unsigned int block_c, size_type c)
    {
      setRow(block_r);
      setCol(block_c);
      return col(c);
    }
    
    Derived& self() { return static_cast<Derived&>(*this); }
    const Derived& self() const { return static_cast<const Derived&>(*this); }
  };

  
  /// Mapper implementation for non-parallel block matrices
  struct BlockMapper : public MapperBase< BlockMapper >
  {
    typedef unsigned int size_type;
    
    /// Default constructor
    BlockMapper()
      : nComp(0), rowOffset(0), colOffset(0), nrow(0), ncol(0), sizes(0)
    {}
    
    BlockMapper(BlockMapper const& other)
      : nComp(other.nComp),
	rowOffset(other.rowOffset),
	colOffset(other.colOffset),
	nrow(other.nrow),
	ncol(other.ncol)
    {
      sizes.resize(nComp);
      std::copy(other.sizes.begin(), other.sizes.end(), sizes.begin());
    }
    
    BlockMapper& operator=(BlockMapper const& other)
    {
      nComp = other.getNumComponents();
      rowOffset = other.row(0);
      colOffset = other.col(0);
      nrow = other.getNumRows();
      ncol = other.getNumCols();
    
      sizes.resize(nComp);
      for (size_t i = 0; i < nComp; ++i)
	sizes[i] = other.getNumRows(i);
	
      return *this;
    }
    
    /// Constructor for block-matrices
    BlockMapper(const SolverMatrix<Matrix<DOFMatrix* > >& sm )
    : nComp(sm.getOriginalMat()->getSize()), 
      rowOffset(0), colOffset(0), nrow(0), ncol(0), sizes(nComp)
    {
      const Matrix<DOFMatrix* >& orMat(*sm.getOriginalMat());
      const int ns = orMat.getNumRows();
      for (int i= 0; i < ns; i++) {
	sizes[i] = orMat[i][i]->getFeSpace()->getAdmin()->getUsedSize();
	nrow += sizes[i];
      }
      ncol = nrow;
    }
    
    /// Constructor for block-matrices
    BlockMapper(const Matrix<DOFMatrix*>& orMat )
    : nComp(orMat.getSize()), 
      rowOffset(0), colOffset(0), nrow(0), ncol(0), sizes(nComp)
    {
      const int ns = orMat.getNumRows();
      for (int i= 0; i < ns; i++) {
	sizes[i] = orMat[i][i]->getFeSpace()->getAdmin()->getUsedSize();
	nrow += sizes[i];
      }
      ncol = nrow;
    }
    
    /// Constructor for single matrix
    BlockMapper(const DOFMatrix* sm )
    : nComp(1), 
      rowOffset(0), colOffset(0), 
      nrow(0), ncol(0), 
      sizes(nComp)
    {
      sizes[0] = sm->getFeSpace()->getAdmin()->getUsedSize();
      nrow += sizes[0];
      ncol = nrow;
    }

    /// Constructor for system with equal components
    BlockMapper(unsigned int nComp, size_type nDOFperComp)
    : nComp(nComp), 
      rowOffset(0), colOffset(0), 
      nrow(nComp*nDOFperComp), ncol(nrow),
      sizes(nComp)
    {
       for (unsigned int i = 0; i < nComp; ++i)
         sizes[i] = nDOFperComp;
    }
      
    /// calculate row offset for row component \param r
    inline void setRow( unsigned int r)
    { 
      assert( r <= sizes.size() );
      rowOffset = sum(r); 
    }

    /// calculate column offset for col component \param c
    inline void setCol( unsigned int c)
    {
      assert( c <= sizes.size() ); 
      colOffset = sum(c);
    }

    inline size_type row(size_type r) const { return r + rowOffset; }
    inline size_type col(size_type c) const { return c + colOffset; }

    inline size_type getNumRows() const { return nrow; }
    inline size_type getNumCols() const { return ncol; }

    inline size_type getNumRows(unsigned int comp) const { return sizes[comp]; }
    inline size_type getNumCols(unsigned int comp) const { return sizes[comp]; }

    inline unsigned int getNumComponents() const { return nComp; }

  private: // methods

    ///compute the sum of sizes from [0, end)
    inline size_type sum(unsigned int end) const 
    {
      unsigned int ret(0);
      for (unsigned int i(0); i < end; ++i)
	ret += sizes[i];
      return ret;
    }
    
  private: // variables

    unsigned int nComp;
    size_type rowOffset;
    size_type colOffset;
    unsigned int nrow;
    unsigned int ncol;
    
    std::vector< size_type > sizes;
  };
  
  

  
  /// Mapper implementation for non-parallel rectangular block matrices
  struct RectangularMapper : public MapperBase< RectangularMapper >
  {
    typedef unsigned int size_type;
    
    /// Constructor for block-matrices
    RectangularMapper(const SolverMatrix<Matrix<DOFMatrix* > >& sm )
    : nRowComp(sm.getOriginalMat()->getNumRows()), 
      nColComp(sm.getOriginalMat()->getNumCols()), 
      rowOffset(0), colOffset(0), nrow(0), ncol(0), 
      sizes_rows(nRowComp), sizes_cols(nColComp)
    {
      const Matrix<DOFMatrix* >& orMat(*sm.getOriginalMat());

      for (unsigned int i = 0; i < nRowComp; i++) {
	for (unsigned int j = 0; j < nColComp; j++) {
	  if (orMat[i][j]) {
	    sizes_rows[i] = orMat[i][j]->getRowFeSpace()->getAdmin()->getUsedSize();
	    sizes_cols[j] = orMat[i][j]->getColFeSpace()->getAdmin()->getUsedSize();
	  }
	}
	nrow += sizes_rows[i];
      }
      for (unsigned int j = 0; j < nColComp; j++)
	ncol += sizes_cols[j];
    }
      
    /// calculate row offset for row component \param r
    inline void setRow( unsigned int r)
    { 
      assert( r <= sizes_rows.size() );
      rowOffset = sum(r, sizes_rows); 
    }

    /// calculate column offset for col component \param c
    inline void setCol( unsigned int c)
    {
      assert( c <= sizes_cols.size() ); 
      colOffset = sum(c, sizes_cols);
    }

    inline size_type row(size_type r) const { return r + rowOffset; }
    inline size_type col(size_type c) const { return c + colOffset; }

    inline size_type getNumRows() const { return nrow; }
    inline size_type getNumCols() const { return ncol; }

    inline size_type getNumRows(unsigned int comp) const { return sizes_rows[comp]; }
    inline size_type getNumCols(unsigned int comp) const { return sizes_cols[comp]; }

    inline unsigned int getNumComponents() const { return nRowComp; }
    inline unsigned int getNumRowComponents() const { return nRowComp; }
    inline unsigned int getNumColComponents() const { return nColComp; }

  private: // methods

    ///compute the sum of sizes from [0, end)
    inline size_type sum(unsigned int end, std::vector< size_type >& sizes) const 
    {
      unsigned int ret(0);
      for (unsigned int i(0); i < end; ++i)
	ret += sizes[i];
      return ret;
    }
    
  private: // variables

    unsigned int nRowComp, nColComp;
    size_type rowOffset;
    size_type colOffset;
    unsigned int nrow;
    unsigned int ncol;
    
    std::vector< size_type > sizes_rows;
    std::vector< size_type > sizes_cols;
  };
}

#endif // AMDIS_MAPPER_H
