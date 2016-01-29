/** \file Mapper.h */

#pragma once

#include "solver/SolverMatrix.h"
#include "MTL4Types.h"
#include <vector>

namespace AMDiS
{
  /**
   * \brief BaseClass for all mapper types
   *
   * A mapper assigned a local row/column index of a cell of
   * a block matrix to a global matrix index.
   * Call \ref setRow and \ref setCol first, to select a block
   * in the block matrix. Then get with \ref row, respective \ref col,
   * the global matrix index to the assigned local index.
   **/
  template <class Derived>
  struct MapperBase
  {
    typedef MTLTypes::size_type size_type;

    template <class BaseInserter>
    struct Inserter
    {
      typedef mtl::matrix::mapped_inserter<BaseInserter, Derived> type;
    };

    /// set the current block row
    void setRow( unsigned int r)
    {
      self().setRow(r);
    }

    /// set the current block columns
    void setCol( unsigned int c)
    {
      self().setCol(c);
    }

    /// return global matrix row, for local row in the current matrix block
    size_type row(size_type r) const
    {
      return self().row(r);
    }

    /// return global matrix column, for local column in the current matrix block
    size_type col(size_type c) const
    {
      return self().col(c);
    }

    /// return overall number of rows
    size_type getNumRows() const
    {
      return self().getNumRows();
    }

    /// return overall number of columns
    size_type getNumCols() const
    {
      return self().getNumCols();
    }

    size_type getNumRows(unsigned int comp) const
    {
      return self().getNumRows(comp);
    }
    size_type getNumCols(unsigned int comp) const
    {
      return self().getNumCols(comp);
    }

    /// return number of components/blocks
    unsigned int getNumComponents() const
    {
      return self().getNumComponents();
    }

    size_type row(unsigned int block_r, unsigned int block_c, size_type r)
    {
      setRow(block_r);
      setCol(block_c);
      return row(r);
    }

    size_type col(unsigned int block_r, unsigned int block_c, size_type c)
    {
      setRow(block_r);
      setCol(block_c);
      return col(c);
    }

    Derived& self()
    {
      return static_cast<Derived&>(*this);
    }
    const Derived& self() const
    {
      return static_cast<const Derived&>(*this);
    }
  };


  /// Mapper implementation for non-parallel block matrices
  struct BlockMapper : public MapperBase<BlockMapper>
  {
    using Super = MapperBase<BlockMapper>;
    using size_type = unsigned long;

    /// Default constructor.
    BlockMapper()
      : nComp(0), rowOffset(0), colOffset(0), nrow(0), ncol(0), sizes(0)
    {}

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

    /// Copy constructor.
    BlockMapper(BlockMapper const& other);

    /// Copy assignment operator.
    BlockMapper& operator=(BlockMapper const& other);

    /// Constructor that takes a system-matrix
    explicit BlockMapper(Matrix<DOFMatrix*> const& orMat);

    /// Constructor that takes a solvermatrix,
    /// calls constructor for system-matrix
    explicit BlockMapper(SolverMatrix<Matrix<DOFMatrix*>> const& sm)
      : BlockMapper(*sm.getOriginalMat())
    {}

    /// Constructor for single matrix.
    explicit BlockMapper(DOFMatrix const* mat);

    /// calculate row offset for row component \param r
    void setRow(unsigned int r)
    {
      assert( r <= sizes.size() );
      rowOffset = sum(r);
    }

    /// calculate column offset for col component \param c
    void setCol( unsigned int c)
    {
      assert( c <= sizes.size() );
      colOffset = sum(c);
    }

    size_type row(size_type r) const
    {
      return r + rowOffset;
    }
    size_type col(size_type c) const
    {
      return c + colOffset;
    }

    size_type getNumRows() const
    {
      return nrow;
    }
    size_type getNumCols() const
    {
      return ncol;
    }

    size_type getNumRows(unsigned int comp) const
    {
      return sizes[comp];
    }
    size_type getNumCols(unsigned int comp) const
    {
      return sizes[comp];
    }

    unsigned int getNumComponents() const
    {
      return nComp;
    }

  private: // methods

    ///compute the sum of sizes from [0, end)
    size_type sum(unsigned int end) const
    {
      unsigned int ret = 0;
      for (unsigned int i = 0; i < end; ++i)
        ret += sizes[i];
      return ret;
    }

  private: // variables
    unsigned int nComp;
    size_type rowOffset;
    size_type colOffset;
    unsigned int nrow;
    unsigned int ncol;

    std::vector<size_type> sizes;
  };




  /// Mapper implementation for non-parallel rectangular block matrices
  struct RectangularMapper : public MapperBase<RectangularMapper>
  {
    using Super = MapperBase<RectangularMapper>;
    using size_type = unsigned long;

    /// Constructor for block-matrices
    explicit RectangularMapper(Matrix<DOFMatrix*> const& sm);

    /// Constructor for block-matrices, calls constructor for system-matrix
    explicit RectangularMapper(SolverMatrix<Matrix<DOFMatrix*>> const& sm)
      : RectangularMapper(*sm.getOriginalMat())
    {}

    /// calculate row offset for row component \param r
    void setRow( unsigned int r)
    {
      assert( r <= sizes_rows.size() );
      rowOffset = sum(r, sizes_rows);
    }

    /// calculate column offset for col component \param c
    void setCol( unsigned int c)
    {
      assert( c <= sizes_cols.size() );
      colOffset = sum(c, sizes_cols);
    }

    size_type row(size_type r) const
    {
      return r + rowOffset;
    }
    size_type col(size_type c) const
    {
      return c + colOffset;
    }

    size_type getNumRows() const
    {
      return nrow;
    }
    size_type getNumCols() const
    {
      return ncol;
    }

    size_type getNumRows(unsigned int comp) const
    {
      return sizes_rows[comp];
    }
    size_type getNumCols(unsigned int comp) const
    {
      return sizes_cols[comp];
    }

    unsigned int getNumComponents() const
    {
      return nRowComp;
    }
    unsigned int getNumRowComponents() const
    {
      return nRowComp;
    }
    unsigned int getNumColComponents() const
    {
      return nColComp;
    }

  private: // methods

    ///compute the sum of sizes from [0, end)
    size_type sum(unsigned int end, std::vector<size_type>& sizes) const
    {
      unsigned int ret = 0;
      for (unsigned int i = 0; i < end; ++i)
        ret += sizes[i];
      return ret;
    }

  private: // variables

    unsigned int nRowComp, nColComp;
    size_type rowOffset;
    size_type colOffset;
    unsigned int nrow;
    unsigned int ncol;

    std::vector<size_type> sizes_rows;
    std::vector<size_type> sizes_cols;
  };

} // end namespace AMDiS
