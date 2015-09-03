/** \file Mapper.cc */

#include "Mapper.h"
#include "FiniteElemSpace.h"
#include "DOFMatrix.h"
#include "FixVec.h"

namespace AMDiS
{

  BlockMapper::BlockMapper(BlockMapper const& other)
    : nComp(other.nComp),
      rowOffset(other.rowOffset),
      colOffset(other.colOffset),
      nrow(other.nrow),
      ncol(other.ncol)
  {
    sizes.resize(nComp);
    std::copy(other.sizes.begin(), other.sizes.end(), sizes.begin());
  }


  BlockMapper& BlockMapper::operator=(BlockMapper const& other)
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


  BlockMapper::BlockMapper(const Matrix<DOFMatrix*>& orMat )
    : nComp(orMat.getSize()),
      rowOffset(0), colOffset(0), nrow(0), ncol(0), sizes(nComp)
  {
    const int ns = orMat.getNumRows();
    for (int i= 0; i < ns; i++)
    {
      sizes[i] = orMat[i][i]->getFeSpace()->getAdmin()->getUsedSize();
      nrow += sizes[i];
    }
    ncol = nrow;
  }


  BlockMapper::BlockMapper(const DOFMatrix* sm )
    : nComp(1),
      rowOffset(0), colOffset(0),
      nrow(0), ncol(0),
      sizes(nComp)
  {
    sizes[0] = sm->getFeSpace()->getAdmin()->getUsedSize();
    nrow += sizes[0];
    ncol = nrow;
  }


  RectangularMapper::RectangularMapper(const Matrix<DOFMatrix*>& orMat)
    : nRowComp(orMat.getNumRows()),
      nColComp(orMat.getNumCols()),
      rowOffset(0), colOffset(0), nrow(0), ncol(0),
      sizes_rows(nRowComp), sizes_cols(nColComp)
  {
    for (unsigned int i = 0; i < nRowComp; i++)
    {
      for (unsigned int j = 0; j < nColComp; j++)
      {
        if (orMat[i][j])
        {
          sizes_rows[i] = orMat[i][j]->getRowFeSpace()->getAdmin()->getUsedSize();
          sizes_cols[j] = orMat[i][j]->getColFeSpace()->getAdmin()->getUsedSize();
        }
      }
      nrow += sizes_rows[i];
    }
    for (unsigned int j = 0; j < nColComp; j++)
      ncol += sizes_cols[j];
  }

} // end namespace AMDiS
