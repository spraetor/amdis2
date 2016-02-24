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


#include <algorithm>
#include "Global.h"
#include "DOFMatrix.h"
#include "parallel/MatrixNnzStructure.hpp"
#include "parallel/ParallelDofMapping.hpp"

using namespace std;

namespace AMDiS
{
  namespace Parallel
  {

    MatrixNnzStructure::~MatrixNnzStructure()
    {
      clear();
    }


    void MatrixNnzStructure::clear()
    {
      if (dnnz)
      {
        delete [] dnnz;
        dnnz = NULL;
      }

      if (onnz)
      {
        delete [] onnz;
        onnz = NULL;
      }
    }

    void MatrixNnzStructure::create(Matrix<DOFMatrix*>& mat,
                                    ParallelDofMapping& rowDofMap,
                                    ParallelDofMapping& colDofMap,
                                    PeriodicMap* perMap,
                                    const ElementObjectDatabase& elObjDb,
                                    bool localMatrix)
    {
      FUNCNAME("MatrixNnzStructure::create()");

      int nRankRows = rowDofMap.getRankDofs();	// Number of DOFs owned by rank.

      int rankStartRowIndex = rowDofMap.getStartDofs();		// Smallest global index of a DOF owned by the rank.

      int nRankCols = colDofMap.getRankDofs();

      int nOverallCols = colDofMap.getOverallDofs();	// Number of global DOFs (this value is thus the same on all ranks).

      int rankStartColIndex = colDofMap.getStartDofs();

      create(nRankRows, (!localMatrix ? nRankRows : -1));

      using mtl::tag::row;
      using mtl::tag::nz;
      using mtl::begin;
      using mtl::end;
      namespace traits = mtl::traits;
      typedef DOFMatrix::base_matrix_type Matrix;
      typedef traits::range_generator<row, Matrix>::type cursor_type;
      typedef traits::range_generator<nz, cursor_type>::type icursor_type;

      typedef vector<pair<int, int>> MatrixNnzEntry;
      //     typedef map<int, DofContainer> RankToDofContainer;

      // Stores to each rank a list of nnz entries (i.e. pairs of row and column
      // index) that this rank will send to. These nnz entries will be assembled
      // on this rank, but because the row DOFs are not DOFs of this rank they
      // will be send to the owner of the row DOFs.
      map<int, MatrixNnzEntry> sendMatrixEntry;

      // Maps to each DOF that must be send to another rank the rank number of the
      // receiving rank.
      map<pair<DegreeOfFreedom, int>, int> sendDofToRank;

      int nComponents = mat.getNumRows();

      // First, create for all ranks, to which we send data to, MatrixNnzEntry
      // object with 0 entries.
      for (int i = 0; i < nComponents; i++)
      {
        const FiniteElemSpace* feSpace = NULL;
        for (int j = 0; j < nComponents; j++)
          if (mat[i][j])
            feSpace = mat[i][j]->getRowFeSpace();

        TEST_EXIT_DBG(feSpace)("No FE space found!\n");

        for (DofComm::Iterator it(rowDofMap.getDofComm(feSpace).getRecvDofs(), feSpace);
             !it.end(); it.nextRank())
        {
          sendMatrixEntry[it.getRank()].resize(0);

          for (; !it.endDofIter(); it.nextDof())
            sendDofToRank[make_pair(it.getDofIndex(), i)] = it.getRank();
        }
      }

      // Create list of ranks from which we receive data from.
      std::set<int> recvFromRank;
      for (int i = 0; i < nComponents; i++)
      {
        const FiniteElemSpace* feSpace = NULL;
        for (int j = 0; j < nComponents; j++)
          if (mat[i][j])
            feSpace = mat[i][j]->getRowFeSpace();

        for (DofComm::Iterator it(rowDofMap.getDofComm(feSpace).getSendDofs(), feSpace);
             !it.end(); it.nextRank())
          recvFromRank.insert(it.getRank());
      }


      // === Traverse matrices to create nnz data. ===

      for (int rowComp = 0; rowComp < nComponents; rowComp++)
      {
        for (int colComp = 0; colComp < nComponents; colComp++)
        {
          DOFMatrix* dofMat = mat[rowComp][colComp];

          if (!dofMat)
            continue;

          const FiniteElemSpace* rowFeSpace = dofMat->getRowFeSpace();
          const FiniteElemSpace* colFeSpace = dofMat->getColFeSpace();

          if (rowDofMap.isDefinedFor(rowComp) == false ||
              colDofMap.isDefinedFor(colComp) == false)
            continue;

          // === Prepare MTL4 iterators for the current DOF matrix. ===

          Matrix baseMat = dofMat->getBaseMatrix(); // TBD: warum hier keine Referenz?

          traits::col<Matrix>::type getCol(baseMat);
          traits::row<Matrix>::type getRow(baseMat);

          traits::const_value<Matrix>::type value(baseMat);
          cursor_type cursor = begin<row>(baseMat);
          cursor_type cend = end<row>(baseMat);


          // === Iterate on all DOFs of the row mapping. ===

          DofMap::iterator rowIt = rowDofMap[rowComp].begin(); // row dofmap DOFIt
          DofMap::iterator rowEndIt = rowDofMap[rowComp].end();
          for (; rowIt != rowEndIt; ++rowIt)
          {

            // Go to the corresponding matrix row (note, both the mapping and the
            // matrix rows are stored in increasing order).
            while (cursor.value() != rowIt->first)
              ++cursor;
            size_t _row = cursor.value(); // DegreeOfFreedom

            // The corresponding global matrix row index of the current row DOF.
            int petscRowIdx = 0;
            if (localMatrix)
            {
              petscRowIdx = rowDofMap.getLocalMatIndex(rowComp, _row);
            }
            else
            {
              if (rowDofMap.isMatIndexFromGlobal())
                petscRowIdx = rowDofMap.getMatIndex(rowComp, rowIt->second.global);
              else
                petscRowIdx = rowDofMap.getMatIndex(rowComp, _row);
            }


            // Set of periodic row associations (is empty, if row DOF is not
            // periodic.
            std::set<int> perRowAsc;
            if (perMap)
              perMap->fillAssociations(rowFeSpace, rowIt->second.global, elObjDb, perRowAsc);


            if (localMatrix || rowDofMap[rowComp].isRankDof(_row))
            {
              // === The current row DOF is a rank DOF, so create the       ===
              // === corresponding nnz values directly on rank's nnz data.  ===

              // This is the local row index of the local PETSc matrix.
              int localPetscRowIdx = petscRowIdx;

              if (localMatrix == false)
                localPetscRowIdx -= rankStartRowIndex;

              TEST_EXIT(localPetscRowIdx >= 0 && localPetscRowIdx < nRankRows)
              ("Should not happen! \n Debug info: DOF = %d   globalRowIndx = %d   petscRowIdx = %d   localPetscRowIdx = %d   rStart = %d   compontens = %d from %d   nRankRows = %d\n",
               _row,
               rowDofMap[rowComp][_row].global,
               petscRowIdx,
               localPetscRowIdx,
               rankStartRowIndex,
               rowComp,
               nComponents,
               nRankRows);


              if (localMatrix)
              {
                for (icursor_type icursor = begin<nz>(cursor),
                     icend = end<nz>(cursor); icursor != icend; ++icursor)
                  if (colDofMap[colComp].isSet(getCol(*icursor)))
                    dnnz[localPetscRowIdx]++;
              }
              else
              {
                MultiIndex colDofIndex;

                // Traverse all non zero entries in this row.
                for (icursor_type icursor = begin<nz>(cursor),
                     icend = end<nz>(cursor); icursor != icend; ++icursor)
                {
                  if (colDofMap[colComp].find(getCol(*icursor), colDofIndex) == false)
                    continue;

                  // Set of periodic row associations (is empty, if row DOF is not
                  // periodic.
                  std::set<int> perColAsc = perRowAsc;
                  if (perMap)
                    perMap->fillAssociations(colFeSpace, colDofIndex.global, elObjDb, perColAsc);

                  if (perColAsc.empty())
                  {
                    if (colDofMap[colComp].isRankDof(getCol(*icursor)))
                      dnnz[localPetscRowIdx]++;
                    else
                      onnz[localPetscRowIdx]++;
                  }
                  else
                  {
                    vector<int> newCols;
                    perMap->mapDof(colFeSpace, colDofIndex.global, perColAsc, newCols);

                    for (size_t aa = 0; aa < newCols.size(); aa++)
                    {
                      int petscColIdx = colDofMap.getMatIndex(colComp, newCols[aa]);

                      // The row DOF is a rank DOF, if also the column is a rank DOF,
                      // increment the diagNnz values for this row,
                      // otherwise the offdiagNnz value.
                      if (petscColIdx >= rankStartColIndex &&
                          petscColIdx < rankStartColIndex + nRankCols)
                        dnnz[localPetscRowIdx]++;
                      else
                        onnz[localPetscRowIdx]++;
                    }
                  }
                }

                if (!perRowAsc.empty())
                {
                  vector<int> newRows;
                  perMap->mapDof(rowFeSpace, rowIt->second.global, perRowAsc, newRows);

                  dnnz[localPetscRowIdx] +=
                    (newRows.size() - 1) * (onnz[localPetscRowIdx] + dnnz[localPetscRowIdx]);
                }
              }
            }
            else
            {
              // === The current row DOF is not a rank DOF, i.e., its values   ===
              // === are also created on this rank, but afterthere they will   ===
              // === be send to another rank. So we need to send also the      ===
              // === corresponding nnz structure of this row to the corres-    ===
              // === ponding rank.                                             ===

              // Send all non zero entries to the member of the row DOF.
              int sendToRank = sendDofToRank[make_pair(cursor.value(), rowComp)];
              MultiIndex colDofIndex;

              for (icursor_type icursor = begin<nz>(cursor),
                   icend = end<nz>(cursor); icursor != icend; ++icursor)
              {
                if (colDofMap[colComp].find(getCol(*icursor), colDofIndex) == false)
                  continue;

                int petscColIdx = (colDofMap.isMatIndexFromGlobal() ?
                                   colDofMap.getMatIndex(colComp, colDofIndex.global) :
                                   colDofMap.getMatIndex(colComp, getCol(*icursor)));

                sendMatrixEntry[sendToRank].
                push_back(make_pair(petscRowIdx, petscColIdx));
              }

            } // if (isRankDof[cursor.value()]) ... else ...
          } // for each row in mat[i][j]
        }
      }


      if (localMatrix == false)
      {
        // === Send and recv the nnz row structure to/from other ranks. ===

        StdMpi<MatrixNnzEntry> stdMpi(rowDofMap.getMpiComm());
        stdMpi.send(sendMatrixEntry);

        for (std::set<int>::iterator it = recvFromRank.begin();
             it != recvFromRank.end(); ++it)
          stdMpi.recv(*it);
        stdMpi.startCommunication();

        // === Evaluate the nnz structure this rank got from other ranks and add ===
        // === it to the PETSc nnz data structure.                               ===

        for (map<int, MatrixNnzEntry>::iterator it = stdMpi.getRecvData().begin();
             it != stdMpi.getRecvData().end(); ++it)
        {
          if (it->second.size() > 0)
          {
            for (unsigned int i = 0; i < it->second.size(); i++)
            {
              int r = it->second[i].first;
              int c = it->second[i].second;

              int localRowIdx = r - rankStartRowIndex;

              TEST_EXIT_DBG(localRowIdx >= 0 && localRowIdx < nRankRows)
              ("Got row index %d/%d (nRankRows = %d) from rank %d. Should not happen!\n",
               r, localRowIdx, nRankRows, it->first);

              if (c < rankStartColIndex || c >= rankStartColIndex + nRankCols)
                onnz[localRowIdx]++;
              else
                dnnz[localRowIdx]++;
            }
          }
        }
      }

      // The above algorithm for calculating the number of nnz per row over-
      // approximates the value, i.e., the number is always equal or larger to
      // the real number of nnz values in the global parallel matrix. For small
      // matrices, the problem may arise, that the result is larger than the
      // number of elements in a row. This is fixed in the following.

      for (int i = 0; i < nRankRows; i++)
      {
        dnnz[i] = std::min(dnnz[i], nRankCols);
        if (onnz)
          onnz[i] = std::min(onnz[i], nOverallCols - nRankCols);
      }

#if (DEBUG != 0)
      int nMax = 0;
      int nSum = 0;
      for (int i = 0; i < nRankRows; i++)
      {
        nMax = std::max(nMax, dnnz[i]);
        nSum += dnnz[i];
      }

      MSG("NNZ in diag block: max = %d, avrg = %.0f\n",
          nMax, (nSum > 0 ? (static_cast<double>(nSum) / nRankRows) : 0));

      if (!localMatrix)
      {
        nMax = 0;
        nSum = 0;

        for (int i = 0; i < nRankRows; i++)
        {
          nMax = std::max(nMax, onnz[i]);
          nSum += onnz[i];
        }

        MSG("NNZ in offdiag block: max = %d, avrg = %.0f\n",
            nMax, (nSum > 0 ? (static_cast<double>(nSum) / nRankRows) : 0));
      }
#endif
    }


    void MatrixNnzStructure::create(int nRows0, int nRows1)
    {
      if (nRows0 == 0)
        return;

      TEST_EXIT_DBG(nRows0 > 0)("Should not happen!\n");

      clear();

      dnnz = new int[nRows0];
      for (int i = 0; i < nRows0; i++)
        dnnz[i] = 0;

      if (nRows1 > 0)
      {
        onnz = new int[nRows1];
        for (int i = 0; i < nRows1; i++)
          onnz[i] = 0;
      }
    }

  }
}
