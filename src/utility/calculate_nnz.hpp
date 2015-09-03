/** \file calculate_nnz.hpp */

#pragma once

namespace AMDiS
{
  inline void calculate_nnz(const DOFMatrix& matrix, std::vector<size_t>& nnz_per_row)
  {
    std::vector<std::set<DegreeOfFreefom>> nnz(nnz_per_row.size());
    std::vector<DegreeOfFreedom> rowIndices, colIndices;
    const FiniteElemSpace* rowFeSpace = matrix.getRowFeSpace(),
                           colFeSpace = matrix.getColFeSpace();

    TraverseStack stack;
    ElInfo* elInfo = stack.traverseFirst(rowFeSpace->getMesh(), -1, Mesh::CALL_LEAF_EL);
    while (elInfo)
    {
      rowFeSpace->getBasisFcts()->getLocalIndices(elInfo->getElement(),
          rowFeSpace->getAdmin(),
          rowIndices);
      if (rowFeSpace == colFeSpace)
      {
        colIndices = rowIndices;
      }
      else
      {
        colFeSpace->getBasisFcts()->getLocalIndices(elInfo->getElement(),
            colFeSpace->getAdmin(),
            colIndices);
      }

      for (size_t r = 0; r < rowIndices.size(); r++)
        for (size_t c = 0; c < rowIndices.size(); c++)
          nnz[rowIndices[r]].insert(colIndices[c]);

      elInfo = stack.traverseNext(elInfo);
    }

    for (size_t r = 0; r < nnz.size(); r++)
      nnz_per_row[r] = nnz[r].size();
  }

} // end namspace AMDiS
