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



/** \file calculate_nnz.hpp */

#ifndef AMDIS_UTILITY_CALCULATE_NNZ_HPP
#define AMDIS_UTILITY_CALCULATE_NNZ_HPP

namespace AMDiS {

  inline void calculate_nnz(const DOFMatrix& matrix, std::vector<size_t>& nnz_per_row)
  {
    std::vector<std::set<DegreeOfFreefom> > nnz(nnz_per_row.size());
    std::vector<DegreeOfFreedom> rowIndices, colIndices;
    const FiniteElemSpace* rowFeSpace = matrix.getRowFeSpace(), 
			   colFeSpace = matrix.getColFeSpace();
      
    TraverseStack stack;
    ElInfo *elInfo = stack.traverseFirst(rowFeSpace->getMesh(), -1, Mesh::CALL_LEAF_EL);
    while (elInfo) {
      rowFeSpace->getBasisFcts()->getLocalIndices(elInfo->getElement(),
						  rowFeSpace->getAdmin(),
						  rowIndices);
      if (rowFeSpace == colFeSpace) {
	colIndices = rowIndices;
      } else {
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

#endif