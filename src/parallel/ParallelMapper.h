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
#ifndef PARALLEL_MAPPER_H
#define PARALLEL_MAPPER_H

#include "parallel/MeshDistributor.h"
#include "solver/Mapper.h"
#include "MTL4Types.h"

namespace AMDiS
{
  namespace Parallel 
  {
    template< typename size_type_ > 
    class ParallelMapper_base : public MapperBase< ParallelMapper_base<size_type_> >
    {
      typedef size_type_ size_type;
      
      ParallelDofMapping& map;
      
      /// the current component row in the problem system
      size_type rowComp;
      
      /// the current component column in the problem system
      size_type colComp;
      
      unsigned int nComponents;
      
      MultiIndex rowMultiIndex;
      MultiIndex colMultiIndex;
      
    public:
      
      ParallelMapper_base(ParallelDofMapping& map)
      : map(map),
	rowComp(0), 
	colComp(0), 
	nComponents(map.getNumberOfComponents()) 
      { FUNCNAME_DBG("ParallelMapper_base::_constructor()");
      
	TEST_EXIT_DBG(nComponents > 0)
	  ("The system must have at least one component\n");
      }
      
      inline ParallelDofMapping& getMap()
      {
	return map;
      }
      
      size_type row(size_type r_) const
      {
// 	if (map[rowComp].find(r_, rowMultiIndex) == false)
// 	  throw std::runtime_error("row-index not found!\n");
// 	int globalRowDof = rowMultiIndex.global;
// 	int rowIndex = map.getMatIndex(rowComp, globalRowDof);
// 	return rowIndex;
	return map.getMatIndex(rowComp, r_);
      }
      
      size_type col(size_type c_) const
      {
// 	if (map[colComp].find(c_, colMultiIndex) == false)
// 	  throw std::runtime_error("col-index not found!\n");
// 	int globalColDof = colMultiIndex.global;
// 	int colIndex = map.getMatIndex(colComp, globalColDof);
// 	return colIndex;
	return map.getMatIndex(colComp, c_);
      }
      
      inline unsigned int getNumComponents() const { return nComponents; }
      inline void setRow(size_type r_) { rowComp = r_ ;}
      inline void setCol(size_type c_) { colComp = c_; }
      
      inline size_type getNumRows() const { return map.getOverallDofs(); }
      inline size_type getNumCols() const { return map.getOverallDofs(); }
    };
  
  } // end namespace Parallel
  
  typedef Parallel::ParallelMapper_base< MTLTypes::size_type > ParallelMapper;
  
} // end namespace AMDiS

#endif
