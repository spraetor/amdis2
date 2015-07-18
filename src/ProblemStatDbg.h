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



/** \file ProblemStatDbg.h */

#ifndef AMDIS_PROBLEM_STAT_DBG_H
#define AMDIS_PROBLEM_STAT_DBG_H

#include "ProblemStat.h"

namespace AMDiS {

  class ProblemStatDbg : public ProblemStatSeq
  {
  public:
    ProblemStatDbg(std::string nameStr,
		   ProblemIterationInterface *problemIteration = NULL)
      : ProblemStatSeq(nameStr, problemIteration)
    {}

    /** \brief
     * This function writes the system matrix and the system vector to a file
     * for debugging. The entries of both, the matrix and the vector, are not
     * indicated by the dof indices, but by the world coordinates of the dofs.
     * This makes it possible to compare dof matrices, where the dof indices
     * are different, e.g., when using domain decomposition parallelization.
     */
    void writeDbgMatrix(std::string filename);

    /** \brief
     * Reads a file, that was created using the function \ref writeDbgMatrix, and
     * compares the date with the system matrix of this problem.
     */
    void readAndCompareDbgMatrix(std::vector<std::string> filenames);

  protected:
    /** \brief
     * Create from the current problem a map from dof indices to world coordinates.
     * This function is used by the debugging function \ref writeDbgMatrix and \ref
     * readAndCompareDbgMatrix.
     */
    void createDofToCoordMap(DofToCoord &dofMap);

    /** \brief
     * Create from the current problem a map from world coordinates of all dofs to
     * the corresponding dof indices. This function is used by the debugging function 
     * \ref writeDbgMatrix and \ref readAndCompareDbgMatrix.
     */
    void createCoordToDofMap(CoordToDof &dofMap);
  };
}

#endif
