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



/** \file BoundaryManager.h */

#ifndef AMDIS_BOUNDARYMANAGER_H
#define AMDIS_BOUNDARYMANAGER_H

#include <map>
#include <vector>

#include "AMDiS_fwd.h"
#include "Boundary.h"
#include "BoundaryCondition.h"

namespace AMDiS {

  typedef std::map<BoundaryType, BoundaryCondition*> BoundaryIndexMap;

  /**
   * \ingroup Assembler
   *
   * \brief
   * A BoundaryManager handles a set of boundary conditions and applies
   * this conditions to DOFVectorBase and DOFMatrix objects. Each DOFVectorBase
   * and each DOFMatrix has its own BoundaryManager.
   */
  class BoundaryManager
  {
  public:
    BoundaryManager(const FiniteElemSpace *feSpace);

    BoundaryManager(BoundaryManager &bm);

    ~BoundaryManager();

    /// Adds a local boundary condition to the list of managed conditions.
    void addBoundaryCondition(BoundaryCondition *localBC);

    void initMatrix(DOFMatrix *matrix);

    void exitMatrix(DOFMatrix *matrix);

    void initVector(DOFVectorBase<double> *vector);

    void exitVector(DOFVectorBase<double> *vector);

    /// Calls DOFVectorBase::fillBoundaryCondition() for each local boundary
    /// condition in \ref localBCs.
    void fillBoundaryConditions(ElInfo *elInfo, DOFVectorBase<double> *vec);

    /// Calls DOFMatrix::fillBoundaryCondition() for each local boundary condition
    /// in \ref localBCs.
    void fillBoundaryConditions(ElInfo *elInfo, DOFMatrix *mat);

    /// Calls BoundaryCondition::boundResidual() for each boundary condition in 
    /// \ref localBCs.
    double boundResidual(ElInfo *elInfo, DOFMatrix *matrix, 
			 const DOFVectorBase<double> *dv);

    BoundaryCondition *getBoundaryCondition(BoundaryType type) const
    {
      return localBCs[type];
    }

    const BoundaryIndexMap& getBoundaryConditionMap() const
    {
      return localBCs;
    }

    void setBoundaryConditionMap(const BoundaryIndexMap& bcs) 
    {
      localBCs = bcs;
    }
    
    /// Returns true, if there is at least one boundary object with the given
    /// boundary id, which implements a periodic boundary.
    static bool isBoundaryPeriodic(BoundaryType b)
    {
      for (auto const& boundary_map : globalBoundaryMap[b])
	if (boundary_map[i]->isPeriodic())
	  return true;

      return false;
    }

  protected:
    /// Map of managed local boundary conditions.
    BoundaryIndexMap localBCs;

    /// Temporary variable for functions fillBoundaryconditions.
    BoundaryType* localBound;

    /// Temporary variable for functions fillBoundaryconditions.
    std::vector<DegreeOfFreedom> dofVec;

    /// Stores the number of byte that were allocated in the constructor for
    /// each localBounds value. Is used to free the memory in the destructor.
    int allocatedMemoryLocalBounds;

    /** \brief
     * For every boundary id we store here all possible boundary object (although
     * it's not clear if it is meaningful to have different boundary conditions on
     * the same boundary id). 
     *
     * We have to use this global variable, because the mesh traverse interface 
     * does not provide more information about traversed boundaries at elements
     * than the boundary id.
     *
     * TODO: Change interface such that mesh traverse returns the boundary objects
     * directly and we can remove this global variable. The biggest problem will be
     * than serialization and deserialization of the mesh.
     */
    static std::map<BoundaryType, std::vector<BoundaryCondition*> > globalBoundaryMap;
  };

}

#endif
