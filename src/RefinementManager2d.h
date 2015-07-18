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



/** \file RefinementManager2d.h */

#ifndef AMDIS_REFINEMENT_MANAGER_2D_H
#define AMDIS_REFINEMENT_MANAGER_2D_H

namespace AMDiS {

  /** \ingroup Adaption
   * \brief
   * Implements a RefinementManager for 2-dimensional meshes.
   */
  class RefinementManager2d : public RefinementManager
  {
  public:
    /// Calls base class constructor.
    RefinementManager2d()
      : RefinementManager()
    {}

    /// Destructor 
    virtual ~RefinementManager2d() {}

  protected:
    /// Used by \ref setNewCoords
    void newCoordsFct(ElInfo *elInfo);

    /// Implements RefinementManager::setNewCoords
    void setNewCoords(int macroEl = -1);

    /** \brief
     *  gets the elements around the refinement edge with vertices  node[0] and
     *  node[1] ; refines those elements at this edge are not compatible 
     *  devisible;  			                                 
     *  dir determines the direction of the first neighbour
     *  get_refine_patch returns  1  if the domain's boundary is reached while
     *  looping around the refinement edge, otherwise  0
     */
    void getRefinePatch(ElInfo **elInfo, DegreeOfFreedom *edge[2], int dir,
			RCNeighbourList &refineList, int *n_neigh);

    /// Refines all elements in the patch.
    DegreeOfFreedom refinePatch(DegreeOfFreedom *edge[2], 
				RCNeighbourList &refineList,
				int n_neigh, bool bound);

    /// Implements RefinementManager::refineFunction.
    ElInfo* refineFunction(ElInfo* elInfo);

    /// Refines one Triangle.
    void bisectTriangle(Triangle *el, DegreeOfFreedom* newDofs[3]);
  };

}

#endif // AMDIS_REFINEMENT_MANAGER_2D_H
