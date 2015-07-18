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



/** \file RefinementManager3d.h */

#ifndef AMDIS_REFINEMENT_MANAGER_3D_H
#define AMDIS_REFINEMENT_MANAGER_3D_H

namespace AMDiS {

  /** \ingroup Adaption 
   * \brief
   * Implements a RefinementManager for 3-dimensional meshes.
   */
  class RefinementManager3d : public RefinementManager
  {
  public:
    /// Calls base class constructor.
    RefinementManager3d() 
      : RefinementManager()
    {}

    /// destructor 
    virtual ~RefinementManager3d() {}

  protected:
    /// Used by \ref setNewCoords
    void newCoordsFct(ElInfo *el_info, RCNeighbourList &refineList);

    /// Implements RefinementManager::setNewCoords
    void setNewCoords(int macroEl = -1);

    /** \brief
     * Gets the elements around the refinement edge with vertices edge[0] and
     * edge[1]. Refines those elements at this edge that are not compatible 
     * devisible. The function returns 1 if the domain's boundary is reached
     * while looping around the refinement edge, otherwise 0.
     * 
     * \param[in] direction   Determines the direction of the first neighbour
     *  
     */
    bool getRefinePatch(ElInfo **el_info, 
			DegreeOfFreedom *edge[2], 
			int direction,
			RCNeighbourList &refineList, 
			int *n_neigh);

    /// Refines all elements in the patch.
    DegreeOfFreedom refinePatch(DegreeOfFreedom *edge[2], RCNeighbourList &refineList,
				int n_neigh, bool bound);

    /// Implements RefinementManager::refineFunction.
    ElInfo* refineFunction(ElInfo* el_info);

    /// Refines one Tetrahedron.
    void bisectTetrahedron(RCNeighbourList &refineList, int index,
			   DegreeOfFreedom *dof[3], DegreeOfFreedom *edge[2]);

    /// Used by \ref bisectTetrahedron
    void fillPatchConnectivity(RCNeighbourList &refineList, int index);
  };

  class FixRefinementPatch 
  {
  public:
    typedef std::pair<Element*, int> EdgeInEl;
    typedef std::pair<EdgeInEl, EdgeInEl> EdgesMap;
    typedef std::vector<EdgesMap> ConnectedEdges;

    static std::map<Mesh*, ConnectedEdges> connectedEdges;

    static void getOtherEl(Mesh* mesh,
                           TraverseStack *stack, 
			   std::vector<EdgeInEl> &refineEdges);
  };

}

#endif // AMDIS_REFINEMENT_MANAGER_3D_H
