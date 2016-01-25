/** \file RefinementManager3d.h */

#pragma once

//#include "RefinementManager.h"

namespace AMDiS
{
  /** \ingroup Adaption
   * \brief
   * Implements a RefinementManager for 3-dimensional meshes.
   */
  class RefinementManager3d : public RefinementManager
  {
  public:
    /// Calls base class constructor.
    RefinementManager3d() = default;

  protected:
    /// Used by \ref setNewCoords
    void newCoordsFct(ElInfo* el_info, RCNeighbourList& refineList);

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
    bool getRefinePatch(ElInfo** el_info,
                        DegreeOfFreedom* edge[2],
                        int direction,
                        RCNeighbourList& refineList,
                        int* n_neigh);

    /// Refines all elements in the patch.
    DegreeOfFreedom refinePatch(DegreeOfFreedom* edge[2], RCNeighbourList& refineList,
                                int n_neigh, bool bound);

    /// Implements RefinementManager::refineFunction.
    ElInfo* refineFunction(ElInfo* el_info);

    /// Refines one Tetrahedron.
    void bisectTetrahedron(RCNeighbourList& refineList, int index,
                           DegreeOfFreedom* dof[3], DegreeOfFreedom* edge[2]);

    /// Used by \ref bisectTetrahedron
    void fillPatchConnectivity(RCNeighbourList& refineList, int index);
  };


  struct FixRefinementPatch
  {
    using EdgeInEl       = std::pair<Element*, int>;
    using EdgesMap       = std::pair<EdgeInEl, EdgeInEl>;
    using ConnectedEdges = std::vector<EdgesMap>;

    static std::map<Mesh*, ConnectedEdges> connectedEdges;

    static void getOtherEl(Mesh* mesh,
                           TraverseStack* stack,
                           std::vector<EdgeInEl>& refineEdges);
  };

} // end namespace AMDiS
