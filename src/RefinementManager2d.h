/** \file RefinementManager2d.h */

#pragma once

namespace AMDiS
{
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

  protected:
    /// Used by \ref setNewCoords
    void newCoordsFct(ElInfo* elInfo);

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
    void getRefinePatch(ElInfo** elInfo, DegreeOfFreedom* edge[2], int dir,
                        RCNeighbourList& refineList, int* n_neigh);

    /// Refines all elements in the patch.
    DegreeOfFreedom refinePatch(DegreeOfFreedom* edge[2],
                                RCNeighbourList& refineList,
                                int n_neigh, bool bound);

    /// Implements RefinementManager::refineFunction.
    ElInfo* refineFunction(ElInfo* elInfo);

    /// Refines one Triangle.
    void bisectTriangle(Triangle* el, DegreeOfFreedom* newDofs[3]);
  };

} // end namespace AMDiS
