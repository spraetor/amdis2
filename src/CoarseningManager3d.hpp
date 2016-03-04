#pragma once

// AMDiS includes
#include "AMDiS_base.hpp"
#include "CoarseningManager.hpp"

namespace AMDiS
{
  /** \ingroup Adaption
   * \brief
   * Implements a CoarseningManager for 3-dimensional meshes.
   */
  class CoarseningManager3d : public CoarseningManager
  {
  public:
    /// Calls base class constructor and checks dimension of mesh.
    CoarseningManager3d()
      : CoarseningManager()
    {}

  protected:
    /// Implements \ref CoarseningManager::coarsenFunction
    void coarsenFunction(ElInfo* el_info);

    /// Coarsens a single Tetrahedron of the coarsening patch. DOFs
    /// in the interior of the element are removed; DOFs for higher order
    /// at the boundary or the coarsening patch still belong to
    /// the parent. Do not remove them form the mesh!!!
    void coarsenTetrahedron(RCNeighbourList& coarsenList, int index);

    /// Gets the patch for coarsening starting on element
    /// el_info->el in direction of neighbour [3-dir]; returns 1 if a boundary
    /// reached and 0 if we come back to the starting element.
    /// We complete the loop also in the case of a incompatible
    /// coarsening patch since then all marks of patch elements are reset by
    /// coarsenPatch() and this minimizes calls of traverseNeighbour();
    /// if we reach a boundary while looping around the edge we loop back to
    /// the starting element before we return
    bool getCoarsenPatch(ElInfo* el_info, DegreeOfFreedom* edge[2],
                         int dir, RCNeighbourList& coarsenList, int* n_neigh);

    /// First rebuild the DOFs on the parents then do restriction
    /// of data (if possible) and finally coarsen the patch elements
    void coarsenPatch(RCNeighbourList& coarsenList, int n_neigh, int bound);

  };

} // end namespace AMDiS
