/** \file CoarseningManager2d.h */

#pragma once

#include "CoarseningManager.h"

namespace AMDiS 
{
  /** \ingroup Adaption  
   * \brief
   * Implements a CoarseningManager for 2-dimensional meshes.
   */
  class CoarseningManager2d : public CoarseningManager
  {
  public:

    /// Calls base class constructor and checks dimension of mesh. 
    CoarseningManager2d() 
      : CoarseningManager() 
    {}

  protected:
    /// Implements \ref CoarseningManager::coarsenFunction
    void coarsenFunction(ElInfo *el_info);

    /** \brief
     * Coarsens a single Triangle of the coarsening patch. DOFs
     * in the interior of the element are removed; DOFs for higher order
     * at the boundary or the coarsening patch still belong to 
     * the parent. Do not remove them form the mesh!!!                  
     */
    void coarsenTriangle(Triangle *el);

    /** \brief
     * First rebuild the DOFs on the parents then do restriction
     * of data (if possible) and finally coarsen the patch elements
     */
    void coarsenPatch(RCNeighbourList &coarsenList, int n_neigh, int bound);
  };

} // end namespace AMDiS
