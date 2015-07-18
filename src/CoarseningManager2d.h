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



/** \file CoarseningManager2d.h */

#ifndef AMDIS_COARSENINGMANAGER_2D_H
#define AMDIS_COARSENINGMANAGER_2D_H

#include "CoarseningManager.h"

namespace AMDiS {

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

    /// destructor
    virtual ~CoarseningManager2d() {}

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


}

#endif // AMDIS_COARSENINGMANAGER_2D_H
