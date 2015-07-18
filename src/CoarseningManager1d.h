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



/** \file CoarseningManager1d.h */

#ifndef AMDIS_COARSENINGMANAGER_1D_H
#define AMDIS_COARSENINGMANAGER_1D_H

#include "CoarseningManager.h"

namespace AMDiS {

  /** \ingroup Adaption 
   * \brief
   * Implements a CoarseningManager for 1-dimensional meshes.
   */
  class CoarseningManager1d : public CoarseningManager
  {
  public:
    /// Calls base class constructor and checks dimension of mesh. 
    CoarseningManager1d() 
      : CoarseningManager() 
    {}

    /// destructor
    virtual ~CoarseningManager1d() {}

    /** \brief
     * Overloads CoarseningManager::coarsenMesh. In 1d a simple recursive
     * coarsening algorithm is implemented which doesn't need coarsenFunction.
     */
    Flag coarsenMesh(Mesh *aMesh);

  protected:
    /// Not needed in this sub class
    void coarsenFunction(ElInfo *) 
    {
      FUNCNAME("CoarseningManager1d::coarsenFunction");
      ERROR_EXIT("not used for dim = 1");
    }
  
    /// Needed instead of coarsenFunction in 1d.
    int coarsenRecursive(Line *parent);
  
  };

}

#endif // AMDIS_COARSENINGMANAGER_1D_H
