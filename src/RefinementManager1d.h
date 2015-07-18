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



/** \file RefinementManager1d.h */

#ifndef AMDIS_REFINEMENT_MANAGER_1D_H
#define AMDIS_REFINEMENT_MANAGER_1D_H

namespace AMDiS {

  /** \ingroup Adaption
   * \brief
   * Implements a RefinementManager for 1-dimensional meshes.
   */
  class RefinementManager1d : public RefinementManager
  {
  public:
    /// Calls base class constructor.
    RefinementManager1d()
      : RefinementManager()
    {}

    /// destructor 
    virtual ~RefinementManager1d() {}

    /// Implements RefinementManager::refineMesh.
    Flag refineMesh(Mesh *aMesh);

    /// Implements RefinementManager::setNewCoords
    void setNewCoords(int macroEl = -1);

  protected:
    /// Used by refineMesh while mesh traversal
    void recursiveRefineFunction(ElInfo* el_info);

    /// Used by \ref setNewCoords
    void newCoordsFct(ElInfo *el_info);
  };

}

#endif // AMDIS_REFINEMENT_MANAGER_1D_H
