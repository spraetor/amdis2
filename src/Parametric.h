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



/** \file Parametric.h */

#ifndef AMDIS_PARAMETRIC_H
#define AMDIS_PARAMETRIC_H

#include "AMDiS_fwd.h"

namespace AMDiS {

  /** \brief
   * Interface for parametric elements. A Parametric object can be given to
   * a Mesh. 
   * Mesh traversal then adds parametric info automatically.
   */
  class Parametric {
  public:
    virtual ~Parametric() {}

    /** \brief
     * Can return a pointer to a new parametric ElInfo or the original pointer 
     * given with elinfo. If the returned ElInfo is newly created, 
     * it must be deallocated in removeParametricInfo.
     * If the element is not parametric, elinfo is returned with unchanged 
     * content:
     * addParametricInfo(elinfo) == elinfo and 
     * *addParametricInfo(elInfo) == *elinfo.
     * If parametric information are added to the returned elinfo 
     * ElInfo::parametric must be set true.
     */
    virtual ElInfo *addParametricInfo(ElInfo *elInfo) = 0;

    /** \brief
     * Must free memory which was allocated in addParametricInfo() for this 
     * ElInfo. Must set ElInfo::parametric to false if the result is no longer
     * parametric.
     * Must return the corresponding original ElInfo pointer given to 
     * addParametricInfo(). 
     * removeParametricInfo(addParametricInfo(elinfo)) == elinfo !!!
     * It is not necessary to reset any data on the original ElInfo to their
     * original values in any case.
     * *(removeParametricInfo(addParametricInfo(elinfo))) != *elinfo is ok!!!
     */
    virtual ElInfo *removeParametricInfo(ElInfo *elInfo) = 0;
  };

  /// Implementation of linear parametric elements.
  class ParametricFirstOrder : public Parametric
  {
  public:
    /// Constructor. \ref dofCoords are set to coords.
    ParametricFirstOrder(WorldVector<DOFVector<double>*> *coords)
      : dofCoords_(coords)
    {}

    /// Implementation of \ref Parametric::addParametricInfo().
    ElInfo *addParametricInfo(ElInfo *elInfo);

    /// Implementation of \ref Parametric::removeParametricInfo().
    ElInfo *removeParametricInfo(ElInfo *elInfo);

  protected:
    /// Pointer to a DOFVector of world coordinates.
    WorldVector<DOFVector<double>*> *dofCoords_;
  };

}

#endif
