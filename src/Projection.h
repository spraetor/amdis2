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



/** \file Projection.h */

#ifndef AMDIS_PROJECTION_H
#define AMDIS_PROJECTION_H

#include "FixVec.h"
#include <map>

namespace AMDiS {

  /// Different possible types for a \ref Projection.
  enum ProjectionType {
    BOUNDARY_PROJECTION = 0, /**< Projection of boundary parts of an element. */
    VOLUME_PROJECTION = 1    /**< Projection of whole elements. */
  };

  /** \brief
   * A Projection is a mapping from world coordinates to world coordinates.
   * It must fullfill the condition \ref project(project(x)) == project(x).
   * Each projection object has its own unique \ref projectionID. This is
   * used to connect the projection indices used in the macro file with the
   * corresponding projection object. See \ref getProjection(int id).
   */
  class Projection
  {
  public:
    /// Constructs a prjection with given id and type.
    Projection(int id, ProjectionType type) 
      : projectionID(id),
	projectionType(type)
    {
      TEST_EXIT(id != 0)("don't use 0 as projection id. is used as no projection\n");
      TEST_EXIT(projectionMap[id] == NULL)
	("there is already a projection with this id\n");
      projectionMap[id] = this;
    }

    virtual ~Projection() {}

    /// Projection method. Must be overriden in sub classes.
    virtual void project(WorldVector<double>& x) = 0;

    /// Returns \ref projectionID.
    inline int getID() 
    { 
      return projectionID; 
    }

    /// Returns \ref projectionType;
    inline ProjectionType getType() 
    { 
      return projectionType; 
    }
    
    /// Returns the projection with the given id, if existing. Returns NULL otherwise.
    static Projection* getProjection(int id) 
    {
      return projectionMap[id];
    }

  protected:
    /// Unique projection id.
    int projectionID;

    /// Type of this projection.
    ProjectionType projectionType;

    /// Static mapping from ids to projection objects. Used in \ref getProjection().
    static std::map<int, Projection*> projectionMap;
  };

}

#include "BallProject.h"
#include "CylinderProject.h"
#endif