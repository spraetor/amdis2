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



/** \file CylinderProject.h */

#ifndef AMDIS_CYLINDERPROJECT_H
#define AMDIS_CYLINDERPROJECT_H

#include "MatrixVectorOperations.h"

namespace AMDiS {

  /** \brief
   * Projects world coordinates to the surface of a cylinder with given center, 
   * radius and direction. Can be used as boundary or volume projection.
   */
  class CylinderProject : public Projection
  {
  public:
    /// Constructor.
    CylinderProject(int id,
		    ProjectionType type,
		    WorldVector<double> &c,
		    WorldVector<double> &d,
		    double r) 
      : Projection(id, type),
	center_(c),
	direction_(d),
	radius_(r)
    {
      double norm = sqrt(direction_*direction_);
      direction_ *= 1.0 / norm;
    }

    /// Destructor.
    virtual ~CylinderProject() {}

    /// Implementation of Projection::project();
    void project(WorldVector<double> &x) 
    {
      x -= center_;
      WorldVector<double> v1 = direction_; v1 *= (x*direction_);
      WorldVector<double> v2 = x; v2 -= v1;
      double norm = sqrt(v2 * v2);
      TEST_EXIT(norm != 0.0)("can't project vector x\n");
      v2 *= 1.0 / norm;
      x = v2; x *= radius_; x += v1;
      x += center_;
    }

  protected:
    /// Center of the cylinder.
    WorldVector<double> center_;

    /// Direction of the cylinder.
    WorldVector<double> direction_;

    /// Radius of the cylinder.
    double radius_;
  };

}

#endif
