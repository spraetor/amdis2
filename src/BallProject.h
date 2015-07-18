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



/** \file BallProject.h */

#ifndef AMDIS_BALLPROJECT_H
#define AMDIS_BALLPROJECT_H

#include "MatrixVectorOperations.h"

namespace AMDiS {

  /** \brief
   * Projects world coordinates to the surface of a ball with given center and 
   * radius. Can be used as boundary or volume projection.
   */
  class BallProject : public Projection
  {
  public:
    /// Constructor.
    BallProject(int id, 
		ProjectionType type,
		WorldVector<double> &center,
		double radius) 
      : Projection(id, type),
	center_(center),
	radius_(radius)
    {}

    /// Destructor.
    virtual ~BallProject() {}

    /// Implementation of Projection::project();
    void project(WorldVector<double> &x) 
    {
      x -= center_;
      double norm = sqrt(x*x);
      TEST_EXIT(norm != 0.0)("can't project vector x\n");
      x *= radius_ / norm;
      x += center_;
    }

  protected:
    /// Center of the ball.
    WorldVector<double> center_;

    /// Radius of the ball.
    double radius_;
  };

}

#endif
