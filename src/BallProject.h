/** \file BallProject.h */

#pragma once

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

} // end namespace AMDiS
