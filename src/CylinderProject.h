/** \file CylinderProject.h */

#pragma once

#include <cmath>
#include "Log.h"
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
      double norm = std::sqrt(direction_*direction_);
      direction_ *= 1.0 / norm;
    }

    /// Implementation of Projection::project();
    void project(WorldVector<double> &x) 
    {
      x -= center_;
      WorldVector<double> v1 = direction_; v1 *= (x*direction_);
      WorldVector<double> v2 = x; v2 -= v1;
      double norm = std::sqrt(v2 * v2);
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

} // end namespace AMDiS
