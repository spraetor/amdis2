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


#include "BoundaryElementTopDist.h"

namespace reinit
{
  using namespace AMDiS;
  int
  BoundaryElementTopDist::calcDistOnBoundaryElement(
    ElInfo* elInfo,
    FixVec<double, VERTEX>& dVec)
  {
    //Cartesian coordinates of the 3 or 4 vertices of a simplex
    FixVec<WorldVector<double>,VERTEX> Vert(dim);
    //normal unit vector in Cartesian coordinates
    WorldVector<double> normalVec;
    //Cartesian coordinates of dim intersection points
    FixVec<WorldVector<double>,DIMEN> planeVecs_DIMEN(dim);
    //Cartesian coordinates of 3 or 4 intersection points (3d)
    FixVec<WorldVector<double>,VERTEX> planeVecs(dim);
    //barycentric coordinates of the intersection point
    //(straight line - normal)
    DimVec<double> SP_Vec(dim);
    //Cartesian coordinates of SP_Vec
    WorldVector<double> SP;

    double lambda;
    double h_dist = 0.0;
    double h_dist2 = 0.0;
    double dist = 0.0;

    // for VelocityExt
    int edgeUpdateInd = -1;
    int sP1NextEdge = -1;
    int sP2NextEdge = -1;
    double lambdaNextPt = 0.0;

    // Get intersection information.
    int  elStatus = elLS->createElementLevelSet(elInfo);
    if (elStatus != ElementLevelSet::LEVEL_SET_BOUNDARY)
      return elStatus;

    //bar. coordinates of the intersection points
    VectorOfFixVecs<DimVec<double>>* elIntersecPoints =
                                   elLS->getElIntersecPoints();

    //Cartesian coordinates of the vertices
    for( int i=0; i<=dim; i++)
    {
      Vert[i] = elInfo->getCoord(i);
    }

    //Cartesian coordinates of the intersection point
    //attention: in 3d it is possible that there are 4 intersection
    //points, but her we save only 3 of them
    for(int i=0; i<dim/*ElementLevelSet::getNumElIntersecPoints()*/; i++)
    {
      elInfo->coordToWorld((*elIntersecPoints)[i], (planeVecs_DIMEN[i]));
    }

    //calculating of the normal unit vector
    calcNormal(planeVecs_DIMEN, normalVec);

    //2D or 3D
    if( dim == 2)
    {
      //loop over all vertices
      for(int i=0; i<=dim; i++)
      {
        //Cartesian coordinates of the intersection point
        //(straight line- normal)
        h_dist = 0.0;
        for (int j=0; j<dim; ++j)
        {
          h_dist += (Vert[i][j]-planeVecs_DIMEN[0][j]) * normalVec[j];
        }
        SP[0]=Vert[i][0]-h_dist*normalVec[0];
        SP[1]=Vert[i][1]-h_dist*normalVec[1];

        //barycentric coordinates of the intersection point
        elInfo->worldToCoord(SP, SP_Vec);

        //calculating of the distance
        //intersection point inside of the 2d-simplex?
        if(SP_Vec[0]>=0 && SP_Vec[1]>=0 && SP_Vec[2]>=0)
        {
          dist = fabs(h_dist);
          dVec[i] = dist;
          //save barycentric coordinates of intersection point,
          //for calculation of the velocity
          if(velExt != NULL)
          {
            velExt->setBarycentricCoords_2D_boundary(SP_Vec[0], SP_Vec[1], SP_Vec[2], i);
          }
        }
        //dist = min of the distance to the two intersection points
        else
        {
          for (int j=0; j<elLS->getNumElIntersecPoints(); j++)
          {
            h_dist = 0.0;
            for(int k=0; k<dim; k++)
            {
              h_dist += (Vert[i][k]-planeVecs_DIMEN[j][k])*
                        (Vert[i][k]-planeVecs_DIMEN[j][k]);
            }
            if (j == 0)
            {
              h_dist2 = h_dist;

              // Store index of intersection point which is next
              // point on interface.
              edgeUpdateInd = j;
            }
            else if (h_dist < h_dist2)
            {
              h_dist2 = h_dist;

              // Store index of intersection point which is next
              // point on interface.
              edgeUpdateInd = j;
            }
          }
          dist = sqrt(h_dist2);
          dVec[i] = dist;
          //save barycentric coordinates of intersection point,
          //for calculation of the velocity
          if(velExt != NULL)
          {
            velExt->setBarycentricCoords_2D_boundary(
              (*elIntersecPoints)[edgeUpdateInd][0],
              (*elIntersecPoints)[edgeUpdateInd][1],
              (*elIntersecPoints)[edgeUpdateInd][2],
              i);
          }
        }//end of the  else-part according to "if(!normalDistCalc)"
      }//end of the loop over all vertices
    }//end of the 2d-case

    else         //dim ==3
    {
      //Cartesian coordinates of the intersection point
      //attention: in 3d it is possible that there are 4 intersection
      //points, here we save them all
      for(int i=0; i<elLS->getNumElIntersecPoints(); i++)
      {
        elInfo->coordToWorld((*elIntersecPoints)[i], (planeVecs[i]));
      }

      //loop over all 4 vertices
      for ( int i=0; i<=dim; i++)
      {
        h_dist = 0.0;

        //intersection point of the normal with the intersection plane
        for ( int j = 0; j<dim; j++)
        {
          h_dist+=(Vert[i][j]-planeVecs[0][j])*normalVec[j];
        }
        for ( int j=0; j<dim; j++)
        {
          SP[j]=Vert[i][j]-h_dist*normalVec[j];
        }

        //barycentric coordinates of the intersection point
        elInfo->worldToCoord(SP, SP_Vec);

        //intersection point inside of the 3d-simplex
        //distance calculated with the unit normal vector of the
        //intersection plane
        if (SP_Vec[0]>=0 && SP_Vec[1]>=0 && SP_Vec[2]>=0 && SP_Vec[3]>=0)
        {
          dist = fabs(h_dist);
          dVec[i] = dist;
          //save barycentric coordinates of intersection point,
          //for calculation of the velocity
          if(velExt != NULL)
          {
            velExt->setBarycentricCoords_3D_boundary(SP_Vec[0], SP_Vec[1], SP_Vec[2], SP_Vec[3], i);
          }
        }
        //dist = min of the distances to the cutting straight lines at
        //the sides of the simplex
        else if(elLS->getNumElIntersecPoints() == 3)
        {
          dist = 1.e15;
          for (int j=0; j<dim; j++)
          {
            h_dist = calcDistOnIntersecEdge(planeVecs[j],
                                            planeVecs[(j+1)%3],
                                            Vert[i],
                                            lambda);
            // dist = (dist < h_dist) ? dist : h_dist;
            if(h_dist<dist)
            {
              dist = h_dist;
              if (velExt != NULL)
              {
                sP1NextEdge = j;
                sP2NextEdge = (j+1)%3;
                lambdaNextPt = lambda;
              }
            }
            //  if (velExt != NULL) {
            // 		    sP1NextEdge = j;
            // 		    sP2NextEdge = (j+1)%3;
            // 		    lambdaNextPt = lambda;

          }
          dVec[i] = dist;
          //save barycentric coordinates of intersection point,
          //for calculation of the velocity
          if(velExt != NULL)
          {
            velExt->calcBarycentricCoords_3D_boundary(
              (*elIntersecPoints)[sP1NextEdge],
              (*elIntersecPoints)[sP2NextEdge],
              lambdaNextPt,
              i);
          }
        }
        else if(elLS->getNumElIntersecPoints() == 4)
        {
          dist = 1.e15;
          for (int j=0; j<=1; j++)
          {
            h_dist = calcDistOnIntersecEdge(planeVecs[j],
                                            planeVecs[(j+2)%3],
                                            Vert[i],
                                            lambda);
            //  dist = (dist < h_dist) ? dist : h_dist;

            // 		  if (velExt != NULL) {
            // 		    sP1NextEdge = j;
            // 		    sP2NextEdge = (j+2)%3;
            // 		    lambdaNextPt = lambda;
            // 		  }
            if(h_dist<dist)
            {
              dist = h_dist;
              if (velExt != NULL)
              {
                sP1NextEdge = j;
                sP2NextEdge = (j+2)%3;
                lambdaNextPt = lambda;
              }
            }
          }

          h_dist = calcDistOnIntersecEdge(planeVecs[1],
                                          planeVecs[2],
                                          Vert[i],
                                          lambda);
          //  dist = (dist < h_dist) ? dist : h_dist;

          // 	      if (velExt != NULL) {
          // 		sP1NextEdge = 1;
          // 		sP2NextEdge = 2;
          // 		lambdaNextPt = lambda;
          // 	      }
          if(h_dist<dist)
          {
            dist = h_dist;
            if (velExt != NULL)
            {
              sP1NextEdge = 1;
              sP2NextEdge = 2;
              lambdaNextPt = lambda;
            }
          }

          h_dist = calcDistOnIntersecEdge(planeVecs[0],
                                          planeVecs[3],
                                          Vert[i],
                                          lambda);
          //  dist = (dist < h_dist) ? dist : h_dist;

          // 	      if (velExt != NULL) {
          // 		sP1NextEdge = 0;
          // 		sP2NextEdge = 3;
          // 		lambdaNextPt = lambda;
          // 	      }
          if(h_dist<dist)
          {
            dist = h_dist;
            if (velExt != NULL)
            {
              sP1NextEdge = 0;
              sP2NextEdge = 3;
              lambdaNextPt = lambda;
            }
          }

          dVec[i] = dist;
          //save barycentric coordinates of intersection point,
          //for calculation of the velocity
          if(velExt != NULL)
          {
            velExt->calcBarycentricCoords_3D_boundary(
              (*elIntersecPoints)[sP1NextEdge],
              (*elIntersecPoints)[sP2NextEdge],
              lambdaNextPt,
              i);
          }
        }
      }//end of the loop over all 4 vertices
    }//end of the else-part according to "if( dim == 2)"

    return elStatus;
  }

  double
  BoundaryElementTopDist::calcDistOnIntersecEdge(
    const WorldVector<double>& sP1,
    const WorldVector<double>& sP2,
    const WorldVector<double>& v,
    double& lambda)
  {
    WorldVector<double> SP_projected;
    double dist = 0.0;

    // Calculate projection of v on straight line through sP1 and sP2.
    projected_on_a_straight_line(sP1, sP2, v, SP_projected, lambda);

    // Get distance and next point on interface (implicitly: lambda).
    cases_for_lambda(v, SP_projected, sP1, sP2, lambda, dist);

    return dist;
  }

  double
  BoundaryElementTopDist::calc_dist_between_two_points (
    const WorldVector<double>& point1,
    const WorldVector<double>& point2)
  {
    double h_dist = 0.0;

    for(int i=0; i<dim; i++)
    {
      h_dist += (point1[i] - point2[i])*(point1[i] - point2[i]);
    }
    h_dist = sqrt(h_dist);

    return h_dist;
  }

  void
  BoundaryElementTopDist::cases_for_lambda(
    const WorldVector<double>& SP,
    const WorldVector<double>& SP_projected,
    const WorldVector<double>& point1,
    const WorldVector<double>& point2,
    double& lambda,
    double& dist)
  {
    if(0 <= lambda && lambda <= 1)
    {
      dist = calc_dist_between_two_points(SP, SP_projected);
    }
    if(lambda < 0.0)
    {
      lambda = 0;
      dist = calc_dist_between_two_points(SP, point1);
    }
    if(1.0 < lambda)
    {
      lambda = 1;
      dist = calc_dist_between_two_points(SP, point2);
    }
  }

  void
  BoundaryElementTopDist::projected_on_a_straight_line(
    const WorldVector<double>& sP1,
    const WorldVector<double>& sP2,
    const WorldVector<double>& v,
    WorldVector<double>& v_p,
    double& lambda)
  // lambda < 0:  v_p---sP1---sP2
  // 0<lambda<1:  sP1---v_p---sP2
  // 1 < lambda:  sP1---sP2---v_p
  {
    lambda = 0.0;
    WorldVector<double> sP1_sP2;
    WorldVector<double> sP1_v;
    double norm_sP1_sP2 = 0.0;

    for (int i=0; i<dim; ++i)
    {
      sP1_sP2[i]=sP2[i]-sP1[i];
    }
    for (int i=0; i<dim; ++i)
    {
      sP1_v[i]=v[i]-sP1[i];
    }
    for (int i=0; i<dim; ++i)
    {
      norm_sP1_sP2 += sP1_sP2[i]*sP1_sP2[i];
    }
    // norm_sP1_sP2 = sqrt(norm_sP1_sP2);

    for (int i=0; i<dim; ++i)
    {
      lambda += sP1_v[i]*sP1_sP2[i];
    }
    lambda = lambda / norm_sP1_sP2;

    for (int i=0; i<dim; i++)
    {
      v_p[i] = sP1[i] + lambda * sP1_sP2[i];
    }
  }

}
