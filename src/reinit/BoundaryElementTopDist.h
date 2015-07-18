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




#ifndef BOUNDARYELEMENTTOPDIST_H
#define BOUNDARYELEMENTTOPDIST_H

#include "ElInfo.h"
#include "FixVec.h"
#include "ElementLevelSet.h"
#include "BoundaryElementDist.h"
#include "VelocityExt.h"

namespace reinit {

using namespace AMDiS;

class BoundaryElementTopDist : public BoundaryElementDist
{
public:
  BoundaryElementTopDist(ElementLevelSet *elLS_, 
			 int dim_,
			 VelocityExt *velExt_ = NULL)
    : BoundaryElementDist(elLS_, dim_),
      velExt(velExt_)
  {}

  ~BoundaryElementTopDist()
  {}

  /**
   * Calculates distance from the interface for all vertices of a boundary 
   * element.
   * Distance is here the topological distance.
   *
   * Return value: Status of element elInfo.
   */
  int calcDistOnBoundaryElement(ElInfo *elInfo,
				FixVec<double, VERTEX> &dVec);

 protected:
  /**
   * Calculates distance of vertex v to intersection edge through
   * sP1 and sP2 within element.
   * lambda implicitly gives the next point on the intersection edge.
   */
  double calcDistOnIntersecEdge(const WorldVector<double> &sP1,        
				const WorldVector<double> &sP2,        
				const WorldVector<double> &v,          
				double &lambda);

  //calculates the distance between two points
  double calc_dist_between_two_points(
       const WorldVector<double> &point1,                   
       const WorldVector<double> &point2);                   
  
  void cases_for_lambda(const WorldVector<double> &SP,           
			const WorldVector<double> &SP_projected, 
			const WorldVector<double> &point1,       
			const WorldVector<double> &point2,       
			double &lambda,                     
			double &dist);                           
  
  //calculates the projection of the point "v" onto a straight line 
  //given by two points "sP1" nad "sP2"
  void projected_on_a_straight_line(const WorldVector<double> &sP1,    
				    const WorldVector<double> &sP2,    
				    const WorldVector<double> &v,      
				    WorldVector<double> &v_p,          
				    double &lambda_out );              

 protected:
  /**
   * Object needed to extrapolate velocity from the interface.
   */  
  VelocityExt *velExt;
 };
 
}

using reinit::BoundaryElementTopDist;

#endif  // BOUNDARYELEMENTTOPDIST_H
