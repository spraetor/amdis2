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


#include "ElementUpdate_2d.h"

namespace reinit {

using namespace AMDiS;

double ElementUpdate_2d::calcElementUpdate(const FixVec<WorldVector<double> *, VERTEX> &vert,
					   FixVec<double, VERTEX> &uhVal)
{
  WorldVector<double> &v0 = *(vert[0]);
  WorldVector<double> &v1 = *(vert[1]);
  WorldVector<double> &v2 = *(vert[2]);

  int dim = Global::getGeo(WORLD);
  for (int i = 0; i < dim; i++) {
    xhminusYh[i] = v2[i] - v0[i];
    zhminusYh[i] = v1[i] - v0[i];
    xhminusZh[i] = v2[i] - v1[i];
  }
  
  double norm_zhminusYh = sqrt(zhminusYh * zhminusYh);
  double norm_xhminusYh = sqrt(xhminusYh * xhminusYh);
  double norm_xhminusZh = sqrt(xhminusZh * xhminusZh);
  
  double delta = (uhVal[1] - uhVal[0]) / norm_zhminusYh; 
  double sP = xhminusYh * zhminusYh;
  double c_alpha = sP / (norm_xhminusYh * norm_zhminusYh);
  sP = xhminusZh * zhminusYh;
  double c_beta = -sP / (norm_xhminusZh * norm_zhminusYh);
  
  double update = 0;
  if (c_alpha <= delta) {
    update = uhVal[0] + norm_xhminusYh;
    //save barycentric coordinates for calculation of the velocity
    if (velExt != NULL)      
      velExt->setBarycentricCoords_2D(1,0,0);      
  } else if (delta <= -c_beta) {
    update = uhVal[1] + norm_xhminusZh;
    //save barycentric coordinates for calculation of the velocity
    if (velExt != NULL)      
      velExt->setBarycentricCoords_2D(0,1,0);      
  } else {
    update = uhVal[0] + 
      (c_alpha * delta + sqrt((1 - c_alpha * c_alpha) * (1 - delta * delta))) *
      norm_xhminusYh;
    //calculate and save barycentric coordinates for calculation of the velocity
    if (velExt != NULL)      
      velExt->calcBarycentricCoords_2D(delta, c_alpha, norm_zhminusYh, norm_xhminusYh);      
  }
  
  return update;
}

}
