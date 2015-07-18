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


#include "Line.h"
#include "DOFAdmin.h"
#include "Mesh.h"
#include "CoarseningManager.h"
#include "FixVec.h"

namespace AMDiS {

  const int Line::vertexOfEdge[1][2] = {{0, 1}};

  const int Line::sideOfChild[2][2] = {{ 0,-1},
				       {-1, 1}};

  const int Line::vertexOfParent[2][2] = {{ 0,-1},
					  {-1, 1}};

  int Line::getVertexOfPosition(GeoIndex position, 
				int positionIndex,
				int vertexIndex) const
  {
    FUNCNAME("Line::getVertexOfPosition()");

    switch(position) {
    case VERTEX:
      return positionIndex;
      break;
    case CENTER:
      return vertexIndex;
      break;
    default:
      ERROR_EXIT("invalid position\n");
      return 0;
    }  
  }

  void Line::sortFaceIndices(int face, FixVec<int,WORLD> &vec) const
  {
    vec[0] = face;
  }

}
