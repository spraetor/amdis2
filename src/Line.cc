#include "Line.h"
#include "DOFAdmin.h"
#include "Mesh.h"
#include "CoarseningManager.h"
#include "FixVec.h"

namespace AMDiS 
{
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

} // end namespace AMDiS
