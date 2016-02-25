#include "Line.hpp"

#include "CoarseningManager.hpp"
#include "DOFAdmin.hpp"
#include "FixVec.hpp"
#include "Mesh.hpp"

namespace AMDiS
{
  constexpr int Line::vertexOfEdge[1][2];
  constexpr int Line::sideOfChild[2][2];
  constexpr int Line::vertexOfParent[2][2];

  int Line::getVertexOfPosition(GeoIndex position,
                                int positionIndex,
                                int vertexIndex) const
  {
    FUNCNAME("Line::getVertexOfPosition()");

    switch(position)
    {
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

  void Line::sortFaceIndices(int face, FixVec<int,WORLD>& vec) const
  {
    vec[0] = face;
  }

} // end namespace AMDiS
