#include <string>
#include <map>

#include "MacroElement.h"
#include "Boundary.h"
#include "FiniteElemSpace.h"
#include "Mesh.h"
#include "FixVec.h"
#include "FixVecConvert.h"

namespace AMDiS 
{
  MacroElement::MacroElement(int dim)
    : element(NULL),
      coord(dim),
      boundary(dim, INTERIOR),
      projection(dim),
      neighbour(dim),
      oppVertex(dim),
      index(-1), 
      elType(0),
      deserializedNeighbourIndices(NULL)
  {
    neighbour.set((MacroElement*)(NULL));
    projection.set((Projection*)(NULL));
  }


  MacroElement::~MacroElement()
  {
    if (element)
      delete element;    
  }


  MacroElement& MacroElement::operator=(const MacroElement &el)
  {
    if (this == &el)
      return *this;

    coord = el.coord;
    boundary = el.boundary;
    projection = el.projection;
    oppVertex = el.oppVertex;
    index = el.index;
    elType = el.elType;  
    
    return *this;
  }

}
