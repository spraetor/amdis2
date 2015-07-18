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
      coord(dim, NO_INIT),
      boundary(dim, DEFAULT_VALUE, INTERIOR),
      projection(dim, NO_INIT),
      neighbour(dim, NO_INIT),
      oppVertex(dim, NO_INIT),
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
