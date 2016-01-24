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
      neighbour_inv(dim),
      oppVertex(dim),
      index(-1),
      elType(0)
  {
    neighbour.set((MacroElement*)(NULL));
    projection.set((Projection*)(NULL));
  }
  
  MacroElement::MacroElement(MacroElement const& that)
    : element(that.element),
      coord(that.coord),
      boundary(that.boundary),
      projection(that.projection),
      neighbour(that.neighbour),
      neighbour_inv(that.neighbour_inv),
      oppVertex(that.oppVertex),
      index(that.index),
      elType(that.elType)
  {
  }


  MacroElement::~MacroElement()
  {
    if (element)
      delete element;
  }


  MacroElement& MacroElement::operator=(MacroElement const& el)
  {
    if (this != &el) {
      coord = el.coord;
      boundary = el.boundary;
      projection = el.projection;
      oppVertex = el.oppVertex;
      index = el.index;
      elType = el.elType;
    }
    return *this;
  }
}
