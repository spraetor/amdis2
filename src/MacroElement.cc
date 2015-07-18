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


#include <string>
#include <map>
#include "MacroElement.h"
#include "Boundary.h"
#include "FiniteElemSpace.h"
#include "Mesh.h"
#include "FixVec.h"
#include "FixVecConvert.h"
#include "Serializer.h"

namespace AMDiS {

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


  void MacroElement::serialize(std::ostream &out)
  {
    // write element-tree
    out << element->getTypeName() << "\n";
    element->serialize(out);

    // write coords
    int size = coord.getSize();
    SerUtil::serialize(out, size);

    for (int i = 0; i < size; i++)
      coord[i].serialize(out);

    // write boundary
    boundary.serialize(out);

    // write projection
    size = projection.getSize();
    SerUtil::serialize(out, size);
    for (int i = 0; i < size; i++) {
      int id = projection[i] ? projection[i]->getID() : -1;
      SerUtil::serialize(out, id);
    }

    // write neighbour
    size = neighbour.getSize();
    SerUtil::serialize(out, size);
    for (int i = 0; i < size; i++) {
      int index = neighbour[i] ? neighbour[i]->getIndex() : -1;
      SerUtil::serialize(out, index);
    }
  
    // write oppVertex
    oppVertex.serialize(out);

    // write index
    SerUtil::serialize(out, index);

    // write elType
    SerUtil::serialize(out, elType);
  }


  void MacroElement::deserialize(std::istream &in)
  {
    FUNCNAME("MacroElement::deserialize()");

    // === Read element-tree. ===

    std::string typeName;
    in >> typeName;
    in.get();

    if (element) {
      TEST_EXIT(typeName == element->getTypeName())("wrong element type name\n");
    } else {
      if (typeName == "Line") 
	element = new Line(NULL);
      if (typeName == "Triangle")
	element = new Triangle(NULL);
      if (typeName == "Tetrahedron") 
	element = new Tetrahedron(NULL);
    }

    element->deserialize(in);

    // === Read coords. ===

    int size;
    SerUtil::deserialize(in, size);
    
    if (coord.getSize())
      TEST_EXIT(coord.getSize() == size)("invalid size\n");
    else
      coord.initSize(size);
    
    for (int i = 0; i < size; i++)
      coord[i].deserialize(in);

    boundary.deserialize(in);

    // === Read projections. ===

    SerUtil::deserialize(in, size);

    if (projection.getSize())
      TEST_EXIT(projection.getSize() == size)("invalid size\n");
    else
      projection.initSize(size);
    
    for (int i = 0; i < size; i++) {
      int id;
      SerUtil::deserialize(in, id);
      projection[i] = (id != -1) ? Projection::getProjection(id) : NULL;
    }

    // === Read neighbour indices. ===

    SerUtil::deserialize(in, size);
  
    TEST_EXIT(deserializedNeighbourIndices)
      ("Neighbour indices for deserializing not set!\n");
    
    deserializedNeighbourIndices->resize(size);

    if (neighbour.getSize())
      TEST_EXIT(neighbour.getSize() == size)("invalid size\n");
    else
      neighbour.initSize(size);
    
    for (int i = 0; i < size; i++)
      SerUtil::deserialize(in, (*deserializedNeighbourIndices)[i]);

    deserializedNeighbourIndices = NULL;

    // read oppVertex
    oppVertex.deserialize(in);

    SerUtil::deserialize(in, index);
    SerUtil::deserialize(in, elType);
  }

}
