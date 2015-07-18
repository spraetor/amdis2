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


#include "ElementData.h"

namespace AMDiS {

  void ElementData::coarsenElementData(Element* parent, 
				       Element* thisChild,
				       Element* otherChild,
				       int elTypeParent) 
  {
    if (decorated) {
      decorated->coarsenElementData(parent, thisChild, otherChild, elTypeParent);
      delete decorated;
      decorated = NULL;
    }
  }

  bool ElementData::deleteDecorated(int typeID)
  {
    if (decorated) {
      if (decorated->isOfType(typeID)) {
	ElementData *tmp = decorated;
	decorated = decorated->decorated;
	delete tmp;
	tmp = NULL;
	return true;
      } else {
	return decorated->deleteDecorated(typeID);
      }
    } 
    return false;    
  }

  void ElementData::deleteDecorated()
  {
    if (decorated) {
      decorated->deleteDecorated();
      delete decorated;
    }
  }
  
  ElementData::~ElementData()
  {
  }

  void ElementData::serialize(std::ostream& out) 
  {
    std::string decoratedType;
    if (decorated) {
      decoratedType = decorated->getTypeName();
      out << decoratedType << "\n";
      decorated->serialize(out);
    } else {
      out << "NULL\n";
    }
  }

  void ElementData::deserialize(std::istream& in) 
  {
    TEST_EXIT(decorated == NULL)
      ("there are already decorated element data\n");
    std::string decoratedType;
    in >> decoratedType; 
    in.get();
    if (decoratedType != "NULL") {
      decorated = 
	CreatorMap<ElementData>::getCreator(decoratedType, "deserialization from file")->create();
      decorated->deserialize(in);
    } else {
      decorated = NULL;
    }
  }

}
