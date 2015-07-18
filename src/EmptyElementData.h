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



/** \file EmptyElementData.h */

#ifndef AMDIS_EMPTYELEMENTDATA_H
#define AMDIS_EMPTYELEMENTDATA_H

#include "Element.h"
#include "ElementData.h"
#include "FixVec.h"

namespace AMDiS {

  const int EMPTY_ED = 6;

  class EmptyElementData : public ElementData
  {
  public:
    bool isOfType(int typeID) const 
    {
      if (typeID == EMPTY_ED) 
	return true;
      return false;
    }

    class Creator : public CreatorInterface<ElementData>
    {
    public:
      ElementData* create() 
      {
	return new EmptyElementData;
      }
    };
    

    EmptyElementData(ElementData *decorated = NULL)
      : ElementData(decorated)
    {}

    
    bool refineElementData(Element* parent, 
			   Element* child1,
			   Element* child2,
			   int elType)
    {
      ElementData::refineElementData(parent, child1, child2, elType);
      child1->setElementData(new EmptyElementData(child1->getElementData()));
      child2->setElementData(new EmptyElementData(child2->getElementData()));
      return false;
    }
    

    ElementData *clone() const 
    { 
      EmptyElementData *newObj = new EmptyElementData;
      newObj->decorated_ = ElementData::clone();
      return newObj; 
    }

    
    std::string getTypeName() const 
    { 
      return "EmptyElementData"; 
    }
    

    const int getTypeID() const 
    { 
      return EMPTY_ED; 
    }
  };

}

#endif
