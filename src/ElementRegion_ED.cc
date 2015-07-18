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


#include "ElementRegion_ED.h"

namespace AMDiS {
  
  bool ElementRegion_ED::refineElementData(Element* parent, 
					   Element* child1,
					   Element* child2,
					   int elType)
  {
    ElementData::refineElementData(parent, child1, child2, elType);
    
    ElementRegion_ED *ep;
    
    ep = new ElementRegion_ED(child1->getElementData());
    ep->region = region;
    child1->setElementData(ep);
    
    ep = new ElementRegion_ED(child2->getElementData());
    ep->region = region;
    child2->setElementData(ep);
    
    return false;
  }

  
  ElementData *ElementRegion_ED::clone() const
  { 
    ElementRegion_ED *newObj = new ElementRegion_ED;
    newObj->region = region;
    newObj->decorated = ElementData::clone();
    return newObj; 
  }
  
} // end namespace AMDiS
