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


#include "SurfaceRegion_ED.h"

namespace AMDiS {
  
  bool SurfaceRegion_ED::refineElementData(Element* parent, 
					   Element* child1,
					   Element* child2,
					   int elType)
  {
    ElementData::refineElementData(parent, child1, child2, elType);
    
    int sideOfChild;
    SurfaceRegion_ED *surfaceRegion;
    
    sideOfChild = parent->getSideOfChild(0, side, elType);
    if (sideOfChild >= 0) {
      surfaceRegion = new SurfaceRegion_ED(child1->getElementData());
      surfaceRegion->setSide(sideOfChild);
      surfaceRegion->setRegion(region);
      child1->setElementData(surfaceRegion);
    }    
    
    sideOfChild = parent->getSideOfChild(1, side, elType);
    if (sideOfChild >= 0) {
      surfaceRegion = new SurfaceRegion_ED(child2->getElementData());
      surfaceRegion->side = sideOfChild;
      surfaceRegion->region = region;
      child2->setElementData(surfaceRegion);
    }
    
    return false;
  }
  
  
  ElementData *SurfaceRegion_ED::clone() const 
  { 
    SurfaceRegion_ED *newObj = new SurfaceRegion_ED;
    newObj->side = side;
    newObj->region = region;
    newObj->decorated = ElementData::clone();
    return newObj; 
  }

} // end namespace AMDiS
