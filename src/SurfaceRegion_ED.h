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



/** \file SurfaceRegion_ED.h */

#ifndef AMDIS_SURFACEREGION_H
#define AMDIS_SURFACEREGION_H

#include "ElementData.h"
#include "FixVec.h"

namespace AMDiS {

  class SurfaceRegion_ED : public ElementData
  {
  public:
    bool isOfType(int typeID) const 
    {
      if (typeID == SURFACE_REGION) 
	return true;
      return false;
    }

    class Creator : public CreatorInterface<ElementData>
    {
    public:
      ElementData* create() 
      {
	return new SurfaceRegion_ED;
      }
    };

    SurfaceRegion_ED(ElementData *decorated = NULL)
      : ElementData(decorated),
	side(-1),
	region(-1)
    {}

    bool refineElementData(Element* parent, 
			   Element* child1,
			   Element* child2,
			   int elType);

    ElementData *clone() const ;

    std::string getTypeName() const 
    { 
      return "SurfaceRegion_ED"; 
    }

    int getTypeID() const 
    { 
      return SURFACE_REGION; 
    }

    void setSide(int s) 
    { 
      side = s; 
    }

    int getSide() const 
    { 
      return side; 
    }

    void setRegion(int r) 
    { 
      region = r; 
    }

    int getRegion() const 
    { 
      return region; 
    }

  protected:
    int side;
    int region;
  };

}

#endif
