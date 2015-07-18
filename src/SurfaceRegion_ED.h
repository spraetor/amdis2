/** \file SurfaceRegion_ED.h */

#pragma once

#include "ElementData.h"
#include "FixVec.h"

namespace AMDiS 
{

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

} // end namespace AMDiS
