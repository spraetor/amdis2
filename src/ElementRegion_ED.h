/** \file ElementRegionED.h */

#pragma once

#include "ElementData.h"
#include "FixVec.h"

namespace AMDiS 
{
  class ElementRegion_ED : public ElementData
  {
  public:
    bool isOfType(int typeID) const 
    {
      if (typeID == ELEMENT_REGION) 
	return true;
      return false;
    }

    class Creator : public CreatorInterface<ElementData>
    {
    public:
      ElementData* create() 
      {
	return new ElementRegion_ED;
      }
    };

    ElementRegion_ED(ElementData *decorated = NULL)
      : ElementData(decorated),
	region(-1)
    {}

    bool refineElementData(Element* parent, 
			   Element* child1,
			   Element* child2,
			   int elType);

    ElementData *clone() const;

    std::string getTypeName() const 
    { 
      return "ElementRegion_ED"; 
    }

    int getTypeID() const 
    { 
      return ELEMENT_REGION; 
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
    int region;
  };

} // end namespace AMDiS
