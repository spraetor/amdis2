#pragma once

#include "ElementData.hpp"
#include "FixVec.hpp"

namespace AMDiS
{

  class SurfaceRegion_ED : public ElementData
  {
  public:
    virtual bool isOfType(int typeID) const override
    {
      if (typeID == SURFACE_REGION)
        return true;
      return false;
    }

    struct Creator : public CreatorInterface<ElementData>
    {
      virtual ElementData* create() override
      {
        return new SurfaceRegion_ED;
      }
    };

    SurfaceRegion_ED(ElementData* decorated = NULL)
      : ElementData(decorated),
        side(-1),
        region(-1)
    {}

    virtual bool refineElementData(Element* parent,
                                   Element* child1,
                                   Element* child2,
                                   int elType) override;

    virtual ElementData* clone() const override;

    virtual int getTypeID() const override
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
