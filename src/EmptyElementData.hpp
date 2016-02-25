/** \file EmptyElementData.h */

#pragma once

#include "Element.hpp"
#include "ElementData.hpp"
#include "FixVec.hpp"

namespace AMDiS
{
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


    EmptyElementData(ElementData* decorated = NULL)
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


    ElementData* clone() const
    {
      EmptyElementData* newObj = new EmptyElementData;
      newObj->decorated_ = ElementData::clone();
      return newObj;
    }


    const int getTypeID() const
    {
      return EMPTY_ED;
    }
  };

} // end namespace AMDiS
