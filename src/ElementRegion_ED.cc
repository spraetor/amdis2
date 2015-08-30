#include "ElementRegion_ED.h"
#include "Element.h"

namespace AMDiS
{
  bool ElementRegion_ED::refineElementData(Element* parent,
      Element* child1,
      Element* child2,
      int elType)
  {
    ElementData::refineElementData(parent, child1, child2, elType);

    ElementRegion_ED* ep;

    ep = new ElementRegion_ED(child1->getElementData());
    ep->region = region;
    child1->setElementData(ep);

    ep = new ElementRegion_ED(child2->getElementData());
    ep->region = region;
    child2->setElementData(ep);

    return false;
  }


  ElementData* ElementRegion_ED::clone() const
  {
    ElementRegion_ED* newObj = new ElementRegion_ED;
    newObj->region = region;
    newObj->decorated = ElementData::clone();
    return newObj;
  }

} // end namespace AMDiS
