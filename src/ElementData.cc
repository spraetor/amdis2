#include "ElementData.h"

namespace AMDiS 
{
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

} // end namespace AMDiS
