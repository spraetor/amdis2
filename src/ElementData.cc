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

  void ElementData::serialize(std::ostream& out) 
  {
    std::string decoratedType;
    if (decorated) {
      decoratedType = decorated->getTypeName();
      out << decoratedType << "\n";
      decorated->serialize(out);
    } else {
      out << "NULL\n";
    }
  }

  void ElementData::deserialize(std::istream& in) 
  {
    TEST_EXIT(decorated == NULL)
      ("there are already decorated element data\n");
    std::string decoratedType;
    in >> decoratedType; 
    in.get();
    if (decoratedType != "NULL") {
      decorated = 
	CreatorMap<ElementData>::getCreator(decoratedType, "deserialization from file")->create();
      decorated->deserialize(in);
    } else {
      decorated = NULL;
    }
  }

} // end namespace AMDiS
