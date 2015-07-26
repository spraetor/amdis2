#include "BoundaryObject.h"
#include "Mesh.h"
#include "FiniteElemSpace.h"
#include "BasisFunction.h"
// #include "MacroElement.h"
#include "Element.h"

namespace AMDiS {

  BoundaryObject::BoundaryObject()
    : elType(0),
      reverseMode(false),
      excludedSubstructures(0)
  { }


  BoundaryObject::BoundaryObject(Element *e, 
                        				 int eType, 
                        				 GeoIndex sObj, 
                        				 int iObj, 
                        				 bool rMode)
    : el(e),
      elIndex(e->getIndex()),
      elType(eType),
      subObj(sObj),
      ithObj(iObj),
      reverseMode(rMode),
      excludedSubstructures(0)
  { }


  bool BoundaryObject::computeReverseMode(BoundaryObject &obj0, 
                              					  BoundaryObject &obj1,
                              					  const FiniteElemSpace *feSpace,
                              					  BoundaryType boundary)
  {
    FUNCNAME("BoundaryObject::computeReverseMode()");

    bool reverseMode = false;
    
    TEST_EXIT_DBG(obj0.el && obj1.el)("BoundaryObject without element pointer.\n");

    switch (feSpace->getMesh()->getDim()) {
    case 2:
      reverseMode = true;
      break;

    case 3:
      TEST_EXIT_DBG(obj1.elType == 0)
	("Only 3D macro elements with level 0 are supported. This element has level %d!\n", 
	 obj1.elType);


      if (obj0.subObj == EDGE) {	
      	int el0_v0 = obj0.el->getVertexOfEdge(obj0.ithObj, 0);
      	int el1_v0 = obj0.el->getVertexOfEdge(obj1.ithObj, 0);
#if DEBUG != 0
      	int el0_v1 = obj0.el->getVertexOfEdge(obj0.ithObj, 1);
      	int el1_v1 = obj0.el->getVertexOfEdge(obj1.ithObj, 1);
#endif

      	const BasisFunction *basFcts = feSpace->getBasisFcts();
      	int nBasFcts = basFcts->getNumber();
      	std::vector<DegreeOfFreedom> localDofs0(nBasFcts), localDofs1(nBasFcts);
      	basFcts->getLocalIndices(obj0.el, feSpace->getAdmin(), localDofs0);
      	basFcts->getLocalIndices(obj1.el, feSpace->getAdmin(), localDofs1);
      
      	Mesh *mesh = feSpace->getMesh();
      
      	if (mesh->isPeriodicAssociation(boundary) == false) {
      	  TEST_EXIT_DBG(localDofs0[el0_v0] == localDofs1[el1_v0] ||
      			localDofs0[el0_v0] == localDofs1[el1_v1])
      	    ("This should not happen!\n");
      	  TEST_EXIT_DBG(localDofs0[el0_v1] == localDofs1[el1_v0] ||
      			localDofs0[el0_v1] == localDofs1[el1_v1])
      	    ("This should not happen!\n");
      
      	  if (localDofs0[el0_v0] != localDofs1[el1_v0])
      	    reverseMode = true; 	
      	} else {
      	  if (mesh->associated(localDofs0[el0_v0], localDofs1[el1_v0]) == false)
      	    reverseMode = true;
      	}
      }

      if (obj0.subObj == FACE && obj0.ithObj != 1) {
      	const BasisFunction *basFcts = feSpace->getBasisFcts();
      	int nBasFcts = basFcts->getNumber();
      	std::vector<DegreeOfFreedom> localDofs0(nBasFcts), localDofs1(nBasFcts);
      	basFcts->getLocalIndices(obj0.el, feSpace->getAdmin(), localDofs0);
      	basFcts->getLocalIndices(obj1.el, feSpace->getAdmin(), localDofs1);
      	
      	if (obj0.ithObj == 2 || obj0.ithObj == 3)
      	  reverseMode = (localDofs0[0] != localDofs1[0]);
      	  
      	if (obj0.ithObj == 0)
      	  reverseMode = (localDofs0[1] != localDofs1[1]);
      }
      break;

    default:
      ERROR_EXIT("This should not happen!\n");
    }

    return reverseMode;
  }


  bool BoundaryObject::operator==(const BoundaryObject& other) const
  {
    return (other.elIndex == elIndex && 
	    other.subObj == subObj && 
	    other.ithObj == ithObj);
  }


  bool BoundaryObject::operator!=(const BoundaryObject& other) const
  {
    return (other.elIndex != elIndex || 
	    other.subObj != subObj || 
	    other.ithObj != ithObj);
  }


  bool BoundaryObject::operator<(const BoundaryObject& other) const
  {
    if (elIndex == other.elIndex) {
      if (subObj == other.subObj)
        return ithObj < other.ithObj;

      return subObj < other.subObj;
    }
     
    return (elIndex < other.elIndex);
  }


  bool AtomicBoundary::operator==(const AtomicBoundary& other) const
  {
    return (rankObj == other.rankObj && 
	    neighObj == other.neighObj && 
	    type == other.type);
  }

  bool AtomicBoundary::operator!=(const AtomicBoundary& other) const
  {
    return (rankObj != other.rankObj ||
	    neighObj != other.neighObj ||
	    type != other.type);
  }

} // end namespace AMDiS
