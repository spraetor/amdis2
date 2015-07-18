#include "LeafData.h"
#include "Element.h"
#include "Mesh.h"

namespace AMDiS 
{
  bool LeafDataEstimatable::refineElementData(Element* parent, 
					      Element* child1,
					      Element* child2,
					      int elType)
 
  {
    ElementData::refineElementData(parent, child1, child2, elType);

    LeafDataEstimatable *child1Data = 
      new LeafDataEstimatable(child1->getElementData());
    LeafDataEstimatable *child2Data = 
      new LeafDataEstimatable(child2->getElementData());

    child1Data->setErrorEstimate(0, errorEstimate / 2.0);
    child2Data->setErrorEstimate(0, errorEstimate / 2.0);

    child1->setElementData(child1Data);
    child2->setElementData(child2Data);

    return true;
  }

  
  void LeafDataEstimatable::coarsenElementData(Element* parent, 
					       Element* thisChild,
					       Element* otherChild,
					       int elTypeParent)
 
  {    
#ifndef NDEBUG
    bool test = 
#endif
       otherChild->deleteElementData(ESTIMATABLE);
    TEST_EXIT_DBG(test)("couldn't delete LeafDataEstimatable at otherChild");

    parent->setElementData(new LeafDataEstimatable(parent->getElementData()));
    ElementData::coarsenElementData(parent, thisChild, otherChild, elTypeParent);
  }
  
  
  bool LeafDataEstimatableVec::refineElementData(Element* parent, 
						 Element* child1,
						 Element* child2,
						 int elType)
 
  {
    ElementData::refineElementData(parent, child1, child2, elType);
    child1->setElementData(new LeafDataEstimatableVec(child1->getElementData()));
    child2->setElementData(new LeafDataEstimatableVec(child2->getElementData()));

    return true;
  }
  

  void LeafDataEstimatableVec::coarsenElementData(Element* parent, 
						  Element* thisChild,
						  Element* otherChild,
						  int elTypeParent)
 
  {
    DBG_VAR(bool test =) otherChild->deleteElementData(ESTIMATABLE);
    TEST_EXIT_DBG(test)("couldn't delete LeafDataEstimatableVec at otherChild");
    parent->setElementData(new LeafDataEstimatableVec(parent->getElementData()));
    ElementData::coarsenElementData(parent, thisChild, otherChild, elTypeParent);
  }
  
  
  bool LeafDataCoarsenable::refineElementData(Element* parent, 
					      Element* child1,
					      Element* child2,
					      int elType)
  {
    ElementData::refineElementData(parent, child1, child2, elType);
    child1->setElementData(new LeafDataCoarsenable(child1->getElementData()));
    child2->setElementData(new LeafDataCoarsenable(child2->getElementData()));
    return true;
  }

  
  void LeafDataCoarsenable::coarsenElementData(Element* parent, 
					       Element* thisChild,
					       Element* otherChild,
					       int elTypeParent)
  {
    DBG_VAR(bool test =) otherChild->deleteElementData(COARSENABLE);
    TEST_EXIT_DBG(test)("couldn't delete LeafDataCoarsenable at otherChild");

    parent->setElementData(new LeafDataCoarsenable(parent->getElementData()));
    ElementData::coarsenElementData(parent, thisChild, otherChild, elTypeParent);
  }
  

  bool LeafDataCoarsenableVec::refineElementData(Element* parent, 
						 Element* child1,
						 Element* child2,
						 int elType)
  {
    ElementData::refineElementData(parent, child1, child2, elType);
    child1->setElementData(new LeafDataCoarsenableVec(child1->getElementData()));
    child2->setElementData(new LeafDataCoarsenableVec(child2->getElementData()));
    return true;
  }
  

  void LeafDataCoarsenableVec::coarsenElementData(Element* parent, 
						  Element* thisChild,
						  Element* otherChild,
						  int elTypeParent) 
  {
    DBG_VAR(bool test =) otherChild->deleteElementData(COARSENABLE);
    TEST_EXIT_DBG(test)("couldn't delete LeafDataCoarsenableVec at otherChild");
    parent->setElementData(new LeafDataCoarsenableVec(parent->getElementData()));
    ElementData::coarsenElementData(parent, thisChild, otherChild, elTypeParent);
  }

  
  bool LeafDataPeriodic::refineElementData(Element* parent, 
					   Element* child1,
					   Element* child2,
					   int elType) 
  {
    ElementData::refineElementData(parent, child1, child2, elType);

    Element* child[2] = {child1, child2};
    int dim = parent->getMesh()->getDim();
    LeafDataPeriodic *ld[2] = {NULL, NULL};
    std::list<LeafDataPeriodic::PeriodicInfo>::iterator it;

    for (it = periodicInfoList.begin(); it != periodicInfoList.end(); ++it) {
      BoundaryType type = it->type;
      int parentSide = it->elementSide;
      int mode = it->periodicMode;

      // for both children
      for (int i = 0; i < 2; i++) {
	// get childs side
	int sideOfChild = parent->getSideOfChild(i, parentSide, elType);
	if (sideOfChild == -1) 
	  continue;
	
	// create new leaf data if necessary
	if (!ld[i]) {
	  ld[i] = new LeafDataPeriodic(child[i]->getElementData());
	  child[i]->setElementData(ld[i]);
	}
    
	  // create new periodic coords
	DimVec<WorldVector<double> > coords(dim - 1, NO_INIT);
	
	// for each vertex of childs side
	for (int j = 0; j < dim; j++) {
	  // get parents vertex nr.
	  int childVertex = child[i]->getVertexOfPosition(INDEX_OF_DIM(dim - 1, dim),
							  sideOfChild, j);
	  int parentVertex = child[i]->getVertexOfParent(i, childVertex, elType);
	  
	  if (parentVertex == -1) {
	    // create coords for new vertex
	    WorldVector<double> newCoords;
	    newCoords = (*(it->periodicCoords))[0];
	    newCoords += (*(it->periodicCoords))[1];
	    newCoords *= 0.5;
	    coords[j] = newCoords;
	  } else {
	    int posAtSide = parent->getPositionOfVertex(parentSide, parentVertex);
	    coords[j] = (*(it->periodicCoords))[posAtSide];
	  }
	}
	
	ld[i]->addPeriodicInfo(mode, type, sideOfChild, &coords);
      }
    }

    return false;
  }
  
  
  LeafDataPeriodic::PeriodicInfo::PeriodicInfo(const PeriodicInfo &rhs) 
  {
    periodicMode = rhs.periodicMode;
    type = rhs.type;
    elementSide = rhs.elementSide;
    if (rhs.periodicCoords) {
      int dim = rhs.periodicCoords->getSize() - 1;
      periodicCoords = new DimVec<WorldVector<double> >(dim, NO_INIT);
      for (int i = 0; i < dim + 1; i++)
	(*periodicCoords)[i] = (*(rhs.periodicCoords))[i];
    } else {
      periodicCoords = NULL;
    }
  }
  
  
  LeafDataPeriodic::PeriodicInfo::PeriodicInfo(int mode, BoundaryType t, int side,
					       const DimVec<WorldVector<double> > *coords)
    : periodicMode(mode),
      type(t),
      elementSide(side),
      periodicCoords(NULL)
  {
    if (coords) {
      int dim = coords->getSize() - 1;
      periodicCoords = new DimVec<WorldVector<double> >(dim, NO_INIT);
      for (int i = 0; i < dim + 1; i++)
	(*periodicCoords)[i] = (*coords)[i];
    }
  }

} // end namespace AMDiS
