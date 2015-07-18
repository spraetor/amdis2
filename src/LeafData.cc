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


#include "LeafData.h"
#include "Element.h"
#include "Mesh.h"
#include "Serializer.h"

namespace AMDiS {

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

  void LeafDataEstimatable::serialize(std::ostream& out)
  {
    ElementData::serialize(out);
    SerUtil::serialize(out, errorEstimate);
  }

  void LeafDataEstimatable::deserialize(std::istream& in)
  {
    ElementData::deserialize(in);
    SerUtil::deserialize(in, errorEstimate);
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

  void LeafDataEstimatableVec::serialize(std::ostream& out) 
  {
    ElementData::serialize(out);
    unsigned int size = errorEstimate.size();
    SerUtil::serialize(out, size);
    
    for (std::map<int, double>::iterator it = errorEstimate.begin(); 
	 it != errorEstimate.end(); ++it) {
      SerUtil::serialize(out, it->first);
      SerUtil::serialize(out, it->second);
    }
  }
  
  void LeafDataEstimatableVec::deserialize(std::istream& in) 
  {
    ElementData::deserialize(in);
    unsigned size;
    SerUtil::deserialize(in, size);
    for (unsigned int i = 0; i < size; i++) {
      int index;
      double estimate;
      SerUtil::deserialize(in, index);
      SerUtil::deserialize(in, estimate);
      errorEstimate[index] = estimate;
    }
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

  void LeafDataCoarsenable::serialize(std::ostream& out) 
  {
    ElementData::serialize(out);
    SerUtil::serialize(out, coarseningError);
  }

  void LeafDataCoarsenable::deserialize(std::istream& in) 
  {
    ElementData::deserialize(in);
    SerUtil::deserialize(in, coarseningError);
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

  void LeafDataCoarsenableVec::serialize(std::ostream& out) 
  {
    ElementData::serialize(out);
    unsigned int size = coarseningError.size();
    SerUtil::serialize(out, size);
    
    for (std::map<int, double>::iterator it = coarseningError.begin(); 
	 it != coarseningError.end(); ++it) {
      SerUtil::serialize(out, it->first);
      SerUtil::serialize(out, it->second);
    }
  }
  
  void LeafDataCoarsenableVec::deserialize(std::istream& in) 
  {
    ElementData::deserialize(in);
    unsigned int size;
    SerUtil::deserialize(in, size);
    for (unsigned int i = 0; i < size; i++) {
      int index;
      double estimate;
      SerUtil::deserialize(in, index);
      SerUtil::deserialize(in, estimate);
      coarseningError[index] = estimate;
    }
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

  void LeafDataPeriodic::serialize(std::ostream& out) 
  {
    ElementData::serialize(out);
    unsigned int size = periodicInfoList.size();
    SerUtil::serialize(out, size);      
    for (std::list<PeriodicInfo>::iterator it = periodicInfoList.begin(); 
	 it != periodicInfoList.end(); ++it)
      it->serialize(out);
  }
  
  void LeafDataPeriodic::deserialize(std::istream& in) 
  {
    ElementData::deserialize(in);
    unsigned int size;
    SerUtil::deserialize(in, size);
    periodicInfoList.resize(size);
    for (std::list<PeriodicInfo>::iterator it = periodicInfoList.begin(); 
	 it != periodicInfoList.end(); ++it)
      it->deserialize(in);
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

  void LeafDataPeriodic::PeriodicInfo::serialize(std::ostream &out) 
  {
    SerUtil::serialize(out, periodicMode);
    SerUtil::serialize(out, type);
    SerUtil::serialize(out, elementSide);
    if (periodicCoords) {
      int size = periodicCoords->getSize();
      SerUtil::serialize(out, size);
      for (int i = 0; i < size; i++)
	(*periodicCoords)[i].serialize(out);
    } else {
      int size = 0;
      SerUtil::serialize(out, size);
    }
  }
  
  void LeafDataPeriodic::PeriodicInfo::deserialize(std::istream &in)
  {
    SerUtil::deserialize(in, periodicMode);
    SerUtil::deserialize(in, type);
    SerUtil::deserialize(in, elementSide);
    
    int size;
    SerUtil::deserialize(in, size);
    if (periodicCoords) 
      delete periodicCoords;
    if (size == 0) {
      periodicCoords = NULL;
    } else {
      periodicCoords = new DimVec<WorldVector<double> >(size-1, NO_INIT);
      for (int i = 0; i < size; i++) 
	(*periodicCoords)[i].deserialize(in);
    }
  }

}
