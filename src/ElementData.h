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



/** \file ElementData.h */

#ifndef AMDIS_ELEMENTDATA_H
#define AMDIS_ELEMENTDATA_H

#include "Serializable.h"
#include "CreatorMap.h"

namespace AMDiS {

  const int ESTIMATABLE = 1;
  const int COARSENABLE = 2;
  const int PERIODIC = 3;
  const int ELEMENT_REGION = 4;
  const int SURFACE_REGION = 5;

  /** \brief
   * Base class for element data. To allow to assign arbitrary data to arbitrary
   * elements at run time, the decorator pattern in combination with the
   * chain-of-responsibility pattern is applied. So only data have to be managed 
   * at each element which are used for this element at this time.
   */
  class ElementData : public Serializable
  {
  public:
    /// constructor
    ElementData(ElementData *dec = NULL) 
      : decorated(dec)
    {}

    /// destructor
    virtual ~ElementData();

    /// Refinement of parent to child1 and child2.
    virtual bool refineElementData(Element* parent, 
				   Element* child1,
				   Element* child2,
				   int elType)
    {
      if (decorated) {
	bool remove = 
	  decorated->refineElementData(parent, child1, child2, elType);

	if (remove) {
	  ElementData *tmp = decorated->decorated;
	  delete decorated;
	  decorated = tmp;
	}
      }
      return false;
    }

    ///
    virtual void coarsenElementData(Element* parent, 
				    Element* thisChild,
				    Element* otherChild,
				    int elTypeParent);

    /// Returns a copy of this ElementData object including all decorated data.
    virtual ElementData *clone() const 
    {
      if (decorated)
	return decorated->clone();

      return NULL;
    }

    /// Returns the id of element data type.
    virtual int getTypeID() const 
    {
      return 0;
    }

    /** \brief
     * Returns whether the ElemnetData object is of the type specified by 
     * typeName. Must return true, even if typeName describes a base class.
     */
    virtual bool isOfType(int typeID) const = 0;

    /// Implements Serializable::serialize().
    virtual void serialize(std::ostream& out);

    /// Implements Serializable::deserialize().
    virtual void deserialize(std::istream& in);

    /// Returns first element data in the chain which is of the spcified type.
    inline ElementData *getElementData(int typeID) 
    {
      if (this->isOfType(typeID)) {
	return this;
      } else {
	if (decorated)
	  return decorated->getElementData(typeID);
      }
      return NULL;
    }

    inline ElementData *getDecorated(int typeID) 
    { 
      if (decorated)
	return decorated->getElementData(typeID);
      
      return NULL;
    }

    /** \ref
     * Search the \ref decorated chain for a specific type ID, and delets
     * this entry.
     */
    bool deleteDecorated(int typeID);

    /// Delets the whole \ref decorated chain.
    void deleteDecorated();

    inline ElementData *getDecorated() 
    { 
      return decorated; 
    }

    inline void setDecorated(ElementData *d) 
    {
      if (getTypeID() == 1)
	if (d != NULL)
	  std::cout << "leafdata decorated with nonzero" << std::endl;

      decorated = d;
    }

  protected:
    /// Pointer to next ElementData object in the chain of responsibility.
    ElementData *decorated;
  };

}

#endif
