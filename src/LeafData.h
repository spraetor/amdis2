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



/** \file LeafData.h */

#ifndef AMDIS_LEAFDATA_H
#define AMDIS_LEAFDATA_H

#include <list>
#include <map>
#include "Serializable.h"
#include "FixVec.h"
#include "ElementData.h"
#include "Boundary.h"

namespace AMDiS {

  class LeafDataEstimatableInterface
  {
  public:
    virtual ~LeafDataEstimatableInterface() {};
    virtual void setErrorEstimate(int, double) = 0;
    virtual double getErrorEstimate(int) = 0;
  };

  class LeafDataEstimatable : public ElementData,
			      public LeafDataEstimatableInterface
  {
  public:
    inline bool isOfType(int type) const 
    {
      if (type == ESTIMATABLE)
	return true;
      
      return false;
    }

    class Creator : public CreatorInterface<ElementData>
    {
    public:
      ElementData* create() 
      {
	return new LeafDataEstimatable;
      }
    };

    /// constructor
    LeafDataEstimatable(ElementData *decorated = NULL)
      : ElementData(decorated), 
	errorEstimate(0.0)
    {}

    /// Refinement of parent to child1 and child2.
    /// @return true: must this ElementData, else not allowed to delete it
    bool refineElementData(Element* parent, 
			   Element* child1,
			   Element* child2,
			   int elType);

    void coarsenElementData(Element* parent, 
			    Element* thisChild,
			    Element* otherChild,
			    int elTypeParent);

    /// Sets \ref errorEstimate
    inline void setErrorEstimate(int, double est) 
    { 
      errorEstimate = est; 
    }

    /// Returns \ref errorEstimate
    inline double getErrorEstimate(int) 
    { 
      return errorEstimate; 
    }

    /// Implements ElementData::clone().
    virtual ElementData *clone() const 
    {
      // create new estimatable leaf data
      LeafDataEstimatable *newObj = new LeafDataEstimatable(NULL);

      newObj->errorEstimate = errorEstimate;

      // clone decorated element data (=> deep copy)
      newObj->decorated = ElementData::clone();

      // return the clone
      return newObj;
    }

    /// Returns the name of element data type.
    inline std::string getTypeName() const 
    { 
      return "LeafDataEstimatable"; 
    }

    inline int getTypeID() const 
    { 
      return ESTIMATABLE; 
    }

    void serialize(std::ostream& out);

    void deserialize(std::istream& in); 

  private:
    double errorEstimate;
  };

  class LeafDataEstimatableVec : public ElementData,
				 public LeafDataEstimatableInterface
  {
  public:
    class Creator : public CreatorInterface<ElementData>
    {
    public:
      ElementData* create() 
      {
	return new LeafDataEstimatableVec;
      }
    };

    inline bool isOfType(int type) const 
    {
      if (type == ESTIMATABLE) 
	return true;
      return false;
    }

    /// constructor
    LeafDataEstimatableVec(ElementData *decorated = NULL)
      : ElementData(decorated)
    {}

    /// Refinement of parent to child1 and child2.
    bool refineElementData(Element* parent, 
			   Element* child1,
			   Element* child2,
			   int elType);

    void coarsenElementData(Element* parent, 
			    Element* thisChild,
			    Element* otherChild,
			    int elTypeParent);

    /// Sets \ref errorEstimate
    inline void setErrorEstimate(int index, double est) 
    { 
      errorEstimate[index] = est; 
    }

    /// Returns \ref errorEstimate
    inline double getErrorEstimate(int index) 
    { 
      return errorEstimate[index];
    }

    /// Implements ElementData::clone().
    virtual ElementData *clone() const 
    {
      // create new estimatable leaf data
      LeafDataEstimatableVec *newObj = 
	new LeafDataEstimatableVec(NULL);

      newObj->errorEstimate = errorEstimate;

      // clone decorated element data (=> deep copy)
      newObj->decorated = ElementData::clone();

      // return the clone
      return newObj;
    }

    void serialize(std::ostream& out);

    void deserialize(std::istream& in);

    std::string getTypeName() const 
    {
      return "LeafDataEstimatableVec";
    }
  
    inline int getTypeID() const 
    { 
      return ESTIMATABLE; 
    }

  private:
    std::map<int, double> errorEstimate;
  };

  class LeafDataCoarsenableInterface
  {
  public:
    virtual ~LeafDataCoarsenableInterface() {}

    /// Sets \ref coarseningError
    virtual void setCoarseningErrorEstimate(int index, double est) = 0;

    /// Returns \ref coarseningError
    virtual double getCoarseningErrorEstimate(int index) = 0;
  };

  class LeafDataCoarsenable : public ElementData,
			      public LeafDataCoarsenableInterface
  {
  public:
    class Creator : public CreatorInterface<ElementData>
    {
    public:
      ElementData* create() 
      {
	return new LeafDataCoarsenable;
      }
    };

    inline bool isOfType(int type) const 
    {
      if(type == COARSENABLE) 
	return true;
      return false;
    }

    /// constructor
    LeafDataCoarsenable(ElementData *decorated = NULL)
      : ElementData(decorated), 
	coarseningError(0.00)
    {}

    ~LeafDataCoarsenable()
    {}

    /// Refinement of parent to child1 and child2.
    bool refineElementData(Element* parent, 
			   Element* child1,
			   Element* child2,
			   int elType);

    /// Refinement of parent to child1 and child2.
    void coarsenElementData(Element* parent, 
			    Element* thisChild,
			    Element* otherChild,
			    int elTypeParent);

    /// Implements ElementData::clone().
    inline ElementData *clone() const 
    {
      // create new estimatable leaf data
      LeafDataCoarsenable *newObj = new LeafDataCoarsenable(NULL);

      // clone decorated element data (=> deep copy)
      newObj->decorated = ElementData::clone();

      // return the clone
      return newObj;
    }

    /// Sets \ref coarseningError
    virtual void setCoarseningErrorEstimate(int , double est) 
    { 
      coarseningError = est; 
    }

    /// Returns \ref coarseningError
    virtual double getCoarseningErrorEstimate(int) 
    { 
      return coarseningError; 
    }

    void serialize(std::ostream& out);

    void deserialize(std::istream& in);

    std::string getTypeName() const 
    {
      return "LeafDataCoarsenable";
    }
  
    inline int getTypeID() const 
    { 
      return COARSENABLE; 
    }

  private:
    double coarseningError;
  };



  class LeafDataCoarsenableVec : public ElementData,
				 public LeafDataCoarsenableInterface
  {
  public:
    class Creator : public CreatorInterface<ElementData>
    {
    public:
      ElementData* create() 
      {
	return new LeafDataCoarsenableVec;
      }
    };

    inline bool isOfType(int type) const 
    {
      if (type == COARSENABLE) 
	return true;
      return false;
    }

    /// constructor
    LeafDataCoarsenableVec(ElementData *decorated = NULL)
      : ElementData(decorated) 
    {}

    /// Implements ElementData::clone().
    inline ElementData *clone() const 
    {
      // create new estimatable leaf data
      LeafDataCoarsenableVec *newObj = 
	new LeafDataCoarsenableVec(NULL);

      newObj->coarseningError = coarseningError; 

      // clone decorated leaf data (=> deep copy)
      newObj->decorated = ElementData::clone();

      // return the clone
      return newObj;
    }

    /// Refinement of parent to child1 and child2.
    bool refineElementData(Element* parent, 
			   Element* child1,
			   Element* child2,
			   int elType);

    /// Refinement of parent to child1 and child2.
    void coarsenElementData(Element* parent, 
			    Element* thisChild,
			    Element* otherChild,
			    int elTypeParent);

    /// Sets \ref coarseningError
    virtual void setCoarseningErrorEstimate(int index, double est) 
    { 
      coarseningError[index] = est; 
    }

    /// Returns \ref coarseningError
    virtual double getCoarseningErrorEstimate(int index) 
    { 
      return coarseningError[index]; 
    }

    void serialize(std::ostream& out);

    void deserialize(std::istream& in);

    std::string getTypeName() const 
    {
      return "LeafDataCoarsenableVec";
    }
  
    inline int getTypeID() const 
    { 
      return COARSENABLE; 
    }

  private:
    std::map<int, double> coarseningError;
  };



  class LeafDataPeriodic : public ElementData
  {
  public:
    class Creator : public CreatorInterface<ElementData>
    {
    public:
      ElementData* create() 
      {
	return new LeafDataPeriodic;
      }
    };

    inline bool isOfType(int type) const 
    {
      if (type == PERIODIC) 
	return true;
      
      return false;
    }


    class PeriodicInfo : public Serializable
    {
    public:
      PeriodicInfo() 
	: periodicCoords(NULL)
      {}

      PeriodicInfo(int mode,
		   BoundaryType t,
		   int side,
		   const DimVec<WorldVector<double> > *coords);

      virtual ~PeriodicInfo() 
      {
	if (periodicCoords) 
	  delete periodicCoords;
      }

      PeriodicInfo(const PeriodicInfo &rhs);

      void serialize(std::ostream &out);

      void deserialize(std::istream &in);

      int periodicMode;

      BoundaryType type;

      int elementSide;

      DimVec<WorldVector<double> > *periodicCoords;
    };

  public:
    /// constructor
    LeafDataPeriodic(ElementData *decorated = NULL)
      : ElementData(decorated)
    {}

    /// Destructor.
    ~LeafDataPeriodic() {}

    /// Implements LeafData::clone().
    inline ElementData *clone() const 
    {
      LeafDataPeriodic *newObj = new LeafDataPeriodic;
      newObj->decorated = ElementData::clone();
      return newObj;
    }

    inline void addPeriodicInfo(int mode,
				BoundaryType type,
				int side,
				const DimVec<WorldVector<double> > *coords)
    {
      PeriodicInfo periodicInfo(mode, type, side, coords);
      periodicInfoList.push_back(periodicInfo);
    }

    inline std::list<PeriodicInfo>& getInfoList() 
    {
      return periodicInfoList;
    }

    void serialize(std::ostream& out);

    void deserialize(std::istream& in);

    std::string getTypeName() const
    {
      return "LeafDataPeriodic";
    }
  
    inline int getTypeID() const 
    { 
      return PERIODIC; 
    }

    bool refineElementData(Element* parent, Element* child1, Element* child2, 
			   int elType);

  private:
    std::list<PeriodicInfo> periodicInfoList;

    friend class LeafDataPeriodicRefinable;
    friend class LeafDataPeriodicCoarsenable;
  };

}

#endif  // AMDIS_LEAFDATA_H

