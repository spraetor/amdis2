/** \file LeafData.h */

#pragma once

#include <list>
#include <map>

#include "FixVec.h"
#include "ElementData.h"
#include "Boundary.h"

namespace AMDiS
{

  class LeafDataEstimatableInterface
  {
  public:
    virtual ~LeafDataEstimatableInterface() {};
    virtual void setErrorEstimate(int, double) = 0;
    virtual double getErrorEstimate(int) const = 0;
  };


  class LeafDataEstimatable : public ElementData,
    public LeafDataEstimatableInterface
  {
  public:
    virtual bool isOfType(int type) const override
    {
      if (type == ESTIMATABLE)
        return true;

      return false;
    }

    struct Creator : public CreatorInterface<ElementData>
    {
      virtual ElementData* create() override
      {
        return new LeafDataEstimatable;
      }
    };

    /// constructor
    LeafDataEstimatable(ElementData* decorated = NULL)
      : ElementData(decorated),
        errorEstimate(0.0)
    {}

    /// Refinement of parent to child1 and child2.
    /// @return true: must this ElementData, else not allowed to delete it
    virtual bool refineElementData(Element* parent,
                                   Element* child1,
                                   Element* child2,
                                   int elType) override;

    virtual void coarsenElementData(Element* parent,
                                    Element* thisChild,
                                    Element* otherChild,
                                    int elTypeParent) override;

    /// Sets \ref errorEstimate
    virtual void setErrorEstimate(int, double est) override
    {
      errorEstimate = est;
    }

    /// Returns \ref errorEstimate
    virtual double getErrorEstimate(int) const override
    {
      return errorEstimate;
    }

    /// Implements ElementData::clone().
    virtual ElementData* clone() const override
    {
      // create new estimatable leaf data
      LeafDataEstimatable* newObj = new LeafDataEstimatable(NULL);

      newObj->errorEstimate = errorEstimate;

      // clone decorated element data (=> deep copy)
      newObj->decorated = ElementData::clone();

      // return the clone
      return newObj;
    }

    virtual int getTypeID() const override
    {
      return ESTIMATABLE;
    }

  private:
    double errorEstimate;
  };


  class LeafDataEstimatableVec : public ElementData,
    public LeafDataEstimatableInterface
  {
  public:
    struct Creator : public CreatorInterface<ElementData>
    {
      virtual ElementData* create() override
      {
        return new LeafDataEstimatableVec;
      }
    };

    virtual bool isOfType(int type) const override
    {
      if (type == ESTIMATABLE)
        return true;
      return false;
    }

    /// constructor
    LeafDataEstimatableVec(ElementData* decorated = NULL)
      : ElementData(decorated)
    {}

    /// Refinement of parent to child1 and child2.
    virtual bool refineElementData(Element* parent,
                                   Element* child1,
                                   Element* child2,
                                   int elType) override;

    virtual void coarsenElementData(Element* parent,
                                    Element* thisChild,
                                    Element* otherChild,
                                    int elTypeParent) override;

    /// Sets \ref errorEstimate
    virtual void setErrorEstimate(int index, double est) override
    {
      errorEstimate[index] = est;
    }

    /// Returns \ref errorEstimate
    virtual double getErrorEstimate(int index) const override
    {
      return errorEstimate.at(index);
    }

    /// Implements ElementData::clone().
    virtual ElementData* clone() const override
    {
      // create new estimatable leaf data
      LeafDataEstimatableVec* newObj =
        new LeafDataEstimatableVec(NULL);

      newObj->errorEstimate = errorEstimate;

      // clone decorated element data (=> deep copy)
      newObj->decorated = ElementData::clone();

      // return the clone
      return newObj;
    }

    virtual int getTypeID() const override
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
    virtual double getCoarseningErrorEstimate(int index) const = 0;
  };


  class LeafDataCoarsenable : public ElementData,
    public LeafDataCoarsenableInterface
  {
  public:
    struct Creator : public CreatorInterface<ElementData>
    {
      virtual ElementData* create() override
      {
        return new LeafDataCoarsenable;
      }
    };

    virtual bool isOfType(int type) const override
    {
      if(type == COARSENABLE)
        return true;
      return false;
    }

    /// constructor
    LeafDataCoarsenable(ElementData* decorated = NULL)
      : ElementData(decorated),
        coarseningError(0.00)
    {}

    /// Refinement of parent to child1 and child2.
    virtual bool refineElementData(Element* parent,
                                   Element* child1,
                                   Element* child2,
                                   int elType) override;

    /// Refinement of parent to child1 and child2.
    virtual void coarsenElementData(Element* parent,
                                    Element* thisChild,
                                    Element* otherChild,
                                    int elTypeParent) override;

    /// Implements ElementData::clone().
    virtual ElementData* clone() const override
    {
      // create new estimatable leaf data
      LeafDataCoarsenable* newObj = new LeafDataCoarsenable(NULL);

      // clone decorated element data (=> deep copy)
      newObj->decorated = ElementData::clone();

      // return the clone
      return newObj;
    }

    /// Sets \ref coarseningError
    virtual void setCoarseningErrorEstimate(int , double est) override
    {
      coarseningError = est;
    }

    /// Returns \ref coarseningError
    virtual double getCoarseningErrorEstimate(int) const override
    {
      return coarseningError;
    }

    virtual int getTypeID() const override
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
    struct Creator : public CreatorInterface<ElementData>
    {
      virtual ElementData* create() override
      {
        return new LeafDataCoarsenableVec;
      }
    };

    virtual bool isOfType(int type) const override
    {
      if (type == COARSENABLE)
        return true;
      return false;
    }

    /// constructor
    LeafDataCoarsenableVec(ElementData* decorated = NULL)
      : ElementData(decorated)
    {}

    /// Implements ElementData::clone().
    virtual ElementData* clone() const override
    {
      // create new estimatable leaf data
      LeafDataCoarsenableVec* newObj =
        new LeafDataCoarsenableVec(NULL);

      newObj->coarseningError = coarseningError;

      // clone decorated leaf data (=> deep copy)
      newObj->decorated = ElementData::clone();

      // return the clone
      return newObj;
    }

    /// Refinement of parent to child1 and child2.
    virtual bool refineElementData(Element* parent,
                                   Element* child1,
                                   Element* child2,
                                   int elType) override;

    /// Refinement of parent to child1 and child2.
    virtual void coarsenElementData(Element* parent,
                                    Element* thisChild,
                                    Element* otherChild,
                                    int elTypeParent) override;

    /// Sets \ref coarseningError
    virtual void setCoarseningErrorEstimate(int index, double est) override
    {
      coarseningError[index] = est;
    }

    /// Returns \ref coarseningError
    virtual double getCoarseningErrorEstimate(int index) const override
    {
      return coarseningError.at(index);
    }

    virtual int getTypeID() const override
    {
      return COARSENABLE;
    }

  private:
    std::map<int, double> coarseningError;
  };


  class LeafDataPeriodic : public ElementData
  {
  public:
    struct Creator : public CreatorInterface<ElementData>
    {
      virtual ElementData* create() override
      {
        return new LeafDataPeriodic;
      }
    };

    virtual bool isOfType(int type) const override
    {
      if (type == PERIODIC)
        return true;

      return false;
    }


    class PeriodicInfo
    {
    public:
      PeriodicInfo()
        : periodicCoords(NULL)
      { }

      PeriodicInfo(int mode,
                   BoundaryType t,
                   int side,
                   const DimVec<WorldVector<double>>* coords);

      virtual ~PeriodicInfo()
      {
        if (periodicCoords)
          delete periodicCoords;
      }

      PeriodicInfo(const PeriodicInfo& rhs);

      int periodicMode;

      BoundaryType type;

      int elementSide;

      DimVec<WorldVector<double>>* periodicCoords;
    };

  public:
    /// constructor
    LeafDataPeriodic(ElementData* decorated = NULL)
      : ElementData(decorated)
    {}

    /// Implements LeafData::clone().
    virtual ElementData* clone() const override
    {
      LeafDataPeriodic* newObj = new LeafDataPeriodic;
      newObj->decorated = ElementData::clone();
      return newObj;
    }

    void addPeriodicInfo(int mode,
                         BoundaryType type,
                         int side,
                         const DimVec<WorldVector<double>>* coords)
    {
      PeriodicInfo periodicInfo(mode, type, side, coords);
      periodicInfoList.push_back(periodicInfo);
    }

    std::list<PeriodicInfo>& getInfoList()
    {
      return periodicInfoList;
    }

    virtual int getTypeID() const override
    {
      return PERIODIC;
    }

    virtual bool refineElementData(Element* parent, Element* child1, Element* child2,
                                   int elType) override;

  private:
    std::list<PeriodicInfo> periodicInfoList;

    friend class LeafDataPeriodicRefinable;
    friend class LeafDataPeriodicCoarsenable;
  };

} // end namespace AMDiS
