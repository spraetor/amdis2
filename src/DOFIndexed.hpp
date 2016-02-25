/** \file DOFIndexed.h */

#pragma once

#include <cassert>
#include <vector>

#include "AMDiS_fwd.hpp"
#include "AMDiS_base.hpp"

namespace AMDiS
{
  /** \ingroup DOFAdministration
   * \brief
   * Interface for objects that stores information indexed by DOF indices
   * (like DOFVector or DOFMatrix). This interface is template type independent,
   * so a DOFAdmin can handle a single list of DOFIndexedBase objects.
   */
  class DOFIndexedBase
  {
  public:
    DOFIndexedBase()
      : coarsenOperation(COARSE_INTERPOL),
        refineOperation(REFINE_INTERPOL)
    {}

    virtual ~DOFIndexedBase() {}

    /// Returns the actual size. Must be overriden by sub classes
    virtual int getSize() const = 0;

    /// Resizes the object to size. Must be overriden by sub classes
    virtual void resize(int size) = 0;

    /// Check that object's size is equal to \p size or 0.  If object's size is 0 change it to \p size.
    virtual void checked_change_dim(int size)
    {
      assert(getSize() == 0 || getSize() == size);
      if (getSize() == 0)
        resize(size);
    }

    /// Used by DOFAdmin::compress. Must be overriden by sub classes
    virtual void compressDOFIndexed(int first, int last,
                                    std::vector<DegreeOfFreedom>& newDOF) = 0;

    /// Performs needed action when a DOF index is freed. Can be overriden in
    /// sub classes. The default behavior is to do nothing.
    virtual void freeDOFContent(DegreeOfFreedom) {}

    /// Interpolation after refinement. Can be overriden in subclasses.
    /// The default behavior is to do nothing.
    virtual void refineInterpol(RCNeighbourList&, int) {}

    /// Restriction after coarsening. Can be overriden in subclasses.
    /// The default behavior is to do nothing.
    virtual void coarseRestrict(RCNeighbourList&, int) {}

    /// Returns the finite element space of this DOFIndexed object. Must be
    /// overriden in sub classes.
    virtual const FiniteElemSpace* getFeSpace() const = 0;

    /// Sets \ref coarsenOperation to op
    void setCoarsenOperation(RefineCoarsenOperation op)
    {
      coarsenOperation = op;
    }

    /// Returns \ref coarsenOperation
    RefineCoarsenOperation getCoarsenOperation() const
    {
      return coarsenOperation;
    }

    /// Sets \ref refineOperation to op
    void setRefineOperation(RefineCoarsenOperation op)
    {
      refineOperation = op;
    }

    /// Returns \ref refineOperation
    RefineCoarsenOperation getRefineOperation() const
    {
      return refineOperation;
    }

  protected:
    /// Specifies what operation should be performed after coarsening
    RefineCoarsenOperation coarsenOperation;
    RefineCoarsenOperation refineOperation;
  };


  /** \ingroup DOFAdministration
   * \brief
   * Templated interface for DOFIndexed objects.
   */
  template <class T>
  class DOFIndexed : public DOFIndexedBase
  {
  public: // typedefs
    typedef T 			value_type;
    typedef DegreeOfFreedom	size_type;
    typedef value_type&		reference;
    typedef value_type const&	const_reference;
    typedef value_type*		pointer;
    typedef value_type const*	const_pointer;

    typedef typename std::vector<value_type>::iterator		iterator;
    typedef typename std::vector<value_type>::const_iterator	const_iterator;

  public:
    virtual ~DOFIndexed() {}

    /// Returns iterator to the begin of container
    virtual iterator begin() = 0;

    /// Returns iterator to the end of container
    virtual iterator end() = 0;

    /// Returns iterator to the begin of container
    virtual const_iterator begin() const = 0;

    /// Returns iterator to the end of container
    virtual const_iterator end() const = 0;

    /// Returns container element at index i
    virtual reference operator[](size_type i) = 0;

    /// Returns container element at index i
    virtual const_reference operator[](size_type i) const = 0;
  };


  template<>
  class DOFIndexed<bool> : public DOFIndexedBase
  {
  public: // typedefs
    typedef bool 				value_type;
    typedef DegreeOfFreedom			size_type;
    typedef std::vector<bool>::reference	reference;
    typedef std::vector<bool>::const_reference	const_reference;
    typedef std::vector<bool>::pointer		pointer;
    typedef std::vector<bool>::const_pointer	const_pointer;

    typedef std::vector<bool>::iterator		iterator;
    typedef std::vector<bool>::const_iterator	const_iterator;

  public:
    virtual ~DOFIndexed() {}

    /// Returns iterator to the begin of container
    virtual iterator begin() = 0;

    /// Returns iterator to the end of container
    virtual iterator end() = 0;

    /// Returns iterator to the begin of container
    virtual const_iterator begin() const = 0;

    /// Returns iterator to the end of container
    virtual const_iterator end() const = 0;

    /// Returns container element at index i
    virtual reference operator[](size_type i) = 0;

    /// Returns container element at index i
    virtual const_reference operator[](size_type i) const = 0;
  };

} // end namespace AMDiS