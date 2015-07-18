/** \file DOFIterator.h */

#pragma once

#include "DOFAdmin.h"
#include "FiniteElemSpace.h"
#include "AMDiS_fwd.h"

namespace AMDiS 
{
  /// Possible types of DOFIterator
  typedef enum {
    USED_DOFS = 0, /**< iterate only used DOFs */
    FREE_DOFS = 1, /**< iterate only free DOFs */
    ALL_DOFS  = 2  /**< iterate all DOFs */
  } DOFIteratorType;


  /** \ingroup DOFAdministration
   * \brief
   * DOFIteratorBase can be the base class of Iterators for DOFIndexed objects
   * or it can be used stand alone. Than it iterates through DOFAdmin's dofFree
   * vector which stores whether a DOF is used or not. If it is used as base
   * class for another Iterator, it provides base functionality, to iterate
   * through the \ref iteratedObject of the sub class. All you have to do is to 
   * override the methods \ref goToBeginOfIteratedObject(), 
   * \ref goToEndOfIteratedObject() and \ref incObjectIterator(). 
   * Normally it is usefull to provide 
   * operator->() and operator*(), to dereference the iterator. But this is not 
   * possible in this base class, because it is template type independent.
   */
  class DOFIteratorBase
  {
  public:
    typedef std::bidirectional_iterator_tag	iterator_category;
  
  public:
    /** \brief
     * Constructs a DOFIteratorBase object of type t which can iterate through 
     * admin's dofFree vector
     */
    DOFIteratorBase(DOFAdmin* admin, DOFIteratorType t) 
      : dofAdmin(admin), 
	dofFree(&(dofAdmin->dofFree)),
	type(t)
    {}

    virtual ~DOFIteratorBase() {}

    /** \brief
     * Resets the iterator to the begin of the iterated object. 
     * Sub classes must
     * implement goToBeginOfIteratedObject() which resets the iterator.
     */
    virtual void reset() 
    {
      position = 0;
      dofFreeIterator = dofFree->begin();
      if (dofFreeIterator == dofFree->end())
	return;

      goToBeginOfIteratedObject();
      if (type != ALL_DOFS)
	if (*dofFreeIterator == (type == USED_DOFS))
	  operator++();
    }

    /** \brief
     * Resets the iterator to the begin of the iterated object. 
     * Sub classes must
     * implement goToBeginOfIteratedObject() which resets the iterator.
     */
    virtual void reset2() 
    {
      position = 0;
      dofFreeIterator = dofFree->begin();
      if (dofFreeIterator == dofFree->end()) 
	return;
      
      goToBeginOfIteratedObject();
      if (type != ALL_DOFS)
	if (*dofFreeIterator == (type == USED_DOFS))
	  operator++();
    }

    /** \brief
     * Prefix operator++.
     * Incrementation depends of the type of the iterator. If type is USED_DOFS,
     * the iterator points to the next used DOF after operator call. If type is
     * FREE_DOFS, it points to the next free DOF and if type is ALL_DOFS, it will
     * point to the next DOF regardless whether it is used or not. Sub classes
     * must implement incObjectIterator() which increments the object
     * iterator.
     */
    inline const DOFIteratorBase& operator++() 
    {
      if (type == ALL_DOFS) {
	incObjectIterator();
	dofFreeIterator++;
	position++;
	return *this;
      }

      if (type == USED_DOFS) {
	if (position >= dofAdmin->getUsedSize()) {
	  position = dofAdmin->getSize();
	  goToEndOfIteratedObject();
	  dofFreeIterator = dofFree->end();
	  return *this;
	}
      }

      do {
	incObjectIterator();
	dofFreeIterator++;
	position++;
      } while ((dofFreeIterator != dofFree->end()) 
	       && (*dofFreeIterator == (type == USED_DOFS)));

      return *this;      
    }

    /// Postfix operator++.
    inline DOFIteratorBase operator++(int) 
    { 
      DOFIteratorBase clone = *this;
      operator++();
      return clone;
    }

    inline const DOFIteratorBase& operator--() 
    {
      if (type == ALL_DOFS) {
	decObjectIterator();
	dofFreeIterator--;
	position--;
	return *this;
      }

      do {
	decObjectIterator();
	dofFreeIterator--;
	position--;
      } while ((dofFreeIterator != dofFree->begin())
	       && (*dofFreeIterator == (type == USED_DOFS)));

      return *this;
    }

    inline DOFIteratorBase operator--(int) 
    {
      DOFIteratorBase clone = *this;
      operator--();
      return clone;
    }
  
    /// Dereferntiation of the \ref dofFreeIterator
    virtual bool isDofFree() 
    {
      return *dofFreeIterator;
    }

    /** \brief
     * Returns whether \ref dofFreeIterator already has reached the end of 
     * \ref dofFree
     */
    bool end() const
    { 
      return (dofFreeIterator == dofFree->end()); 
    }

    bool begin() const
    { 
      return (dofFreeIterator == dofFree->begin()); 
    }

    /// Returns the current position index of this iterator
    int getDOFIndex() const
    { 
      return position; 
    }

  protected:
    /// Override this to enable iteration through the object
    virtual void goToBeginOfIteratedObject() {}

    /// Override this to enable iteration through the object
    virtual void goToEndOfIteratedObject() {}

    /// Override this to enable iteration through the object
    virtual void incObjectIterator() {}
    
    virtual void decObjectIterator() {}

  protected:
    /// DOFAdmin which contains the dofFree vector.
    DOFAdmin *dofAdmin;
    
    /// Current position index.
    int position;

    /// Stores which DOFs are used.
    std::vector<bool> *dofFree; 

    /// Iterator for dofFree.
    std::vector<bool>::iterator dofFreeIterator;

    /// Type of this iterator.
    const DOFIteratorType type;
  };


  /** \ingroup DOFAdministration
   * \brief
   * Implements a DOFIterator for a DOFIndexed<T> object
   */
  template <class T>
  class DOFIterator : public DOFIteratorBase
  {
  public:
    typedef T 					value_type;
    typedef typename DOFIndexed<T>::reference	reference;
    typedef typename DOFIndexed<T>::pointer	pointer;
  
  public:
    /// Constructs a DOFIterator for cont of type t
    DOFIterator(DOFIndexed<T> *obj, DOFIteratorType t) 
      : DOFIteratorBase(dynamic_cast<DOFAdmin*>(obj->getFeSpace()->getAdmin()), t),
	iteratedObject(obj)
    {}

    /// Constructs a DOFIterator for cont of type t
    DOFIterator(DOFAdmin *admin,
		DOFIndexed<T> *obj, 
		DOFIteratorType t) 
      : DOFIteratorBase(admin, t),
	iteratedObject(obj)
    {}

    /// Dereference operator
    reference operator*() 
    { 
      return *it; 
    }

    /// Dereference operator
    pointer operator->() 
    { 
      return &(*it); 
    }

    bool operator!=(const DOFIterator<T>& rhs) const
    {
      if (this->iteratedObject != rhs.iteratedObject) 
	return true;

      if (this->it != rhs.it) 
	return true;

      return false;
    }

    bool operator==(const DOFIterator<T>& rhs) const
    {
      return !(this->operator==(rhs));
    }

  protected:
    /// Implementation of DOFIteratorBase::goToBeginOfIteratedObject()
    void goToBeginOfIteratedObject() 
    { 
      it = iteratedObject->begin(); 
    }

    /// Implementation of DOFIteratorBase::goToEndOfIteratedObject()
    void goToEndOfIteratedObject() 
    { 
      it = iteratedObject->end(); 
    }

    /// Implementation of DOFIteratorBase::incObjectIterator()
    void incObjectIterator() 
    { 
      ++it; 
    }

    /// Implementation of DOFIteratorBase::incObjectIterator()
    void decObjectIterator() 
    { 
      --it; 
    }

  protected:
    /// Object that is iterated
    DOFIndexed<T> *iteratedObject;

    /// Iterator for \ref iteratedObject
    typename std::vector<T>::iterator it;
  };
  
  

  /** \ingroup DOFAdministration
    * \brief
    * Implements a DOFIterator for a const DOFIndexed<T> object
    */
  template <class T>
  class DOFConstIterator : public DOFIteratorBase
  {
  public:
    /// Constructs a DOFIterator for cont of type t
    DOFConstIterator(const DOFIndexed<T> *obj, DOFIteratorType t)
      : DOFIteratorBase(dynamic_cast<DOFAdmin*>(obj->getFeSpace()->getAdmin()), t),
	iteratedObject(obj)
    {}

    /// Constructs a DOFIterator for cont of type t
    DOFConstIterator(DOFAdmin *admin,
		const DOFIndexed<T> *obj,
		DOFIteratorType t)
      : DOFIteratorBase(admin, t),
	iteratedObject(obj)
    {}

    /// Dereference operator
    inline const T& operator*() const
    {
      return *it;
    }

    /// Dereference operator
    inline const T* operator->() const
    {
      return &(*it);
    }

    inline bool operator!=(const DOFIterator<T>& rhs) const
    {
      if (this->iteratedObject != rhs.iteratedObject)
	return true;

      if (this->it != rhs.it)
	return true;

      return false;
    }

    inline bool operator==(const DOFIterator<T>& rhs) const
    {
      return !(this->operator==(rhs));
    }

  protected:
    /// Implementation of DOFIteratorBase::goToBeginOfIteratedObject()
    inline void goToBeginOfIteratedObject()
    {
      it = iteratedObject->begin();
    }

    /// Implementation of DOFIteratorBase::goToEndOfIteratedObject()
    inline void goToEndOfIteratedObject()
    {
      it = iteratedObject->end();
    }

    /// Implementation of DOFIteratorBase::incObjectIterator()
    inline void incObjectIterator()
    {
      ++it;
    }

    /// Implementation of DOFIteratorBase::incObjectIterator()
    inline void decObjectIterator()
    {
      --it;
    }

  protected:
    /// Object that is iterated
    const DOFIndexed<T> *iteratedObject;

    /// Iterator for \ref iteratedObject
    typename std::vector<T>::const_iterator it;
  };



  /** \ingroup DOFAdministration
    * \brief
    * Implements a DOFIterator for a vector of DOFVector<T> objects
    */
  template <class T>
  class DOFVectorIterator : public DOFIteratorBase
  {
  public:
    typedef typename std::vector<T>::iterator VectorIterator;
    
    /// Constructs a DOFIterator for cont of type t
    DOFVectorIterator(std::vector<DOFVector<T>*> &obj, DOFIteratorType t) 
      : DOFIteratorBase(dynamic_cast<DOFAdmin*>(obj[0]->getFeSpace()->getAdmin()), t),
	iteratedObject(obj)
    {
      for (size_t i = 0; i < obj.size(); i++)
	it.push_back(new VectorIterator);
    }

    /// Constructs a DOFIterator for cont of type t
    DOFVectorIterator(DOFAdmin *admin,
		std::vector<DOFVector<T>*> &obj, 
		DOFIteratorType t) 
      : DOFIteratorBase(admin, t),
	iteratedObject(obj)
    {
      for (size_t i = 0; i < obj.size(); i++)
	it.push_back(new VectorIterator);
    }
    
    ~DOFVectorIterator()
    {
      for (size_t i = 0; i < it.size(); i++)
	delete it[i];
    }

    /// Dereference operator
    inline std::vector<T> operator*() 
    { 
      std::vector<T> result(it.size());
      for (size_t i = 0; i < it.size(); i++)
	result[i] = *(*it[i]);
      return result;
    }

    /// Dereference operator
    inline T* operator->() 
    { 
      throw std::runtime_error("operator-> not available for DOFVectorIterator!");
      return &(*(*it[0])); 
    }

    inline bool operator!=(const DOFVectorIterator<T>& rhs) 
    {
      if (this->iteratedObject != rhs.iteratedObject) 
	return true;

      if (this->it != rhs.it) 
	return true;

      return false;
    }

    inline bool operator==(const DOFVectorIterator<T>& rhs) 
    {
      return !(this->operator==(rhs));
    }

  protected:
    /// Implementation of DOFIteratorBase::goToBeginOfIteratedObject()
    inline void goToBeginOfIteratedObject() 
    { 
      for (size_t i = 0; i < it.size(); i++)
	*(it[i]) = iteratedObject[i]->begin(); 
    }

    /// Implementation of DOFIteratorBase::goToEndOfIteratedObject()
    inline void goToEndOfIteratedObject() 
    { 
      for (size_t i = 0; i < it.size(); i++)
	*(it[i]) = iteratedObject[i]->end(); 
    }

    /// Implementation of DOFIteratorBase::incObjectIterator()
    inline void incObjectIterator() 
    { 
      for (size_t i = 0; i < it.size(); i++)
	++(*(it[i])); 
    }

    /// Implementation of DOFIteratorBase::incObjectIterator()
    inline void decObjectIterator() 
    { 
      for (size_t i = 0; i < it.size(); i++)
	--(*(it[i]));
    }

  protected:
    /// Object that is iterated
    std::vector<DOFVector<T>*> iteratedObject;

    /// Iterator for \ref iteratedObject
    std::vector<typename std::vector<T>::iterator*> it;
  };
  
} // end namespace AMDiS
