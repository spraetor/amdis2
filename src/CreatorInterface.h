/** \file CreatorInterface.h */

#pragma once

#include "DOFMatrix.h"

namespace AMDiS
{

  /** \ingroup Common
   * \brief
   * Interface for the implementation of the factory method pattern.
   * The creation of an object of a sub class of BaseClass is deligated
   * to a corresponding sub class of Creator<BaseClass>. So it is possible to
   * manage a CreatorMap, which can be extended at run-time. An example is
   * the LinearSolverInterfaceMap: If you write your own LinearSolverInterface sub class and a
   * corresponding Creator<LinearSolverInterface<T> >, you can add the creator together
   * with a key string to the LinearSolverInterfaceMap. Then you can create an LinearSolverInterface
   * depending of a key string read from the init file, which can also be
   * your own new solver.
   */
  template <class BaseClass>
  class CreatorInterface
  {
  public:
    virtual ~CreatorInterface() {}

    /** \brief
     * Must be implemented by sub classes of CreatorInterface.
     * Creates a new instance of the sub class of BaseClass.
     */
    virtual BaseClass* create() = 0;

    virtual BaseClass* create(const DOFMatrix::base_matrix_type& A)
    {
      return 0;
    }

    /// Can be implemented by sub classes.
    virtual void free(BaseClass*) {}

    ///
    virtual bool isNullCreator()
    {
      return false;
    }
  };

  /** \brief
   * A Creator which creates no object and returns NULL instead.
   * Used together with the key word 'no' in CreatorMap.
   */
  template <class BaseClass>
  class NullCreator : public CreatorInterface<BaseClass>
  {
    /// Creates no object.
    virtual BaseClass* create() override
    {
      return NULL;
    }

    /// Implementation of \ref CreatorInterface::isNullCreator()
    virtual bool isNullCreator() override
    {
      return true;
    }
  };


  /**
   * \ingroup Common
   *
   * \brief
   * Interface for creators with name.
   */
  template <class BaseClass>
  class CreatorInterfaceName : public CreatorInterface<BaseClass>
  {
  public:
    /// Sets \ref name
    void setName(std::string str)
    {
      name = str;
    }

  protected:
    std::string name;
  };

} // end namespace AMDiS
