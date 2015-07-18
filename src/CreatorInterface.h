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



/** \file CreatorInterface.h */

#ifndef AMDIS_CREATORINTERFACE_H
#define AMDIS_CREATORINTERFACE_H

#include "DOFMatrix.h"

namespace AMDiS {

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
  template<typename BaseClass>
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
    virtual void free(BaseClass *) {}

    ///
    virtual bool isNullCreator() 
    { 
      return false; 
    }
  };

  /** \brief
   * A Creator which creates no object abd returns NULL instead.
   * Used together with the key word 'no' in CreatorMap.
   */
  template<typename BaseClass>
  class NullCreator : public CreatorInterface<BaseClass>
  {
    /// Creates no object.
    BaseClass* create() 
    {
      return NULL; 
    }

    virtual ~NullCreator() {}

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
  template< typename BaseClass >
  class CreatorInterfaceName : public CreatorInterface<BaseClass>
  { 
  public:
    virtual ~CreatorInterfaceName() {}

    /// Sets \ref name
    void setName(std::string str) 
    { 
      name = str; 
    }

  protected:
    std::string name;
  };
}

#endif
