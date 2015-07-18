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



/** \file CreatorMap.h */

#ifndef AMDIS_CREATORMAP_H
#define AMDIS_CREATORMAP_H

#include <map>
#include "Global.h"
#include "CreatorInterface.h"

namespace AMDiS {

  /** \ingroup Common
   * \brief
   * A CreatorMap is used to construct objects, which types depends on key words
   * determined at run time. For example the LinearSolverInterfaceMap can create the different
   * solver types depending on the solver parameter of the init file. The benefit
   * of such creator maps is, that you can extend them only by writing an creator
   * class for your own new class and give the creator together with a key word
   * to the map. 
   */
  template <class BaseClass>
  class CreatorMap
  {
  public:
    typedef std::map< std::string, CreatorInterface<BaseClass>* > CreatorMapType;
    
  public:
    /// Adds a new creator together with the given key to the map.
    static void addCreator(std::string key, CreatorInterface<BaseClass>* creator) 
    {
      FUNCNAME("CreatorMap::addCreator()");
      init();
      TEST_EXIT(creatorMap[key] == NULL)
	("there is already a creator for key %s\n",key.c_str());
      creatorMap[key] = creator;
    }
    
    static void addCreator(std::string backend, std::string key, CreatorInterface<BaseClass>* creator)
    {
      addCreator(backend + "_" + key, creator);
    }

    /// Creates a object of the type corresponding to key.
    static CreatorInterface<BaseClass>* getCreator(std::string key,
						   std::string initFileStr) 
    {
      FUNCNAME("CreatorMap::getCreator()");
      init();
      CreatorInterface<BaseClass> *creator = creatorMap[key];
      TEST_EXIT(creator)
	("No creator for key \"%s\" defined in init file for parameter \"%s\"\n", 
	  key.c_str(), initFileStr.c_str());
	
      return creator;
    }

    static void clear();

    static void addDefaultCreators();

  protected:
    /// Constructor is protected because derived maps should be singleton.
    static void init() 
    {
      if (!initialized) {
	initialized = true;
	NullCreator<BaseClass> *nullCreator = new NullCreator<BaseClass>;
	addCreator("0", nullCreator);
	addDefaultCreators();
      }
    }

  protected:
    /// STL map containing the pairs of keys and creators.
    static CreatorMapType creatorMap;

    static bool initialized;
  };

  template <class BaseClass>
  typename CreatorMap<BaseClass>::CreatorMapType CreatorMap<BaseClass>::creatorMap;

  template <class BaseClass>
  bool CreatorMap<BaseClass>::initialized = false;

}

#include "CreatorMap.hh"

#endif
