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



/** \file Serializable.h */

#ifndef AMDIS_SERIALIZABLE_H
#define AMDIS_SERIALIZABLE_H

#include <iostream>
#include <string>

namespace AMDiS {

  /** \brief
   * Interface for JAVA-like serialization of objects. Serializable objects can be
   * written to an out stream and read from an in stream. So i.e. they can be stored
   * on hard disk.
   */
  class Serializable
  {
  public:
    /// Streams the object to output stream out.
    virtual void serialize(std::ostream &out) = 0;

    /// Reads the object from input stream in.
    virtual void deserialize(std::istream &in) = 0;

    /// Returns the type name for this serializable object.
    virtual std::string getTypeName() const 
    { 
      return ""; 
    }

    virtual ~Serializable() {}
  };

}

#endif
