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



/** \file TimedObject.h */

#ifndef AMDIS_TIMEDOBJECT_H
#define AMDIS_TIMEDOBJECT_H

namespace AMDiS
{

  /** \brief
   * This class can be used as base class for time dependent objects where
   * different objects refer to the same time. Herefore a pointer to
   * a double value is stored, pointing to a time value, which can be
   * managed in one central object, maybe the problem class.
   */
  class TimedObject
  {
  public:
    /// Constructor.
    TimedObject()
      : timePtr(NULL)
    {}

    /// Sets the time pointer.
    inline void setTimePtr(double* ptr)
    {
      timePtr = ptr;
    }

    /// Returns the time pointer.
    inline double* getTimePtr()
    {
      return timePtr;
    }

  protected:
    /// Pointer to the externally managed time value.
    double* timePtr;
  };

}

#endif // AMDIS_TIMEDOBJECT_H
