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



/** \file ElInfo.h */

#ifndef AMDIS_ELINFOSTACK_H
#define AMDIS_ELINFOSTACK_H

#include <vector>
#include "ElInfo.h"
#include "AMDiS_fwd.h"

namespace AMDiS {

  /** \ingroup Traverse
   * \brief 
   * Stores a stack of ElInfo object. Is used by meshes for recursive mesh
   * traverse. The use of a stack is cheaper than allocating the ElInfo objects
   * at every recursive step.
   */
  class ElInfoStack
  {
  public:
    /// Constructer, creates the stack.
    ElInfoStack(Mesh *mesh);
    
    /// Destructor, deletes all ElInfos on the stack.
    ~ElInfoStack();

    /// Get a new element from the stack an increase the stack position.
    ElInfo* getNextElement();

    /// Decrease the stack position.
    void getBackElement();

    /// Returns a pointer to the currently used element of the stack.
    ElInfo* getCurrentElement();

  protected:
    /// The mesh on which the traverse is done.
    Mesh *mesh_;

    /// The stack of pointers to ElInfo objects.
    std::vector<ElInfo*> elInfoStack_;

    /// Current position (depth) of the recursive mesh traverse.
    int stackPosition_;
  };
  
} // end namespace AMDiS

#endif  // AMDIS_ELINFOSTACK_H
