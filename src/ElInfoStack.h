/** \file ElInfo.h */

#pragma once

#include <vector>

#include "AMDiS_fwd.h"

namespace AMDiS
{

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
    ElInfoStack(Mesh* mesh);

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
    Mesh* mesh_;

    /// The stack of pointers to ElInfo objects.
    std::vector<ElInfo*> elInfoStack_;

    /// Current position (depth) of the recursive mesh traverse.
    int stackPosition_;
  };

} // end namespace AMDiS
