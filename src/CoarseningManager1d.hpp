#pragma once

// AMDiS includes
#include "CoarseningManager.hpp"
#include "Log.hpp"

namespace AMDiS
{
  /** \ingroup Adaption
   * \brief
   * Implements a CoarseningManager for 1-dimensional meshes.
   */
  class CoarseningManager1d : public CoarseningManager
  {
  public:
    /// Calls base class constructor and checks dimension of mesh.
    CoarseningManager1d()
      : CoarseningManager()
    {}

    /** \brief
     * Overloads CoarseningManager::coarsenMesh. In 1d a simple recursive
     * coarsening algorithm is implemented which doesn't need coarsenFunction.
     */
    Flag coarsenMesh(Mesh* aMesh);

  protected:
    /// Not needed in this sub class
    void coarsenFunction(ElInfo*)
    {
      FUNCNAME("CoarseningManager1d::coarsenFunction");
      ERROR_EXIT("not used for dim = 1");
    }

    /// Needed instead of coarsenFunction in 1d.
    int coarsenRecursive(Line* parent);

  };

} // end namespace AMDiS
