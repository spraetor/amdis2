#pragma once

// AMDiS includes
#include "AMDiS_fwd.hpp"
#include "Flag.hpp"

namespace AMDiS
{
  /** \ingroup Adaption
   * \brief
   * Base class of CoarseningManager1d, CoarseningManager2d, CoarseningManager3d.
   * A CoarseningManager contains all functionality to perform coarsening
   * operations on the mesh.
   */
  class CoarseningManager
  {
  public:
    /// Constructs a CoarseningManager which belongs to aMesh
    CoarseningManager()
      : mesh(NULL),
        stack(NULL),
        doMore(0)
    {}

    /// destructor
    virtual ~CoarseningManager() {}

    /// Returns the Mesh the CoarseningManager belongs to.
    Mesh* getMesh() const
    {
      return mesh;
    }

    /** \brief
     * Tries to coarsen every element of mesh at least mark times. First
     * all elements are marked for coarsening and then coarsenMesh will
     * be called.
     */
    Flag globalCoarsen(Mesh* aMesh, int mark);

    /** \brief
     * Traversal routine for recursiv coarsening of a triangulation. It has
     * a default definition in CoarseningManager but it can be overriden
     * by sub classes (like in CoarseningManager1d), if another implementation
     * is needed.
     */
    virtual Flag coarsenMesh(Mesh* aMesh);


  protected:
    /// Defines the way how one element of the mesh is coarsen.
    virtual void coarsenFunction(ElInfo*) {}

    /** \brief
     *  Propagate coarsening information over the whole hierarchy
     *  by POSTORDER traversal of the hierarchy tree.
     *  leaves:      'increment' coarsening mark,
     *  inner nodes: set coarsening mark to
     *               min(0,child[0].mark+1,child[1].mark+1)
     */
    void spreadCoarsenMark();

    /// Resets the element marks
    void cleanUpAfterCoarsen();

  protected:
    /// The Mesh this CoarseningManager belongs to.
    Mesh* mesh;

    /// Used for non recursive mesh traversal.
    TraverseStack* stack;

    /// Spezifies whether the coarsening operation is still in progress
    bool doMore;

    /// Spezifies how many DOFVectors should restricted while coarsening
    int callCoarseRestrict;

    friend class RCNeighbourList;
  };

} // end namespace AMDiS

#include "CoarseningManager1d.hpp"
#include "CoarseningManager2d.hpp"
#include "CoarseningManager3d.hpp"
