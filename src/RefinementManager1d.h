/** \file RefinementManager1d.h */

#pragma once

namespace AMDiS {

  /** \ingroup Adaption
   * \brief
   * Implements a RefinementManager for 1-dimensional meshes.
   */
  class RefinementManager1d : public RefinementManager
  {
  public:
    /// Calls base class constructor.
    RefinementManager1d()
      : RefinementManager()
    {}

    /// destructor 
    virtual ~RefinementManager1d() {}

    /// Implements RefinementManager::refineMesh.
    Flag refineMesh(Mesh *aMesh);

    /// Implements RefinementManager::setNewCoords
    void setNewCoords(int macroEl = -1);

  protected:
    /// Used by refineMesh while mesh traversal
    void recursiveRefineFunction(ElInfo* el_info);

    /// Used by \ref setNewCoords
    void newCoordsFct(ElInfo *el_info);
  };

} // end namespace AMDiS
