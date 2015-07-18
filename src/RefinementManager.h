/** \file RefinementManager.h */

#pragma once

#include "Global.h"
#include "Flag.h"
#include "AMDiS_fwd.h"

namespace AMDiS 
{
  /** \ingroup Adaption 
   * \brief
   * Base class of RefinementManager1d, RefinementManager2d, RefinementManager3d.
   * A RefinementManager contains all functionality to perform refinement
   * operations on the mesh.
   */
  class RefinementManager
  {
  public:
    /// Constructs a RefinementManager which belongs to aMesh 
    RefinementManager()
      : mesh(NULL), 
	newCoords(false),
	stack(NULL)
    {}

    /// Destructor
    virtual ~RefinementManager() {}

    /** \brief
     * Generates new coordinates on curved boundaries. Can be overriden by
     * sub classes if used.
     */
    virtual void setNewCoords(int macroEl = -1) 
    {
      FUNCNAME("RefinementManager::setNewCoords");
      ERROR_EXIT("called for base class!\n");
    }

    /** \brief
     * Fulfills the refinement for all positive marked elements of the mesh.
     * The default implementation of the base class uses \ref refineFunction
     * for the refinement of a single element while traversing the mesh.
     * Sub classes can overload refineMesh and/or refineFunction to implement
     * their own refinement routines. 
     */
    virtual Flag refineMesh(Mesh *aMesh);

    void refineMacroElement(Mesh *amesh, int macroElIndex);

    /// All elements of the mesh will be refined.
    Flag globalRefine(Mesh *aMesh, int mark);

    /// Set \ref newCoords
    void newCoord(bool nc) 
    { 
      newCoords = nc; 
    }

    bool newCoord() 
    { 
      return newCoords; 
    }

    /** \brief
     * Implements the refinement of el_info->el. Can be overriden by sub
     * classes if used.
     */
    virtual ElInfo* refineFunction(ElInfo*) 
    {
      FUNCNAME("RefinementManager::refineFunction()");
      ERROR_EXIT("called for base class!\n");
      return NULL;
    }

    void setMesh(Mesh *m) 
    {
      mesh = m;
    }

    void setStack(TraverseStack *s) 
    {
      stack = s; 
    }

    TraverseStack *getStack() 
    { 
      return stack; 
    }

  protected:
    /// The Mesh to be refined
    Mesh *mesh;

    /// Number of new vertices on a boundary edge
    bool newCoords;

    /// Still more refinement to do?
    static bool doMoreRecursiveRefine;

    /// Number of DOFVectors which must be interpolated during refinement
    static int callRefineInterpol;
  
    /// Used for non recursive traversal
    TraverseStack* stack;
  };

} // end namespace AMDiS

#include "RefinementManager1d.h"
#include "RefinementManager2d.h"
#include "RefinementManager3d.h"
