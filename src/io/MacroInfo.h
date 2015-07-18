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



/** \file MacroInfo.h */

#ifndef AMDIS_MACROINFO_H
#define AMDIS_MACROINFO_H

#include <deque>
#include "AMDiS_fwd.h"
#include "Global.h"

namespace AMDiS {

  /** \ingroup Input
   * \brief
   * Used for reading a macro triangulation
   */
  class MacroInfo
  {
  public:
    /// Pointer to the Mesh
    Mesh *mesh;

    /// list of macro elements
    std::deque<MacroElement*> mel;

    /// vector of all vertex dofs
    DegreeOfFreedom **dof;

    /// coords[j][k]: kth coordinate of global vertex j
    WorldVector<double> *coords;

    /// mel_vertex[i][k]: global index of kth vertex of element i
    int **mel_vertex;

    /// true, if neighbour information is in macro file
    bool neigh_set;
    
    /// true, if inverse neighbour information is in macro file (see MacroReader.h)
    bool neigh_inverse_set;

    /// true, if boundary information is in macro file
    bool bound_set;
    
    int nElements;
    
    int nVertices;

  public:
    MacroInfo() : initialized(false) {};
    /** \brief
     * Reads macro triangulation from ascii file in AMDiS format.
     * Fills MacroInfo structure.
     * Called by Mesh::readMacro(), fills missing information  
     */
    void readAMDiSMacro(std::string filename, Mesh* mesh);

    /// Fills MacroInfo structure and some pointers in mesh 
    void fill(Mesh *mesh, int nElements, int nVertices);  

    /// Frees memory of MacroInfo
    void clear();

    /** \brief
     * Sets the boundary of all edges/faces with no neigbour to a straight  
     * line/face with dirichlet boundary type
     */
    void dirichletBoundary();

    void fillBoundaryInfo(Mesh *mesh);
    
    bool isInitialized()
    {
      return initialized;
    }

  protected:
    /// Reads indices from macro file
    int read_indices(FILE *file, Vector<int> &id);

    bool initialized;
  };

}

#endif
