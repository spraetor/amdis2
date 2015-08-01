/** \file MacroInfo.h */

#pragma once

#include <deque>

#include <AMDiS_fwd.h>
#include <AMDiS_base.h>
#include <MatrixVector_fwd.h>
#include <Global.h>

namespace AMDiS 
{

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
    
    bool isInitialized() const
    {
      return initialized;
    }

  protected:
    /********************************************************************************/
    /*  read_indices()  reads dim + 1 indices from  file  into  id[0 - dim],        */
    /*    returns true if dim + 1 inputs arguments could be read successfully by    */
    /*    fscanf(), else false                                                      */
    /********************************************************************************/  
    template <class VectorType>
    bool read_indices(FILE *file, VectorType& id);

    bool initialized;
  };

} // end namespace AMDiS
