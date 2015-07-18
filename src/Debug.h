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



/** \file Debug.h */

#ifndef AMDIS_DEBUG_H
#define AMDIS_DEBUG_H

#include <set>
#include "AMDiS_fwd.h"
#include "Global.h"
#include "FixVec.h"

namespace AMDiS {
  
  namespace debug {

    struct DofPtrSortFct {
      bool operator() (const DegreeOfFreedom *dof0, const DegreeOfFreedom *dof1) 
      {
	return (*dof0 < *dof1);
      }
    };

    typedef std::map<int, DofContainer> ElementIdxToDofs;
    typedef std::map<int, FixVec<WorldVector<double>, VERTEX> > ElementIdxToCoords;

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    void writeLocalElementDofs(int rank, 
			       int elIdx, 
			       const FiniteElemSpace *feSpace);
    
    void writeMesh(const FiniteElemSpace *feSpace, 
		   int rank, 
		   std::string filename);
    
    /** \brief
     * Writes a vtu file with the mesh, where all DOFs are set to zero, and only
     * one given DOF is set to one. This can be used to easily identify DOFs in
     * a mesh.
     *
     * \param[in]  rank     If set to -1, the vtu files are written on all ranks.
     *                      Otherwise, only on the given rank the mesh is written. 
     * \param[in]  dof      Defines the DOF, which value is set to one in the mesh file.
     * \param[in]  feSpace  The FE space to be used.
     */
    void writeDofMesh(int rank, 
		      DegreeOfFreedom dof, 
		      const FiniteElemSpace *feSpace);
#endif

    /** \brief
     * Create a vtu file with name 'dofindex.vtu'. All nodes in the mesh are colored
     * by the global DOF index.
     *
     * \param[in]  feSpace   The FE space to be used.
     * \param[in]  filename  Name of the VTU file
     */
    void writeDofIndexMesh(const FiniteElemSpace *feSpace, 
			   std::string filename = "dofindex.vtu");

    void colorEdgeInMesh(const FiniteElemSpace *feSpace,
			 Element *el, 
			 int localEdgeNo, 
			 std::string filename);

    /** \brief
     * Creates a vtu file where all elements in the mesh are colored by the
     * global element indices.
     *
     * \param[in]  feSpace   The FE space to be used.
     * \param[in]  filename  Name of the file.
     * \param[in]  level     If level is -1, all leaf elements will be put to
     *                       the output file, otherwise the elements with the
     *                       given level.
     */
    void writeElementIndexMesh(Mesh *mesh, 
			       std::string filename, 
			       int level = -1);

    /** \brief
     * Creates a vtu file where all elements in the mesh are colored by the
     * macro element indices.
     *
     * \param[in]  feSpace   The FE space to be used.
     * \param[in]  filename  Name of the file.
     */
    void writeMacroElementIndexMesh(Mesh *mesh,
				    std::string filename);

    void highlightElementIndexMesh(Mesh *mesh, int idx, std::string filename);

    void colorMeshByMacroIndex(Mesh *mesh, std::string filename);

    void colorDofVectorByLocalElementDofs(DOFVector<double>& vec, Element *el);
    
    bool colorDofVectorByLocalElementDofs(DOFVector<double>& vec, 
					  Mesh *mesh, 
					  int elIndex);
    
    Element* getDofIndexElement(const FiniteElemSpace *feSpace, 
				DegreeOfFreedom dof);
    
    Element* getLevel0ParentElement(Mesh *mesh, Element *el);

    Element* getLevel0ParentElement(Mesh *mesh, int elIndex);

    Element* getParentElement(Mesh *mesh, Element *el);

    Element* getParentElement(Mesh *mesh, int elIndex);

    Element* getElement(Mesh *mesh, int elIndex);

    void printElementInfo(Element *el);

    void printElementCoords(const FiniteElemSpace *feSpace, Element *el);
    
    void printInfoByDof(const FiniteElemSpace *feSpace, DegreeOfFreedom dof);

    void printMatValuesStatistics(Matrix<DOFMatrix*> *mat);

    void printAllDofCoords(const FiniteElemSpace *feSpace);

    void getAllDofs(const FiniteElemSpace *feSpace, 
		    std::set<const DegreeOfFreedom*>& dofs);

    /** \brief
     * Creates a text file storing the value of a sparse matrix. Each line of the file
     * has three columns:
     *    row  col  value
     * This file can be used in Matlab using the command "spconvert".
     *
     * \param[in]  mat       The matrix which is used the write the text file.
     * \param[in]  filename  Name of the file to be created.
     */
    void writeMatlabMatrix(DOFMatrix &mat, std::string filename);

    /** \brief
     * Creates a text file storing the value of a sparse matrix. Each line of the file
     * has three columns:
     *    row  col  value
     * This file can be used in Matlab using the command "spconvert". The function 
     * works only for a matrix of DOFMatrices, that are all defined from the same
     * FE spaces.
     *
     * \param[in]  mat       The matrix which is used the write the text file.
     * \param[in]  filename  Name of the file to be created.     
     */
    void writeMatlabMatrix(Matrix<DOFMatrix*> &mat, std::string filename);

    /** \brief
     *
     */
    void writeMatlabVector(DOFVector<double> &vec, std::string filename);

    void writeMatlabVector(SystemVector &vec, std::string filename);

    void writeCoordsFile(const FiniteElemSpace *feSpace, std::string filename);

    void printElementHierarchie(Mesh *mesh, int elIndex);    

    void printElementRefinementSequence(Mesh *mesh, Element *el);

    int getLocalNeighbourIndex(Mesh *mesh, int elIndex, int neighIndex);

    void importDofVectorByCoords(DOFVector<double>* vec, std::string filename);

    void exportDofVectorByCoords(const DOFVector<double>* vec, 
				 std::string filename);

    void createNodeCoords(Mesh *mesh, ElementIdxToCoords& coords);
    
    void testNodeCoords(Mesh* mesh, ElementIdxToCoords& coords);
    
    /** \brief
     * Traverse a mesh and store for each element all its vertex DOFs in local sorted 
     * order (by values).
     *
     * \param[in]   mesh    Mesh to be traversed.
     * \param[out]  elMap   Stores to each element the vertex DOFs in sorted order.
     */
    void createSortedDofs(Mesh *mesh, ElementIdxToDofs &elMap);

    /** \brief
     * Takes a map from element indices to lists of DOFs. Checks, if for each element
     * in the mesh the vertex value order is still valid.
     *
     * The element index map must be created by the function \createSortedDofs. Using
     * both functions it can be checked if a renumbering of dofs does not changes the
     * local vertex value order (which is required by AMDiS to be always equal on each
     * element).
     *
     * If the test fails, the function prints some debug information to screen and
     * terminates the programm.
     *
     * \param[in]  mesh   Mesh to be traversed.
     * \param[in]  elMap  Map from element indices to lists of DOFs. It is used to check
     *                    the validaty as described above.
     */
    void testSortedDofs(Mesh *mesh, ElementIdxToDofs &elMap);

    /// Takes tree dofs and returns a list with the dofs sorted by their values.
    void sortDofs(const DegreeOfFreedom* dof0,
		  const DegreeOfFreedom* dof1,
		  const DegreeOfFreedom* dof2,
		  DofContainer &vec);

    /// Takes four dofs and returns a list with the dofs sorted by their values.
    void sortDofs(const DegreeOfFreedom* dof0,
		  const DegreeOfFreedom* dof1,
		  const DegreeOfFreedom* dof2,
		  const DegreeOfFreedom* dof3,
		  DofContainer &vec);    

    /** \brief
     * Takes to vectors of DOF indices and tests if they pairwise equal. To test
     * for equality, coordinates are checked. This makes it possible to test 
     * internal algorithms that manipulate the mesh having two DOFs at the same
     * geometrical position.
     *
     * This function does not return any value but abort the execution if it
     * find two DOFs not to be equal.
     *
     * \param[in]   feSpace   Finite element space which is used to create
     *                        coordinates.
     * \param[in]   dofs0     First DOF container.
     * \paran[in]   dofs1     Second DOF container.
     */
    void testDofsByCoords(const FiniteElemSpace *feSpace,
			  DofContainer &dofs0, 
			  DofContainer &dofs1);

    /** \brief
     * Works in the same way as described above, but the calling function
     * must provide a DOFVector which contains all DOF coordinates. This is
     * more efficient if the function is called multiple times.
     *
     * \param[in]   coords    DOFVector of DOF coordinates.
     * \param[in]   dofs0     First DOF container.
     * \paran[in]   dofs1     Second DOF container.
     */
    void testDofsByCoords(DOFVector<WorldVector<double> > &coords,
			  DofContainer &dofs0, 
			  DofContainer &dofs1); 
  }
}

#endif
