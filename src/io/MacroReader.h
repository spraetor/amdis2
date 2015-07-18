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



/** \file MacroReader.h */

#ifndef AMDIS_MACROREADER_H
#define AMDIS_MACROREADER_H

#include <deque>
#include "AMDiS_fwd.h"
#include "Global.h"
#include "FixVec.h"
#include "Boundary.h"

namespace AMDiS { namespace io {

  /** \defgroup Input Input module */

  /** \ingroup Input
   *  
   * \brief
   * Static class which reads a macro triangulation file and creates
   * the corresponding macro mesh.
   */
  class MacroReader
  {
    // Is used to read periodic macros
    class PeriodicMap
    {
    public:
      void setEntry(DegreeOfFreedom key, DegreeOfFreedom entry) 
      {
	// no trivial entries!
	if (key == entry) 
	  return;
	
	// if a key equal to entry exists ...
	if (getEntry(entry) >= 0) {
	  if (getEntry(entry) == key)
	    return;
	  
	  // ... let it be equal entries
	  setEntry(key, getEntry(entry));
	  return;
	}
	
	// replace entries equal to key
	for (std::map<DegreeOfFreedom, DegreeOfFreedom>::iterator it = 
	       periodicMap.begin(); it != periodicMap.end(); ++it)
	  if (it->second == key)
	    it->second = entry;
	
	// if key exists already, insert new entry with old entry as key
	if (getEntry(key) >= 0)
	  setEntry(getEntry(key), entry);
	
	// set entry
	periodicMap[key] = entry;
      }
            
      DegreeOfFreedom getEntry(DegreeOfFreedom key) 
      {
	std::map<DegreeOfFreedom, DegreeOfFreedom>::iterator it = periodicMap.find(key);
	if (it != periodicMap.end())
	  return it->second;
	return -1;
      }
    
    protected:
      std::map<DegreeOfFreedom, DegreeOfFreedom> periodicMap;
    };


  public:
    /// Creates a Mesh by reading the macro file with the given filename.
    static MacroInfo* readMacro(std::string filename, 
				Mesh* mesh,
				std::string periodicFile,
				int check);
    
    /** \brief
     * Used only for mesh repartition when the mesh is in initial state only 
     * containing macro elements without refinement even dofs. Then this
     * function will reset all the dofs of macros.
     */
    static void restoreMacroDofs(MacroInfo& macroInfo, int check = 1);
 
  protected:
    static void computeNeighbours(Mesh *mesh);

    static void boundaryDOFs(Mesh *mesh);

    static void umb(int *ele, Mesh *mesh,
		    void (*umbvk)(Mesh *mesh, MacroElement*,int k, int *el));

    static int macrotest(Mesh *mesh);

    static void macroTest(Mesh *mesh);

    static bool newEdge(Mesh *mesh, MacroElement *mel,
			int mel_edge_no, int *n_neigh);

    static void fillMelBoundary(Mesh *, MacroElement *mel,
				FixVec<BoundaryType,NEIGH> const&);

    static void fillMelNeigh(MacroElement *mel,
			     std::deque<MacroElement*>& macro_elements,
			     FixVec<int,NEIGH> const&);

    // for graph-structured meshes that have elements with >1 neighbors
    // (0)--[a]--(1)---[b]--(2)
    //            \____[c]__(3)
    //
    //  neigh:     [a]->[b], [b]->[c], [c]->[a]
    //  neigh_inv: [b]->[a], [c]->[b], [a]->[c]
    static void fillMelNeighInv(MacroElement *mel,
			     std::deque<MacroElement*>& macro_elements,
			     FixVec<int,NEIGH> const&);

    static void umbVkantMacro(Mesh *mesh,
			      MacroElement *,
			      int ka,
			      int *ele);

    static void recumb(Mesh *mesh, 
		       MacroElement *mel, MacroElement *macroalt,
		       int *test, double lg, int ka, int *ele, 
		       void (*umbvk)(Mesh *mesh, MacroElement*,int k, int *el));

    static void laengstekante(FixVec<WorldVector<double>,VERTEX> coord, 
			      double *l, int *v);

    static void checkMesh(Mesh *mesh);

    static int basicCheckFct(ElInfo* elInfo, Mesh *mesh);

    static void basicDOFCheckFct(ElInfo* elInfo, Mesh *mesh, int iadmin);

    static void basicNodeFct(ElInfo* elInfo, Mesh *mesh);

    friend class ::AMDiS::MacroInfo;
  };

} } // end namespace io, AMDiS

#endif
