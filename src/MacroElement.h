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



/** \file MacroElement.h */

#ifndef AMDIS_MACROELEMENT_H
#define AMDIS_MACROELEMENT_H

#include <deque>
#include <stdio.h>
#include "AMDiS_fwd.h"
#include "Boundary.h"
#include "Global.h"
#include "Projection.h"
#include "FixVec.h"
#include "Serializable.h"

namespace AMDiS {

  /** \ingroup Triangulation
   * \brief
   * MacroElements form the macro triangulation of a Mesh. In a MacroElement
   * geometrical information are stored, which are used in mesh traversal,
   * to calculate the desired information and fill it in an ElInfo object. 
   */ 
  class MacroElement : public Serializable
  {
  public:
    /// Creates a new MacroElement. The mesh is needed only to get the dimension
    MacroElement(int dim); 

    /// Destructor.
    virtual ~MacroElement();

    ///
    MacroElement& operator=(const MacroElement &el);

    /** \name getting methods
     * \{
     */

    /// Returns \ref index.
    inline int getIndex() const 
    {
      return index; 
    }

    /// Returns ref projection[i]. 
    inline Projection *getProjection(int i) const 
    {
      return projection[i];
    }

    /// Returns \ref el
    inline Element* getElement() const 
    {
      return element; 
    }

    /// Returns the i-th neighbour of this MacroElement \ref neighbour[i]
    inline MacroElement* getNeighbour(int i) const 
    {
      return neighbour[i];
    }

    /// Returns the i-th inverse neighbour of this MacroElement \ref neighbour_inv[i]
    /// Uses the inverse neighbour relation in graph-structured meshes where
    /// elements can have more than one neighbour
    /// If [b] is neighbour of [a], then is [a] neighbour_inv of [b]
    inline MacroElement* getNeighbourInv(int i) const 
    {
      return neighbour_inv[i];
    }

    /// Returns the i-th opp-vertex of this MacroElement \ref oppVertex[i]
    inline int getOppVertex(int i) const 
    {
      return oppVertex[i];
    }

    /// Returns \ref coord[i]
    inline WorldVector<double>& getCoord(int i) 
    {
      return coord[i];    
    }

    /// Returns \ref coord
    inline FixVec<WorldVector<double>, VERTEX>& getCoord() 
    {
      return coord;    
    }

    /// Returns \ref boundary[i]
    inline BoundaryType getBoundary(int i) const 
    {
      return boundary[i]; 
    }

    /// Returns \ref elType
    inline int getElType() const 
    {
      return elType; 
    }

    /** \} */

    /** \name setting methods
     * \{
     */

    /// Sets \ref index
    inline void setIndex(int n) 
    {
      index = n ; 
    }

    /// Sets \ref element if not yet set.
    inline void setElement(Element* el) 
    {
      if (!element) {
	element = el; 
      } else {
	if (element != el) 
	  ERROR("Trying to change element in MacroElement\n");   
      }
    }

    /// Sets \ref elType
    inline void setElType(int typ) 
    {
      elType = typ; 
    }

    /// Sets \ref projection[i] = p.
    inline void setProjection(int i, Projection *p) 
    {
      projection[i] = p;
    }

    /// Sets the i-th Neighbour to n
    inline void setNeighbour(int i, MacroElement *n) 
    {
      neighbour[i] = n;
    }

    /// Sets the i-th opp vertex to c
    inline void setOppVertex(int i, int c)
    {
      oppVertex[i] = c;
    }

    /// Sets \ref boundary[i] to b
    inline void setBoundary(int i, BoundaryType b) 
    {
      boundary[i] = b; 
    }

    ///
    inline void setCoord(int i, const WorldVector<double> c) 
    {
      coord[i] = c;
    }

    /** \} */


    /// Writes the macro element to a file.
    void serialize(std::ostream &out);

    /// Reads a macro element from a file.
    void deserialize(std::istream &in);

    ///
    inline void writeNeighboursTo(std::vector<int> *indices) 
    {
      deserializedNeighbourIndices = indices;
    }

  protected:
    /// Element of this MacroElement.
    Element *element;

    /// Coordinates of the vertices.
    FixVec<WorldVector<double>, VERTEX> coord;

    /// Boundary type of each boundary part (face/edge/vertex)
    FixVec<BoundaryType, BOUNDARY> boundary;

    /// Boundary projection to curved boundaries or element projections.
    FixVec<Projection*, PROJECTION> projection;

    /// Pointers to all neighbours of this MacroElement 
    FixVec<MacroElement*, NEIGH> neighbour;

    /// Pointers to all neighbours of this MacroElement 
    FixVec<MacroElement*, NEIGH> neighbour_inv;
    
    /// opp vertices of this MacroElement
    FixVec<int, NEIGH> oppVertex;

    /// index of this MacroElement
    int index;  

    /// Element type of the MacroElement
    int elType;  

    ///
    std::vector<int> *deserializedNeighbourIndices;


    friend class MacroInfo;
    friend class io::MacroReader;
    friend class ElInfo;
    friend class ElInfo1d;
    friend class ElInfo2d;
    friend class ElInfo3d;
  };

}


#endif  // AMDIS_MACROELEMENT_H 
