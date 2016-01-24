/** \file MacroElement.h */

#pragma once

#include <deque>
#include <stdio.h>

#include "AMDiS_fwd.h"
#include "Boundary.h"
#include "Global.h"
#include "Projection.h"
#include "FixVec.h"

namespace AMDiS
{
  /** \ingroup Triangulation
   * \brief
   * MacroElements form the macro triangulation of a Mesh. In a MacroElement
   * geometrical information are stored, which are used in mesh traversal,
   * to calculate the desired information and fill it in an ElInfo object.
   */
  class MacroElement
  {
  public:
    /// Creates a new MacroElement. The mesh is needed only to get the dimension
    MacroElement(int dim);
    
    /// Copy constructor
    MacroElement(MacroElement const& that);

    /// Destructor.
    virtual ~MacroElement();

    /// Copy assignment operator
    MacroElement& operator=(MacroElement const& that);

    /** \name getting methods
     * \{
     */

    /// Returns \ref index.
    int getIndex() const
    {
      return index;
    }

    /// Returns ref projection[i].
    Projection* getProjection(int i) const
    {
      return projection[i];
    }

    /// Returns \ref el
    Element* getElement() const
    {
      return element;
    }

    /// Returns the i-th neighbour of this MacroElement \ref neighbour[i]
    MacroElement* getNeighbour(int i) const
    {
      return neighbour[i];
    }

    /// Returns the i-th inverse neighbour of this MacroElement \ref neighbour_inv[i]
    /// Uses the inverse neighbour relation in graph-structured meshes where
    /// elements can have more than one neighbour
    /// If [b] is neighbour of [a], then is [a] neighbour_inv of [b]
    MacroElement* getNeighbourInv(int i) const
    {
      return neighbour_inv[i];
    }

    /// Returns the i-th opp-vertex of this MacroElement \ref oppVertex[i]
    int getOppVertex(int i) const
    {
      return oppVertex[i];
    }

    /// Returns \ref coord[i]
    WorldVector<double>& getCoord(int i)
    {
      return coord[i];
    }

    WorldVector<double> const& getCoord(int i) const
    {
      return coord[i];
    }

    /// Returns \ref coord
    FixVec<WorldVector<double>, VERTEX>& getCoord()
    {
      return coord;
    }

    FixVec<WorldVector<double>, VERTEX> const& getCoord() const
    {
      return coord;
    }

    /// Returns \ref boundary[i]
    BoundaryType getBoundary(int i) const
    {
      return boundary[i];
    }

    /// Returns \ref elType
    int getElType() const
    {
      return elType;
    }

    /** \} */

    /** \name setting methods
     * \{
     */

    /// Sets \ref index
    void setIndex(int n)
    {
      index = n ;
    }

    /// Sets \ref element if not yet set.
    void setElement(Element* el)
    {
      if (!element)
      {
        element = el;
      }
      else
      {
        if (element != el)
          ERROR("Trying to change element in MacroElement\n");
      }
    }

    /// Sets \ref elType
    void setElType(int typ)
    {
      elType = typ;
    }

    /// Sets \ref projection[i] = p.
    void setProjection(int i, Projection* p)
    {
      projection[i] = p;
    }

    /// Sets the i-th Neighbour to n
    void setNeighbour(int i, MacroElement* n)
    {
      neighbour[i] = n;
    }

    /// Sets the i-th opp vertex to c
    void setOppVertex(int i, int c)
    {
      oppVertex[i] = c;
    }

    /// Sets \ref boundary[i] to b
    void setBoundary(int i, BoundaryType b)
    {
      boundary[i] = b;
    }

    ///
    void setCoord(int i, const WorldVector<double> c)
    {
      coord[i] = c;
    }

    /** \} */

  protected:
    /// Element of this MacroElement.
    Element* element;

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


    friend class MacroInfo;
    friend class io::MacroReader;
  };

} // end namespace AMDiS
