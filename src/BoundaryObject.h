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



/** \file BoundaryObject.h */

#ifndef AMDIS_BOUNDARY_OBJECT_H
#define AMDIS_BOUNDARY_OBJECT_H

#include <vector>

#include "MacroElement.h"
#include "Element.h"
#include "Boundary.h"

namespace AMDiS {

  typedef std::vector<std::pair<GeoIndex, int> > ExcludeList;

  /// Defines the geometrical objects that forms the boundary;
  struct BoundaryObject {

    BoundaryObject();

    BoundaryObject(Element *e, 
		   int eType, 
		   GeoIndex sObj, 
		   int iObj, 
		   bool rMode = false);

    static bool computeReverseMode(BoundaryObject &obj0,
				   BoundaryObject &obj1,
				   const FiniteElemSpace *feSpace,
				   BoundaryType boundary);
    
    bool operator==(const BoundaryObject& other) const;
    
    bool operator!=(const BoundaryObject& other) const;

    bool operator<(const BoundaryObject& other) const;

    /// The macro element to which the boundary element corresponds to.
    Element* el;

    /// Index of the macro element.
    int elIndex;

    /// Element type index, only used in 3d.
    int elType;

    /// Defines the geometrical object at the boundary. It must be "a part" of
    /// the macro element \ref el, i.e., either 1 (a vertex), 2 (an edge) or 3 
    /// (a face).
    GeoIndex subObj;

    /** \brief
     * Defines which of vertex, edge or face of the macro element is part of the
     * boundary.
     *
     * Example: If the macro element is a triangle, than \ref subObj may be either
     * 1 (vertex) or 2 (edge). Assume its the last one. So this variable defines
     * which of the three possible edges of the triangle is at the interior
     * boundary.
     */
    int ithObj;

    bool reverseMode;

    /** \brief
     * In many situations it may be necessary to exclude some parts of the 
     * element to be part of the boundary. In 3d, when a face is part of the 
     * boundary, an edge or an vertex may be exludeded. In 2d only vertices may
     * be exluded to be part of an edge boundary. This list contains pairs of 
     * exludeded structures. The first component of every pair denotes if it is
     * a vertex or an edge, and the second component denotes the local index of
     * the structure.
     */
    ExcludeList excludedSubstructures;
  };



  /** \brief 
   * Defines one atomic part of the boundary, i.e., two boundary objects where
   * the boundary goes through.
   */
  struct AtomicBoundary {
    AtomicBoundary()
    : type(INTERIOR),
      maxLevel(0)
    {}

    bool operator==(const AtomicBoundary& other) const;

    bool operator!=(const AtomicBoundary& other) const;

    /// The rank's part of the boundary.
    BoundaryObject rankObj;

    /// The object on the other side of the boundary.
    BoundaryObject neighObj;

    /// Integer flag that is used to distinguish between different types of 
    /// boundaries. Till now it is used only for periodic boundaries, which are
    /// also handles as interior boundaries.
    BoundaryType type;

    int maxLevel;
  };


}

#endif
