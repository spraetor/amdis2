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



/** \file ParallelTypes.h */

#include <vector>
#include <set>
#include <map>
#include <boost/container/flat_map.hpp>
#include "BoundaryObject.h"
#include "Global.h"
#include "parallel/InteriorBoundary.h"

#ifndef AMDIS_PARALLEL_TYPES_H
#define AMDIS_PARALLEL_TYPES_H

namespace AMDiS
{
  namespace Parallel
  {

    /// Defines set of DOF indices.
    typedef std::set<DegreeOfFreedom> DofIndexSet;

    /// Defines a mapping type from DOFs to rank numbers.
    typedef std::map<const DegreeOfFreedom*, int> DofToRank;

    /// Defines a mapping type from DOFs to a set of rank numbers.
    typedef std::map<const DegreeOfFreedom*, std::set<int>> DofToPartitions;

    /// Defines a mapping type from DOF indices to a set of rank numbers.
    typedef std::map<DegreeOfFreedom, std::set<int>> DofIndexToPartitions;

    /// Defines a mapping type from rank numbers to sets of DOFs.
    typedef std::map<int, DofContainer> RankToDofContainer;

    /// Defines a mapping type from DOFs to boolean values.
    typedef std::map<const DegreeOfFreedom*, bool> DofToBool;

    /// Defines a mapping type from DOF indices to boolean values.
    typedef std::map<DegreeOfFreedom, bool> DofIndexToBool;

    typedef std::map<const DegreeOfFreedom*, DegreeOfFreedom> DofIndexMap;

    /// Maps a boundary type, i.e., a boundary identifier index, to a periodic
    /// DOF mapping.
    typedef std::map<BoundaryType, std::map<DegreeOfFreedom, DegreeOfFreedom>> PeriodicDofMap;

    /// Different FE spaces may have different DOFs on the same mesh. Thus we
    /// need to have a periodic DOF mapping for each FE space.
    typedef std::map<const FiniteElemSpace*, PeriodicDofMap> PeriodicDofMapFeSpace;

    /// Is used if a DOF index is mapped to multiple indices, i.e., to both, a local
    /// and a global one.
    struct MultiIndex
    {
      int local, global;
    };

    typedef boost::container::flat_map<DegreeOfFreedom, MultiIndex> DofMap;

    typedef std::vector<MeshStructure> MeshCodeVec;

    typedef std::map<int, std::vector<AtomicBoundary>> RankToBoundMap;

    typedef boost::container::flat_map<int, boost::container::flat_map<Mesh*, Element*>> MacroElIndexMap;

    typedef std::map<Mesh*, std::vector<const FiniteElemSpace*>> MeshToFeSpaces;
  }
} // end namespaces

#endif
