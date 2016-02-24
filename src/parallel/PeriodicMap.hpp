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



/** file PeriodicMap.h */

#ifndef AMDIS_PERIODIC_MAP
#define AMDIS_PERIODIC_MAP

#include <map>
#include <set>
#include <vector>
#include "AMDiS_fwd.h"
#include "parallel/ParallelTypes.hpp"
#include "Boundary.h"
#include "Serializer.h"

namespace AMDiS
{
  namespace Parallel
  {

    /** \brief
     * This class stores information about the periodic DOFs in the (sub)domain.
     * To each DOF on a periodic boundary there is the information to which DOF
     * is is periodic. Here we use global DOF indices. Furtheremore, a DOF can have
     * different periodic mapping. Assume we have a 2D box with all four edges
     * being periodic. Than, the four corner vertices have all two periodic
     * mapping. So, a periodic mapping is identified by the global DOF index to be
     * mapped and the boundary index.
     */
    class PeriodicMap
    {
    public:
      PeriodicMap() {}

      /// Reset all data.
      void clear()
      {
        periodicDofMap.clear();
        periodicDofAssociations.clear();
      }


      /// Get a periodic DOF mapping for a given FE space.
      inline PeriodicDofMap& getPeriodicMap(const FiniteElemSpace* feSpace)
      {
        return periodicDofMap[feSpace];
      }


      /** \brief
       * Map a DOF
       *
       * \param[in]  feSpace         FE space from which the DOF comes from.
       * \param[in]  type            Index of the periodic boundary. Is used to get
       *                             the correct mapping if the DOF has multiple
       *                             periodic associations.
       * \param[in]  globalDofIndex  Global DOF index.
       *
       * \return     Mapping of the global DOF index. The function fails if the
       *             the DOF is not periodic in the given FE space and periodic
       *             boundary type.
       */
      inline int map(const FiniteElemSpace* feSpace,
                     BoundaryType type,
                     int globalDofIndex)
      {
        FUNCNAME("PeriodicMap::map()");

        TEST_EXIT_DBG(periodicDofMap.count(feSpace))("Should not happen!\n");
        TEST_EXIT_DBG(periodicDofMap[feSpace][type].count(globalDofIndex) == 1)
        ("There is no periodic association for global DOF %d for boundary type %d!\n",
         globalDofIndex, type);

        return periodicDofMap[feSpace][type][globalDofIndex];
      }


      /// Adds a new periodic mapping. Fails if there is already a mapping for
      /// this DOFs that maps to a different DOF index than the given one.
      inline void add(const FiniteElemSpace* feSpace,
                      BoundaryType type,
                      DegreeOfFreedom dof0, DegreeOfFreedom dof1)
      {
        FUNCNAME("PeriodicMap::map()");

        TEST_EXIT_DBG(periodicDofMap[feSpace][type].count(dof0) == 0 ||
                      periodicDofMap[feSpace][type][dof0] == dof1)
        ("Should not happen!\n");

        periodicDofMap[feSpace][type][dof0] = dof1;
        periodicDofAssociations[feSpace][dof0].insert(type);
      }


      /// Adds a whole periodic mapping to the current one.
      void add(const FiniteElemSpace* feSpace, PeriodicDofMap& newMap);


      /// For a given global DOF index, this function returns the set of periodic
      /// associations, i.e., the boundary types the DOF is associated to, for
      /// this DOF.
      inline std::set<BoundaryType>& getAssociations(const FiniteElemSpace* feSpace,
          int globalDofIndex)
      {
        FUNCNAME("PeriodicMap::getAssociations()");

        TEST_EXIT_DBG(periodicDofAssociations.count(feSpace))
        ("Should not happen!\n");
        TEST_EXIT_DBG(periodicDofAssociations[feSpace].count(globalDofIndex))
        ("Should not happen!\n");

        return periodicDofAssociations[feSpace][globalDofIndex];
      }


      /** \brief
       * Given a global DOF index in some given finite element space, this
       * function creates all valid periodic associations of this DOF. If the
       * DOF is not periodic, this function does not fail but has no effect on
       * the output data.
       *
       * \param[in]  feSpace         feSpace on which the function should work on.
       * \param[in]  globalDofIndex  global index of a DOF.
       * \param[in]  elObjDb         Element object database that is used to check
       *                             if a given periodic index is valid.
       * \param[out] perAsc          set of periodic associations, is not deleted
       *                             by this function
       */
      void fillAssociations(const FiniteElemSpace* feSpace,
                            int globalDofIndex,
                            const ElementObjectDatabase& elObjDb,
                            std::set<int>& perAsc);


      /** \brief
       * Maps a given DOF index for all given periodic DOF associations.
       *
       * \param[in]   feSpace         feSpace on which the function should work on.
       * \param[in]   globalDofIndex  global index of a DOF.
       * \param[in]   perAsc          set of periodic associations.
       * \param[out]  mappedDofs      set of global DOF indices.
       */
      void mapDof(const FiniteElemSpace* feSpace,
                  int globalDofIndex,
                  const std::set<int>& perAsc,
                  std::vector<int>& mappedDofs);

      /** \brief
       * Maps a given DOF index pair for all given periodic DOF associations.
       *
       * \param[in]   rowFeSpace      feSpace of the DOFs on the first component.
       * \param[in]   colFeSpace      feSpace of the DOFs on the second component.
       * \param[in]   globalDofIndex  pair of global index of a DOF.
       * \param[in]   perAsc          set of periodic associations.
       * \param[out]  mappedDofs      set of pairs of global DOF indices.
       */
      void mapDof(const FiniteElemSpace* rowFeSpace,
                  const FiniteElemSpace* colFeSpace,
                  std::pair<int, int> globalDofIndex,
                  const std::set<int>& perAsc,
                  std::vector<std::pair<int, int>>& mappedDofs);


      /// Returns true, if the DOF (global index) is a periodic DOF.
      inline bool isPeriodic(const FiniteElemSpace* feSpace, int globalDofIndex)
      {
        return (periodicDofAssociations[feSpace].count(globalDofIndex) > 0 &&
                periodicDofAssociations[feSpace][globalDofIndex].size() > 0);
      }


      inline bool isPeriodicOnBound(const FiniteElemSpace* feSpace,
                                    BoundaryType type,
                                    int globalDofIndex)
      {
        return periodicDofAssociations[feSpace][globalDofIndex].count(type);
      }


      /// Returns true, if the DOF (global index) of a given FE space is a
      /// periodic DOF for the given boundary type.
      inline bool isPeriodic(const FiniteElemSpace* feSpace,
                             BoundaryType type,
                             int globalDofIndex)
      {
        return (periodicDofMap[feSpace][type].count(globalDofIndex) > 0);
      }

      /// Write the state of the object to serialization file.
      void serialize(std::ostream& out,
                     std::vector<const FiniteElemSpace*> feSpaces);

      /// Read the state of the object from serialization file.
      void deserialize(std::istream& in,
                       std::vector<const FiniteElemSpace*> feSpaces);

    private:
      /// Write \ref periodicDofMap to serialization file.
      void serialize(std::ostream& out, PeriodicDofMap& data);

      /// Write \ref periodicDofAssociations to serialization file
      void serialize(std::ostream& out, std::map<int, std::set<int>>& data);

      /// Read \ref periodicDofMap from serialization file.
      void deserialize(std::istream& in, PeriodicDofMap& data);

      /// Read \ref periodicDofAssociations from serialization file.
      void deserialize(std::istream& in, std::map<int, std::set<int>>& data);


    private:
      /// This map stores, for each FE space, a mapping from boundary indices
      /// to DOF mappings. So, if an FE space and a boundary index is provided,
      /// we eventually get a mapping from global DOF indices to global DOF
      /// indices.
      PeriodicDofMapFeSpace periodicDofMap;

      /// This map stores to each periodic DOF the set of periodic boundaries the
      /// DOF is associated to. In 2D, most DOFs are only on one periodic boundary.
      // Only, e.g., in a box with all boundaries being periodic, the four corners
      /// are associated by two different boundaries.
      std::map<const FiniteElemSpace*, std::map<DegreeOfFreedom, std::set<BoundaryType>>> periodicDofAssociations;

      friend class ParallelDebug;
    };
  }
}

#endif
