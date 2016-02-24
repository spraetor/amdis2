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



/** \file InteriorBoundary.h */

#ifndef AMDIS_INTERIORBOUNDARY_H
#define AMDIS_INTERIORBOUNDARY_H

#include <vector>
#include <map>

#include "AMDiS_fwd.h"
#include "BoundaryObject.h"
#include "parallel/ParallelTypes.hpp"

namespace AMDiS
{
  namespace Parallel
  {

    /** \brief
     * Defines the interior boundary, i.e. a bound within the domain. It is used for
     * the classical domain decomposition parallelization.
     */
    class InteriorBoundary
    {
    public:
      void create(MeshLevelData& levelData,
                  int level,
                  ElementObjectDatabase& elObjDb);

      RankToBoundMap& getOwn()
      {
        return own;
      }

      inline MacroElIndexMap& getElIndexMap()
      {
        return macroElIndexMap;
      }

      inline Element* getElementPtr(int index, Mesh* mesh)
      {
        FUNCNAME_DBG("ElementObjectDatabase::getElementPtr()");
        TEST_EXIT_DBG(macroElIndexMap[index][mesh])
        ("No element pointer in macroElIndex map. Something is wrong.\n");
        return macroElIndexMap[index][mesh];
      }

      RankToBoundMap& getOther()
      {
        return other;
      }

      RankToBoundMap& getPeriodic()
      {
        return periodic;
      }

      std::vector<AtomicBoundary>& getOwn(int rank)
      {
        return own[rank];
      }

      std::vector<AtomicBoundary>& getOther(int rank)
      {
        return other[rank];
      }

      std::vector<AtomicBoundary>& getPeriodic(int rank)
      {
        return periodic[rank];
      }

      bool hasPeriodic()
      {
        return static_cast<bool>(periodic.size());
      }

      int getDegreeOwn(BoundaryObject& bObj);

      /// Writes this object to a file.
      void serialize(std::ostream& out);

      /// Reads the state of an interior boundary from a file.
      void deserialize(std::istream& in, Mesh* mesh);

    private:
      /// In this function, we put some verification procedures to check for
      /// consistency if the created interior boundary data. This function is
      /// enabled even for optimized mode to check all possible meshes. Thus,
      /// everything herein should be fast. If the code is stable, we can think
      /// to remove this function.
      void verifyBoundary();

      AtomicBoundary& getNewOwn(int rank);

      AtomicBoundary& getNewOther(int rank);

      AtomicBoundary& getNewPeriodic(int rank);

      /// Checks whether the given boundary is already a other boundary with
      /// given rank number. Returns true, if the bound is already in the other
      /// boundary database.
      bool checkOther(AtomicBoundary& bound, int rank);

      /// Removes the given boundary object from all owned boundaries. Returns
      /// true, if at least one object has been removed, otherwise false.
      bool removeOwn(BoundaryObject& bound);

      /// Removes the given boundary object from all other boundaries. Returns
      /// true, if at least one object has been removed, otherwise false.
      bool removeOther(BoundaryObject& bound);

      void serialize(std::ostream& out, RankToBoundMap& boundary);

      void deserialize(std::istream& in,
                       RankToBoundMap& boundary,
                       std::map<int, Element*>& elIndexMap);

      void serializeExcludeList(std::ostream& out, ExcludeList& list);

      void deserializeExcludeList(std::istream& in, ExcludeList& list);

    private:
      RankToBoundMap own, other, periodic;

      MacroElIndexMap macroElIndexMap;

      friend class ParallelDebug;

    public:
      /// Iterator for the interior boundary object.
      class iterator
      {
      public:
        iterator(RankToBoundMap& b)
          : bound(b),
            level(0)
        {
          reset();
        }

        /// Set the iterator to the first position.
        void reset()
        {
          mapIt = bound.begin();
          nextNonempty();
        }

        /// Test if iterator is at the final position.
        bool end() const
        {
          return (mapIt == bound.end());
        }

        /// Move iterator to the next position.
        void operator++()
        {
#if 0
          do
          {
            ++vecIt;
          }
          while (vecIt != mapIt->second.end() && vecIt->maxLevel < level);
#else
          ++vecIt;
#endif

          if (vecIt == mapIt->second.end())
          {
            ++mapIt;
            nextNonempty();
          }
        }

        inline AtomicBoundary& operator*()
        {
          return *vecIt;
        }

        inline AtomicBoundary* operator->()
        {
          return &(*vecIt);
        }

        void nextRank()
        {
          ++mapIt;
          nextNonempty();
        }

        inline int getRank()
        {
          return mapIt->first;
        }

      protected:

        inline void nextNonempty()
        {
          do
          {
            // Return, we are at the end.
            if (mapIt == bound.end())
              return;

            // Search for the next non empty boundary map.
            while (mapIt->second.size() == 0)
            {
              ++mapIt;
              if (mapIt == bound.end())
                return;
            }

            vecIt = mapIt->second.begin();

            // Search for the next atomic boundary on the mesh level
#if 0
            while (vecIt != mapIt->second.end() && vecIt->maxLevel < level)
              ++vecIt;
#endif

            // If vector iterator is not at the end, we have found one and
            // can return.
            if (vecIt != mapIt->second.end())
              return;

            // In this case, no boundary on the given level is found, continue
            // with next rank.
            ++mapIt;
          }
          while (true);
        }

      protected:
        RankToBoundMap::iterator mapIt;

        std::vector<AtomicBoundary>::iterator vecIt;

        RankToBoundMap& bound;

        int level;
      };
    };


    class MultiLevelInteriorBoundary
    {
    public:
      void create(MeshLevelData& levelData,
                  ElementObjectDatabase& elObjDb);

      inline InteriorBoundary& operator[](int level)
      {
        TEST_EXIT_DBG(levelIntBound.count(level))("Should not happen!\n");

        return levelIntBound[level];
      }

      /// Writes this object to a file.
      void serialize(std::ostream& out);

      /// Reads the state of an interior boundary from a file.
      void deserialize(std::istream& in, Mesh* mesh);

    private:
      std::map<int, InteriorBoundary> levelIntBound;
    };

  }
}

#endif // AMDIS_INTERIORBOUNDARY_H
