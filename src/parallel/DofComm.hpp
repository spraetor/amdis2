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



/** \file DofComm.h */

#ifndef AMDIS_DOF_COMM_H
#define AMDIS_DOF_COMM_H

#include <mpi.h>
#include <map>
#include "parallel/ParallelTypes.hpp"
#include "FiniteElemSpace.h"
#include "Global.h"

namespace AMDiS
{
  namespace Parallel
  {

    class DofComm
    {
    public:
      DofComm()
      {}

      typedef std::map<const FiniteElemSpace*, DofContainer> FeMapType;
      typedef FeMapType::iterator FeMapIter;
      typedef std::map<int, FeMapType> DataType;
      typedef DataType::iterator DataIter;

      void init(std::vector<const FiniteElemSpace*>& fe);

      void create(Mesh* mesh, InteriorBoundary& boundary);

      Mesh* getMesh()
      {
        return mesh;
      }

      DataType& getSendDofs()
      {
        return sendDofs;
      }

      DataType& getRecvDofs()
      {
        return recvDofs;
      }

      DataType& getPeriodicDofs()
      {
        return periodicDofs;
      }

      // Writes all data of this object to an output stream.
      void serialize(std::ostream& out);

      int getNumberDofs(DataType& data,
                        const FiniteElemSpace* feSpace,
                        bool countDouble = false);

      int getDegree(const FiniteElemSpace* feSpace,
                    const DegreeOfFreedom* dof);

    protected:
      void createContainer(Mesh* mesh,
                           InteriorBoundary& boundary,
                           RankToBoundMap& rankToBoundMap,
                           DataType& data);

    protected:
      /// This map contains for each rank the list of DOFs the current rank must
      /// send to exchange solution DOFs at the interior boundaries.
      DataType sendDofs;

      /// This map contains on each rank the list of DOFs from which the current
      /// rank will receive DOF values (i.e., this are all DOFs at an interior
      /// boundary). The DOF indices are given in rank's local numbering.
      DataType recvDofs;

      /// This map contains on each rank a list of DOFs along the interior bound-
      /// aries to communicate with other ranks. The DOF indices are given in rank's
      /// local numbering. Periodic boundaries within one subdomain are not
      /// considered here.
      DataType periodicDofs;

      std::vector<const FiniteElemSpace*> feSpaces;

      Mesh* mesh;

      friend class Iterator;

    public:
      class Iterator
      {
      public:
        Iterator(DataType& d,
                 const FiniteElemSpace* fe = NULL)
          : data(d),
            dofCounter(-1),
            traverseFeSpace(fe),
            removedDof(false)
        {
          goFirst();
        }

        inline bool end()
        {
          return (dataIter == data.end());
        }

        inline void nextRank()
        {
          do
          {
            ++dataIter;
          }
          while (setNextFeMap() == false);
        }

        inline void nextFeSpace()
        {
          ++feMapIter;
        }

        inline void next()
        {
          ++feMapIter;
          if (feMapIter == dataIter->second.end())
          {
            do
            {
              ++dataIter;
            }
            while (setNextFeMap() == false);
          }
          else
          {
            dofIter = feMapIter->second.begin();
            dofCounter = 0;
          }
        }

        inline void beginDofIter(const FiniteElemSpace* fe = NULL)
        {
          if (fe != NULL)
          {
            feMapIter = dataIter->second.begin();

            while (feMapIter->first != fe &&
                   feMapIter != dataIter->second.end())
              ++feMapIter;
          }

          if (feMapIter != dataIter->second.end())
          {
            dofIter = feMapIter->second.begin();
            dofCounter = 0;
          }
        }

        inline bool endDofIter()
        {
          if (feMapIter == dataIter->second.end())
            return true;

          return (dofIter == feMapIter->second.end());
        }

        inline void nextDof()
        {
          if (removedDof)
          {
            removedDof = false;
          }
          else
          {
            ++dofIter;
            ++dofCounter;
          }
        }

        inline void removeDof()
        {
          dofIter = feMapIter->second.erase(dofIter);
          removedDof = true;
        }

        inline int getRank()
        {
          return dataIter->first;
        }

        inline const FiniteElemSpace* getFeSpace()
        {
          return feMapIter->first;
        }

        inline DofContainer& getDofs()
        {
          return feMapIter->second;
        }

        inline const DegreeOfFreedom* getDof()
        {
          return *dofIter;
        }

        inline DegreeOfFreedom getDofIndex()
        {
          return **dofIter;
        }

        inline int getDofCounter()
        {
          return dofCounter;
        }

      protected:
        void goFirst()
        {
          dataIter = data.begin();

          while (setNextFeMap() == false)
            ++dataIter;
        }

        bool setNextFeMap();

      protected:
        DataType& data;

        DofComm::DataIter dataIter;

        DofComm::FeMapIter feMapIter;

        DofContainer::iterator dofIter;

        int dofCounter;

        const FiniteElemSpace* traverseFeSpace;

        bool removedDof;
      };

    };


    class MultiLevelDofComm
    {
    public:
      void init(MeshLevelData& levelData,
                std::vector<const FiniteElemSpace*>& fe);

      void create(Mesh* mesh, MultiLevelInteriorBoundary& boundary);


      inline DofComm& operator[](int level)
      {
        return levelDofComm[level];
      }

    private:
      std::map<int, DofComm> levelDofComm;
    };

  }
}

#endif // AMDIS_DOF_COMM_H
