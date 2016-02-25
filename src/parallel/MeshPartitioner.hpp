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



/** \file MeshPartitioner.h */

#ifndef AMDIS_MESH_PARTITIONER_H
#define AMDIS_MESH_PARTITIONER_H

#include <map>
#include <set>
#include <mpi.h>

#include "AMDiS_fwd.h"
#include "Mesh.h"
#include "parallel/MpiHelper.hpp"
#include "parallel/ParallelDofMapping.hpp"


namespace AMDiS
{
  namespace Parallel
  {

    enum PartitionMode
    {
      INITIAL = 0,          // initial partitioning of a unpartitioned mesh
      ADAPTIVE_REPART = 1,  // repartitioning of a adaptively refined mesh
      REFINE_PART = 2       // quality improvement of the current partitioning
    };


    /**
     * Abstract class for mesh partitioning. This class provides only a function
     * for a random initial partitioning. A concrete partition must override the
     * functions \ref MeshPartitioner::partition and
     * \ref MeshPartitioner::createPartitionMap.
     */
    class MeshPartitioner
    {
    public:
      MeshPartitioner(std::string name, MPI::Intracomm* comm)
        : initFileStr(name),
          mpiComm(comm),
          mesh(NULL),
          boxPartitioning(false),
          mapLocalGlobal(NULL)
      {}

      virtual ~MeshPartitioner() {}

      /** \brief
       * Creates an initial partitioning of the AMDiS mesh. This partitioning
       * can be arbitrary, the only requirement is that each macro element
       * must be uniquely assign to a rank.
       *
       * \return   If the function returns true, the initial partitioning should
       *           be used for the very first computations. This can be e.g. due
       *           to the usage of an initial partitioning file.
       */
      virtual bool createInitialPartitioning();

      /** \brief
       * Function the takes a weighted set of macro elements and returns a
       * macro mesh partitioning. This function is virtual and must be implemented
       * for a specific algorithm or an external partitioning library.
       *
       * \param[in]  elemWeights   Maps to each macro element in rank's subdomain
       *                           a weight, which is usually the number of leaf
       *                           elements in this macro element.
       * \param[in]  mode          Most external partitioning libraries can make
       *                           a difference whether we want to create a
       *                           first partitioning or we already have created
       *                           one using this library but due to some mesh
       *                           adaptivity we want to repartition the mesh. In
       *                           the later case, the libraries also consider the
       *                           time for redistribution of the new partitioning.
       * \return     Returns a boolean value if the partitioning algorithm created
       *             a correct partitioning. If it is so, the partitioning is
       *             stored in \ref elementInRank and \ref partitionMap.
       */
      virtual bool partition(std::map<int, double>& elemWeights,
                             PartitionMode mode = INITIAL) = 0;

      virtual void createPartitionMap(std::map<int, int>& partitionMap) = 0;

      /// Write partitioner state to disk.
      void serialize(std::ostream& out);

      /// Read partitioner state from disk.
      void deserialize(std::istream& in);

      void setMesh(Mesh* m)
      {
        mesh = m;
      }

      void setBoxPartitioning(bool b)
      {
        boxPartitioning = b;
      }

      void setLocalGlobalDofMap(DofMap* m)
      {
        mapLocalGlobal = m;
      }

      Mesh* getMesh()
      {
        return mesh;
      }

      std::map<int, bool>& getElementInRank()
      {
        return elementInRank;
      }

      std::map<int, std::vector<int>>& getRecvElements()
      {
        return recvElements;
      }

      std::map<int, std::vector<int>>& getSendElements()
      {
        return sendElements;
      }

      /// After mesh repartition this function returns true if the mesh must be
      /// redistributed on at least one rank.
      bool meshChanged()
      {
        int nChanges = recvElements.size() + sendElements.size();
        mpi::globalAdd(nChanges);

        return static_cast<bool>(nChanges);
      }

    protected:
      /// Prefix for reading parameters from init file.
      std::string initFileStr;

      /// Pointer to the MPI communicator the mesh partitioner should make use of.
      MPI::Intracomm* mpiComm;

      /// Pointer to the AMDiS mesh.
      Mesh* mesh;

      /// The mesh partitioner can be used in to different modes, the standard
      /// mode and the so called "box partitioning". The standard mode assigns
      /// macro elements to ranks. If box partitioning is enabled, which makes
      /// only sense if the macro mesh results from meshconv's "lego mesher",
      /// then in 2D boxed (2 macro elements) and in 3D cubes (6 macro
      /// elements) are assigned as a union to ranks.
      bool boxPartitioning;

      /// In box partitioning mode this map stores for each box number the set
      /// of macro element indices the box consists of.
      std::map<int, std::set<int>> boxSplitting;

      /// In box partitioning mode this map stores to each box number the set
      /// of neighbouring boxes.
      std::map<int, std::set<int>> boxNeighbours;

      /// Is the reverse of the map \ref boxSplitting. Thus, it stores for each
      /// macro element index the box number it belongs to.
      std::map<int, int> elInBox;

      DofMap* mapLocalGlobal;

      std::map<int, std::vector<int>> elNeighbours;

      /// Maps to each macro element index (or box index in box
      /// partitioning mode) if it is in rank's partition or not.
      std::map<int, bool> elementInRank;

      /// Maps to each macro element index (or box index in box
      /// partitioning mode) the rank number the element belongs to.
      std::map<int, int> partitionMap;

      /// After mesh repartitioning these maps stores which elements are communicated
      /// from this rank to other ranks.
      std::map<int, std::vector<int>> recvElements, sendElements;
    };
  }
}

#endif
