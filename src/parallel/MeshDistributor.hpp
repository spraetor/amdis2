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



/** \file MeshDistributor.h */

#ifndef AMDIS_MESHDISTRIBUTOR_H
#define AMDIS_MESHDISTRIBUTOR_H


#include <mpi.h>
#include "parallel/DofComm.hpp"
#include "parallel/ElementObjectDatabase.hpp"
#include "parallel/ParallelTypes.hpp"
#include "parallel/MeshLevelData.hpp"
#include "parallel/MeshPartitioner.hpp"
#include "parallel/InteriorBoundary.hpp"
#include "parallel/ParallelDofMapping.hpp"
#include "parallel/PeriodicMap.hpp"
#include "parallel/StdMpi.hpp"
#include "AMDiS_fwd.h"
#include "Containers.h"
#include "Global.h"
#include "ProblemTimeInterface.h"
#include "ProblemIterationInterface.h"
#include "FiniteElemSpace.h"
#include "Serializer.h"
#include "BoundaryManager.h"
#include <string>

#include "operations/functors.hpp"

namespace AMDiS
{
  namespace Parallel
  {


    struct BoundaryDofInfo
    {
      std::map<GeoIndex, DofContainerSet> geoDofs;
    };


    class MeshDistributor
    {
    private:
      MeshDistributor();

    public:
      ~MeshDistributor();

      /// Initialization of mesh distributor.
      void initParallelization();

      /// Clean up procedure for the mesh distributor and attached objects.
      void exitParallelization();

      /** \brief
       * Register a parallel DOF mapping. This DOF mapping object will than
       * automatically updated by the mesh distributer after mesh changes.
       *
       * \param[in]  dofMap   Parallel DOF mapping object.
       */
      void registerDofMap(ParallelDofMapping& dofMap);

      /** \brief
       * Removes a registered DOF mapping from the mesh distributor.
       *
       * \param[in] dofMap   Parallel DOF mapping object to be removed.
       */
      void removeDofMap(ParallelDofMapping& dofMap);

      /// Adds a DOFVector to the set of \ref interchangeVecs. Thus, this vector
      /// will be automatically interchanged between ranks when mesh is
      /// repartitioned.
      template<typename T>
      void addInterchangeVector(DOFVector<T>* vec) {}
      void addInterchangeVector(DOFVector<double>* vec)
      {
        interchangeVectors.push_back(vec);
      }

      /// Removes the pointer to DOFVector @param vec from the
      /// set of interchange vectors.
      template<typename T>
      void removeInterchangeVector(DOFVector<T>* vec) {}
      void removeInterchangeVector(DOFVector<double>* vec)
      {
        std::vector<DOFVector<double>*>::iterator it;
        it = std::find(interchangeVectors.begin(), interchangeVectors.end(), vec);
        if ( it != interchangeVectors.end())
          interchangeVectors.erase(it);
      }

      /// Adds all DOFVectors of a SystemVector to \ref interchangeVecs.
      void addInterchangeVector(SystemVector* vec);

      /// The same as for DOFVectors
      void removeInterchangeVector(SystemVector* vec);

      /** \brief
       * This function checks if the mesh has changed on at least one rank. In
       * this case, the interior boundaries are adapted on all ranks such that
       * they fit together on all ranks. Furthermore the function
       * \ref updateLocalGlobalNumbering() is called to update the DOF numberings
       * and mappings on all rank due to the new mesh structure.
       *
       * \param[in]  tryRepartition   If this parameter is true, repartitioning
       *                              may be done. This depends on several other
       *                              parameters. If the parameter is false, the
       *                              mesh is only checked and adapted but never
       *                              repartitioned.
       */
      void checkMeshChange(bool tryRepartition = true);

      /// Checks if is required to repartition the mesh. If this is the case, a new
      /// partition will be created and the mesh will be redistributed between the
      /// ranks.
      bool repartitionMesh();


      void getImbalanceFactor(double& imbalance,
                              int& minDofs,
                              int& maxDofs,
                              int& sumDofs);

      double getImbalanceFactor();

      /// Calculates the imbalancing factor and prints it to screen.
      void printImbalanceFactor();

      /// Test, if the mesh consists of macro elements only. The mesh partitioning
      /// of the parallelization works for macro meshes only and would fail, if the
      /// mesh is already refined in some way. Therefore, this function will exit
      /// the program if it finds a non macro element in the mesh.
      void testForMacroMesh();

      inline std::string getName()
      {
        return name;
      }

      inline Mesh* getMacroMesh()
      {
        return macroMesh;
      }

      inline Mesh* getMesh(int i = 0)
      {
        return meshes[i];
      }

      inline int getNumberOfMeshes()
      {
        return meshes.size();
      }

      /// Returns the periodic mapping handler, \ref periodicMap.
      inline PeriodicMap& getPeriodicMap()
      {
        return periodicMap;
      }

      //     DofComm& getDofComm(int level)
      //     {
      //       return dofComm[level];
      //     }

      DofComm& getDofComm(Mesh* mesh, int level)
      {
        return dofComms[mesh][level];
      }

      std::map<Mesh*, MultiLevelDofComm>& getDofComms()
      {
        return dofComms;
      }

      InteriorBoundary& getIntBoundary(int level)
      {
        return intBoundary[level];
      }

      std::map<int, int>& getPartitionMap()
      {
        return partitionMap;
      }

      inline long getLastMeshChangeIndex()
      {
        int overallMeshChangeIndex = 0;
        for(size_t i = 0; i < meshes.size(); i++)
        {
          overallMeshChangeIndex += lastMeshChangeIndexs[meshes[i]];
        }
        return overallMeshChangeIndex;
      }

      inline long getLastMeshChangeIndex(Mesh* m)
      {
        return lastMeshChangeIndexs[m];
      }

      inline int getMpiRank()
      {
        return mpiRank;
      }

      inline int getMpiSize(int level)
      {
        return levelData.getMpiComm(level).Get_size();
      }

      inline MPI::Intracomm& getMpiComm(int level)
      {
        return levelData.getMpiComm(level);
      }

      inline bool isInitialized()
      {
        return initialized;
      }

      // Writes all data of this object to an output stream.
      void serialize(std::ostream& out);

      // Reads the object data from an input stream.
      void deserialize(std::istream& in);

      /// Works quite similar to the function \ref synchVector, but instead the
      /// values of subdomain vectors are combined along the boundaries, by a
      /// binary functor.
      // minorRank => majorRank
      template<typename T, typename Operator>
      void synchVector(DOFVector<T>& vec, Operator op)
      {
        const FiniteElemSpace* fe = vec.getFeSpace();
        MultiLevelDofComm& dofComm = dofComms[fe->getMesh()];

        int nLevels = levelData.getNumberOfLevels();
        for (int level = nLevels - 1; level >= 0; level--)
        {
          StdMpi<std::vector<T>> stdMpi(levelData.getMpiComm(level));

          for (DofComm::Iterator it(dofComm[level].getRecvDofs(), fe);
               !it.end(); it.nextRank())
          {
            std::vector<T> dofs;
            dofs.reserve(it.getDofs().size());

            for (; !it.endDofIter(); it.nextDof())
              dofs.push_back(vec[it.getDofIndex()]);

            stdMpi.send(it.getRank(), dofs);
          }

          for (DofComm::Iterator it(dofComm[level].getSendDofs());
               !it.end(); it.nextRank())
            stdMpi.recv(it.getRank());

          stdMpi.startCommunication();

          for (DofComm::Iterator it(dofComm[level].getSendDofs(), fe);
               !it.end(); it.nextRank())
            for (; !it.endDofIter(); it.nextDof())
              op(vec[it.getDofIndex()],
                 stdMpi.getRecvData(it.getRank())[it.getDofCounter()]);
        }
        synchVector(vec);
      }

      /** \brief
       * This function must be used if the values of a DOFVector must be
       * synchronised over all ranks. That means, that each rank sends the
       * values of the DOFs, which are owned by the rank and lie on an interior
       * boundary, to all other ranks also having these DOFs.
       *
       * This function must be used, for example, after the linear system is
       * solved, or after the DOFVector is set by some user defined functions,
       * e.g., initial solution functions.
       */
      // majorRank => minorRank
      template<typename T>
      void synchVector(DOFVector<T>& vec)
      {
        const FiniteElemSpace* fe = vec.getFeSpace();
        MultiLevelDofComm& dofComm = dofComms[fe->getMesh()];

        int nLevels = levelData.getNumberOfLevels();
        for (int level = nLevels - 1; level >= 0; level--)
        {
          StdMpi<std::vector<T>> stdMpi(levelData.getMpiComm(level));

          for (DofComm::Iterator it(dofComm[level].getSendDofs(), fe);
               !it.end(); it.nextRank())
          {

            std::vector<T> dofs;
            dofs.reserve(it.getDofs().size());

            for (; !it.endDofIter(); it.nextDof())
              dofs.push_back(vec[it.getDofIndex()]);

            stdMpi.send(it.getRank(), dofs);
          }

          for (DofComm::Iterator it(dofComm[level].getRecvDofs());
               !it.end(); it.nextRank())
            stdMpi.recv(it.getRank());

          stdMpi.startCommunication();

          for (DofComm::Iterator it(dofComm[level].getRecvDofs(), fe);
               !it.end(); it.nextRank())
            for (; !it.endDofIter(); it.nextDof())
              vec[it.getDofIndex()] =
                stdMpi.getRecvData(it.getRank())[it.getDofCounter()];
        }
      }

      /// Works in the same way as the function above defined for DOFVectors. Due
      /// to performance, this function does not call \ref synchVector for each
      /// DOFVector, but instead sends all values of all DOFVectors all at once.
      void synchVector(SystemVector& vec);

      /// Works quite similar to the function \ref synchVector, but instead the
      /// values of subdomain vectors are add along the boundaries.
      // minorRank => majorRank
      template<typename T>
      void synchAddVector(DOFVector<T>& vec)
      {
        const FiniteElemSpace* fe = vec.getFeSpace();
        MultiLevelDofComm& dofComm = dofComms[fe->getMesh()];

        int nLevels = levelData.getNumberOfLevels();
        for (int level = nLevels - 1; level >= 0; level--)
        {
          StdMpi<std::vector<T>> stdMpi(levelData.getMpiComm(level));

          for (DofComm::Iterator it(dofComm[level].getRecvDofs(), fe);
               !it.end(); it.nextRank())
          {
            std::vector<T> dofs;
            dofs.reserve(it.getDofs().size());

            for (; !it.endDofIter(); it.nextDof())
              dofs.push_back(vec[it.getDofIndex()]);

            stdMpi.send(it.getRank(), dofs);
          }

          for (DofComm::Iterator it(dofComm[level].getSendDofs());
               !it.end(); it.nextRank())
            stdMpi.recv(it.getRank());

          stdMpi.startCommunication();

          for (DofComm::Iterator it(dofComm[level].getSendDofs(), fe);
               !it.end(); it.nextRank())
            for (; !it.endDofIter(); it.nextDof())
              vec[it.getDofIndex()] +=
                stdMpi.getRecvData(it.getRank())[it.getDofCounter()];
        }

        synchVector(vec);
      }

      /// In 3D, a subdomain may not be a valid AMDiS mesh if it contains two
      /// parts which are only connected by an edge. In this case, the standard
      /// refinement algorithm does not work correctly, as two elements connected
      /// only on one edge are not neighours by definition. This functions checks
      /// for this situation and fix the problem. For this, the mesh is search for
      /// all edges connecting two elements that are otherwise not connected.
      void fix3dMeshRefinement();

      /** \brief Is used only within \ref fix3dMeshRefinement.
       *
       * \param[in]  elems            Set of macro element indices.
       * \param[out] disconnectedEls  On output, this vector contains sets of
       *                              element indices. The union is equal to elems.
       *                              Each set contains all element indices, which
       *                              are reachable among each other by neighbour
       *                              relations. Elements within two different sets
       *                              cannot be reached via neigbourhood relation.
       */
      void helpToFix(std::set<int>& elems,
                     std::vector<std::set<int>>& disconnectedEls);

      void setBoundaryDofRequirement(Flag flag)
      {
        createBoundaryDofFlag |= flag;
      }

      BoundaryDofInfo& getBoundaryDofInfo(const FiniteElemSpace* feSpace,
                                          int level)
      {
        FUNCNAME("MeshDistributor::getBoundaryDofInfo()");

        TEST_EXIT_DBG(level < static_cast<int>(boundaryDofInfo.size()))
        ("Wrong level number: %d, whereas array size is %d!\n",
         level, boundaryDofInfo.size());

        return boundaryDofInfo[level][feSpace];
      }

      void getAllBoundaryDofs(const FiniteElemSpace* feSpace,
                              int level,
                              DofContainer& dofs);

      ElementObjectDatabase& getElementObjectDb()
      {
        return elObjDb;
      }

      /// Adds a stationary problem to the global mesh distributor objects.
      static void addProblemStatGlobal(ProblemStatSeq* probStat);

      MeshLevelData& getMeshLevelData()
      {
        return levelData;
      }

      /// Update dof communicators, boundary dof info and the parallel dof mappings.
      /// If it is called for all meshes, \ref updateLocalGlobalNumbering is automatically
      /// called inside. If it is used for each mesh seperately, please don't forget to
      /// add \ref updateLocalGlobalNumbering to update the global matrix index.
      void updateDofRelatedStruct();

      void updateDofRelatedStruct(Mesh* mesh);

      void updateLocalGlobalNumbering();

      /// set variable \ref repartitioningAllowed
      void setRepartitioningAllowed(bool allowed)
      {
        repartitioningAllowed = allowed;
      }

      void setElementWeights(std::map<int, double>& elWgts)
      {
        elemWeights = elWgts;
      }

    protected:
      /// Rebuild only part of the mesh domain, which is necessary
      void quickRepartition(Mesh* mesh);

      /// Rebuild whole mesh domain
      void fullRepartition(Mesh* mesh);

      /// Updates all registered parallel DOF mappings, see \ref dofMaps.
      void updateDofsToDofMapping(Mesh* mesh = NULL);

      /// Updates the DOF after the mesh has been changed, see \ref dofMaps.
      void updateDofsToDofMapping(ParallelDofMapping& dmap,
                                  const FiniteElemSpace* feSpace);

      /// Checks if repartition is needed.
      bool isRepartitionNecessary();

      /// Creates an initial partitioning of the mesh.
      void createInitialPartitioning();

      /// Set for each element on the partitioning level the number of
      /// leaf elements.
      void setInitialElementWeights();

      /// Calculates \ref elemWeights with the gloabl max weight and
      /// global sum of weight.
      void calculateElemWeights();

      ///
      void addProblemStat(ProblemStatSeq* probStat);

      /// Determines the interior boundaries, i.e. boundaries between ranks, and
      /// stores all information about them in \ref interiorBoundary.
      void createInteriorBoundary(bool firstCall);

      ///
      void createBoundaryDofs(Mesh* mesh = NULL);

      /// Removes all macro elements from the mesh that are not part of ranks
      /// partition.
      void removeMacroElements();

      /// Calls \ref createPeriodicMap(feSpace) for all FE spaces that are
      /// handled by the mesh distributor.
      void createPeriodicMap();

      /// Creates, for a specific FE space, to all DOFs in rank's partition that
      /// are on a periodic boundary the mapping from dof index to the other
      /// periodic dof indices. This information is stored in \ref periodicDofMap.
      void createPeriodicMap(const FiniteElemSpace* feSpace);

      /// This function is called only once during the initialization when the
      /// whole macro mesh is available on all cores. It copies the pointers of all
      /// macro elements to \ref allMacroElements and stores all neighbour
      /// information based on macro element indices (and not pointer based) in
      /// \ref macroElementNeighbours. These information are then used to
      /// reconstruct macro elements during mesh redistribution.
      void createMacroElementInfo();

      void updateMacroElementInfo();

      /** \brief
       * Checks for all given interior boundaries if the elements fit together on
       * both sides of the boundaries. If this is not the case, the mesh is
       * adapted. Because refinement of a certain element may forces the
       * refinement of other elements, it is not guaranteed that all rank's meshes
       * fit together after this function terminates. Hence, it must be called
       * until a stable mesh refinement is reached.
       *
       * \param[in] allBound   Defines a map from rank to interior boundaries
       *                       which should be checked.
       * \param[in] mesh       The mesh the interior boundaries belong to.
       *
       * \return    If the mesh has  been changed by this function, it returns
       *            true. Otherwise, it returns false, i.e., the given interior
       *            boundaries fit together on both sides.
       */
      bool checkAndAdaptBoundary(RankToBoundMap& allBound, Mesh* mesh);

      /// Removes all periodic boundary condition information from all matrices and
      /// vectors of all stationary problems and from the mesh itself.
      void removePeriodicBoundaryConditions();

      /// Removes all periodic boundary condition information from all matrices and
      /// vector of a given stationary problem.
      void removePeriodicBoundaryConditions(ProblemStatSeq* probStat);

      // Removes all periodic boundaries from a given boundary map.
      void removePeriodicBoundaryConditions(BoundaryIndexMap& boundaryMap);

      void createMeshLevelStructure();

      /// Writes a vector of dof pointers to an output stream.
      void serialize(std::ostream& out, DofContainer& data);

      /// Writes a \ref RankToDofContainer to an output stream.
      void serialize(std::ostream& out,
                     std::map<int, std::map<const FiniteElemSpace*, DofContainer>>& data);

      /// Reads a vector of dof pointers from an input stream.
      void deserialize(std::istream& in, DofContainer& data,
                       std::map<int, const DegreeOfFreedom*>& dofIndexMap);

      /// Reads a \ref RankToDofContainer from an input stream.
      void deserialize(std::istream& in,
                       std::map<int, std::map<const FiniteElemSpace*, DofContainer>>& data,
                       std::map<const FiniteElemSpace*, std::map<int, const DegreeOfFreedom*>>& dofIndexMap);

      /// Writes a mapping from dof pointers to some values to an output stream.
      template<typename T>
      void serialize(std::ostream& out, std::map<const DegreeOfFreedom*, T>& data)
      {
        FUNCNAME("ParallelDomainBase::serialize()");

        int mapSize = data.size();
        SerUtil::serialize(out, mapSize);
        for (typename std::map<const DegreeOfFreedom*, T>::iterator it = data.begin();
             it != data.end(); ++it)
        {
          int v1 = (*(it->first));
          T v2 = it->second;
          SerUtil::serialize(out, v1);
          SerUtil::serialize(out, v2);
        }
      }

      /// Reads a mapping from dof pointer to some values from an input stream.
      template<typename T>
      void deserialize(std::istream& in, std::map<const DegreeOfFreedom*, T>& data,
                       std::map<int, const DegreeOfFreedom*>& dofIndexMap)
      {
        FUNCNAME("ParallelDomainBase::deserialize()");

        int mapSize = 0;
        SerUtil::deserialize(in, mapSize);
        for (int i = 0; i < mapSize; i++)
        {
          int v1 = 0;
          T v2;
          SerUtil::deserialize(in, v1);
          SerUtil::deserialize(in, v2);

          TEST_EXIT_DBG(dofIndexMap.count(v1) != 0)
          ("Cannot find DOF %d in map!\n", v1);

          data[dofIndexMap[v1]] = v2;
        }
      }

    protected:
      /// List of all stationary problems that are managed by this mesh
      /// distributor.
      std::vector<ProblemStatSeq*> problemStat;

      /// If true, the mesh distributor is already initialized;
      bool initialized;

      /// The rank of the current process.
      int mpiRank;

      /// Name of the problem (as used in the init files)
      std::string name;

      /// Set of all different FE spaces.
      std::vector<const FiniteElemSpace*> feSpaces;

      /// Always equal to meshes[0] which is used as macro
      /// mesh. For example, passed to \ref meshPartitioner.
      Mesh* macroMesh;

      /// Meshes to be managed for parallelization. Currently only two meshes
      /// are allowed since multi mesh method is limited to two meshes.
      std::vector<Mesh*> meshes;

      /// Stores the map of meshes and the corresponding FE spaces defined on them
      MeshToFeSpaces meshToFeSpaces;

      /// A refinement manager that should be used on the mesh. It is used to
      /// refine elements at interior boundaries in order to fit together with
      /// elements on the other side of the interior boundary.
      RefinementManager* refineManager;

      /// Pointer to a mesh partitioner that is used to partition the mesh to
      /// the ranks.
      MeshPartitioner* partitioner;

      /// Pointer to a mesh partitioner that is used for the very first
      /// partitioning of the mesh. In most cases, this pointer points to the
      /// same object as \ref partitioner, but this must not be the case in
      /// general.
      MeshPartitioner* initialPartitioner;

      /// Weights for the elements, i.e., the number of leaf elements within
      /// this element.
      std::map<int, double> elemWeights;

      /// Stores to every macro element index the number of the rank that owns this
      /// macro element.
      std::map<int, int> partitionMap;

      /// Database to store and query all sub-objects of all elements of the
      /// macro mesh.
      ElementObjectDatabase elObjDb;

      /// Defines the interior boundaries of the domain that result from
      /// partitioning the whole mesh.
      MultiLevelInteriorBoundary intBoundary;

      /// Dof communicator objects for each mesh
      std::map<Mesh*, MultiLevelDofComm> dofComms;

      PeriodicMap periodicMap;

      /// This set of values must be interchanged between ranks when the mesh is
      /// repartitioned.
      std::vector<DOFVector<double>*> interchangeVectors;

      /// If the problem definition has been read from a serialization file, this
      /// variable is true, otherwise it is false. This variable is used to stop the
      /// initialization function, if the problem definition has already been read
      /// from a serialization file.
      bool deserialized;

      /// Denotes whether there exists a filewriter for this object.
      bool writeSerializationFile;

      /// If true, it is possible to repartition the mesh during computations.
      bool repartitioningAllowed;

      /// repartition the mesh (only) the first time repartitionMesh() is called
      bool repartitionOnlyOnce;

      /// Stores the number of mesh changes that must lie in between two
      /// repartitionings.
      int repartitionIthChange;

      ///
      int repartitioningWaitAfterFail;

      /// Counts the number of mesh changes after the last mesh repartitioning
      /// was done.
      int nMeshChangesAfterLastRepartitioning;

      /// Countes the number of mesh repartitions that were done. Till now, this
      /// variable is used only for debug outputs.
      int repartitioningCounter;

      /// If repartitioning of the mesh fail, this variable has a positive value
      /// that gives the number of mesh changes the mesh distributer will wait
      /// before trying new mesh repartitioning.
      int repartitioningFailed;

      /// Directory name where all debug output files should be written to.
      std::string debugOutputDir;

      /// Stores the mesh change index. This is used to recognize changes in the
      /// mesh structure (e.g. through refinement or coarsening managers).
      std::map<Mesh*, long> lastMeshChangeIndexs;

      /// Stores for all macro elements of the original macro mesh the
      /// neighbourhood information based on element indices. Thus, each macro
      /// element index is mapped to a vector containing all indices of
      /// neighbouring macro elements.
      std::map<int, std::vector<int>> macroElementNeighbours;

      /// Store all macro elements of the overall mesh, i.e., before the
      /// mesh is redistributed for the first time.
      /// Store all macro elements of the overall mesh, i.e., before the
      /// mesh is redistributed for the first time.
      std::map<Mesh*, std::vector<MacroElement*>> allMacroElements;

      Flag createBoundaryDofFlag;

      /// Stores on each mesh level for all FE spaces the information about
      /// all boundary DOFs.
      std::vector<std::map<const FiniteElemSpace*, BoundaryDofInfo>> boundaryDofInfo;

      /// Stores information about hierarchical decomposition of the mesh levels.
      /// Is used to specify multi level solver methods.
      MeshLevelData levelData;

      /// If there is no mesh adaptivity, the mesh distributor can remove some
      /// data structures which are only used if mesh changes or it must be
      /// redistributed due to some local adaptivity. By default, this variable
      /// is set to true, and thus no special assumption are made.
      bool meshAdaptivity;

      /// Specifies whether the global domain has periodic boundaries. Thus, this
      /// variable is not related to rank's subdomain but to the global problem
      /// and therefore the value if the same on all ranks.
      bool hasPeriodicBoundary;

      /// Set of all parallel DOF mapping object that are registered by parallel
      /// solver objects and must be updated automatically after mesh change.
      std::vector<ParallelDofMapping*> dofMaps;

      /// If true, detailed timings for benchmarking are printed.
      bool printTimings;

      /// If true, detailed information about memory usage are printed.
      bool printMemoryUsage;

    public:
      /// The boundary DOFs are sorted by subobject entities, i.e., first all
      /// face DOFs, edge DOFs and to the last vertex DOFs will be set to
      /// communication structure vectors, \ref sendDofs and \ref recvDofs.
      static const Flag BOUNDARY_SUBOBJ_SORTED;

      /// When boundary DOFs are created, \ref boundaryDofInfo is filled for
      /// all DOFs that this rank will send to other ranks (thus, rank
      /// owned DOFs.
      static const Flag BOUNDARY_FILL_INFO_SEND_DOFS;

      /// When boundary DOFs are created, \ref boundaryDofInfo is filled for
      /// all DOFs that this rank will receive from other ranks (thus, DOFs
      /// that are owned by another rank).
      static const Flag BOUNDARY_FILL_INFO_RECV_DOFS;

      static MeshDistributor* globalMeshDistributor;

      friend class ParallelDebug;
    };
  }
}

#endif // AMDIS_MESHDISTRIBUTOR_H
