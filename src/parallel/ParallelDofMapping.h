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



/** \file FeSpaceMapping.h */

#ifndef AMDIS_FE_SPACE_MAPPING_H
#define AMDIS_FE_SPACE_MAPPING_H

#include <mpi.h>
#include <vector>
#include <map>
#include <set>
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>

#ifndef HAVE_PARALLEL_MTL4
#include <petsc.h>
#include <petscis.h>
#endif

#include "AMDiS_fwd.h"
#include "parallel/DofComm.h"
#include "parallel/MpiHelper.h"
#include "parallel/ParallelTypes.h"
#include "parallel/StdMpi.h"

namespace AMDiS
{
  namespace Parallel
  {

    /** \brief
     * Is used to store matrix indices to all DOFs in rank's subdomain. Thus, the
     * class defines a mapping from component number and DOF index to a global
     * matrix index. This class does not calculate the indices by itself!
     */
    class DofToMatIndex
    {
    public:
      typedef boost::container::flat_map<DegreeOfFreedom, int> MapType;

      DofToMatIndex() {}

      /// Reset the data structure.
      inline void clear()
      {
        data.clear();
      }

      inline void clear(int component)
      {
        data[component].clear();
      }

      /** Add a new mapping for a given DOF.
       *
       * \param[in]   component       Component number for which the mapping
       *                              is defined.
       * \param[in]   dof             DOF index
       * \param[in]   globalMatIndex  Global matrix index.
       */
      inline void add(int component, DegreeOfFreedom dof, int globalMatIndex)
      {
        data[component][dof] = globalMatIndex;
      }

      /// Maps a global DOF index to the global matrix index for a specific
      /// system component number.
      inline int get(int component, DegreeOfFreedom dof) const
      {
        FUNCNAME_DBG("DofToMatIndex::get()");

        std::map<int, MapType>::const_iterator it_component = data.find(component);
        TEST_EXIT_DBG(it_component != data.end())
        ("No mapping data for component %d available!\n", component);

        MapType::const_iterator it_dof = (it_component->second).find(dof);
        TEST_EXIT_DBG(it_dof != (it_component->second).end())
        ("Mapping for DOF %d in component %d does not exist!\n",
         dof, component);

        return it_dof->second;
      }

      /// Returns the number of DOF mappings in one component
      inline int getSize(int component) const
      {
        FUNCNAME_DBG("DofToMatIndex::getSize()");

        std::map<int, MapType>::const_iterator it_component = data.find(component);
        TEST_EXIT_DBG(it_component != data.end())
        ("No mapping data for component %d available!\n", component);
        return (it_component->second).size();
      }

      /// Returns the whole mapping for one component
      inline MapType& getData(int component)
      {
        FUNCNAME_DBG("DofToMatIndex::getData()");

        std::map<int, MapType>::iterator it_component = data.find(component);
        TEST_EXIT_DBG(it_component != data.end())
        ("No mapping data for component %d available!\n", component);
        return it_component->second;
      }

      /// Returns for a given matrix index the component and (local or global) DOF
      /// index. As the data structure is not made for this kind of reverse
      /// search, this is very slow and should be only used for debugging.
      void getReverse(int rowIndex, int& component, int& dofIndex) const;

    private:
      /// The mapping data. For each system component there is a specific map that
      /// maps global DOF indices to global matrix indices.
      std::map<int, MapType> data;
    };


    /// This class defines the parallel mapping of DOFs for one FE space. It is
    /// used by the class \ref ParallelDofMapping to specifiy the mapping for a
    /// set of  FE spaces.
    class ComponentDofMap
    {
    public:
      /// Constructor.
      ComponentDofMap();

      /// Clears all data of the mapping.
      void clear();

      /// Maps a DOF index to both, the local and global index of the mapping. The
      /// global index must not be set.
      MultiIndex& operator[](DegreeOfFreedom d)
      {
        TEST_EXIT_DBG(dofMap.count(d))("DOF %d is not in map!\n", d);

        return dofMap[d];
      }

      /** \brief
       * Searchs the map for a given DOF. It does not fail, if the DOF is not
       * mapped by this mapping. In this case, it returns false. If the DOF is
       * mapped, the result is stored and the function returns true.
       *
       * \param[in]    dof     DOF index for which a mapping is searched.
       * \param[out]   index   In the case that the DOF is mapped, the result
       *                       is stored here.
       */
      inline bool find(DegreeOfFreedom dof, MultiIndex& index)
      {
        DofMap::iterator it = dofMap.find(dof);
        if (it == dofMap.end())
          return false;

        index = it->second;
        return true;
      }

      /// Inserts a new DOF to rank's mapping. The DOF is assumed to be owend by
      /// the rank.
      inline void insertRankDof(DegreeOfFreedom dof0,
                                DegreeOfFreedom dof1 = -1)
      {
        FUNCNAME("ComponentDofMap::insertRankDof()");

        TEST_EXIT_DBG(dofMap.count(dof0) == 0)
        ("DOF %d is already defined in mapping!\n", dof0);

        dofMap[dof0].local = dof1;
        nLocalDofs++;
        if (dof1 != -1)
          nRankDofs++;
      }

      /// Inserts a new DOF to rank's mapping. The DOF exists in rank's subdomain
      /// but is owned by a different rank, thus it is part of an interior boundary.
      inline void insertNonRankDof(DegreeOfFreedom dof0,
                                   DegreeOfFreedom dof1 = -1)
      {
        FUNCNAME("ComponentDofMap::insertNonRankDof()");

        TEST_EXIT_DBG(dofMap.count(dof0) == 0)
        ("DOF %d is already in mapping!\n", dof0);

        dofMap[dof0].local = dof1;
        nLocalDofs++;
        nonRankDofs.insert(dof0);
      }

      /// Checks if a given DOF is in the DOF mapping.
      bool isSet(DegreeOfFreedom dof)
      {
        return static_cast<bool>(dofMap.count(dof));
      }

      /// Checks if a given DOF is a rank owned DOF of the DOF mapping. The DOF
      /// must be a DOF of the mapping (this is not checked here), otherwise the
      /// result is meaningsless.
      bool isRankDof(DegreeOfFreedom dof)
      {
        return !(static_cast<bool>(nonRankDofs.count(dof)));
      }

      bool isRankGlobalDof(int dof)
      {
        return (dof >= rStartDofs && dof < rStartDofs + nRankDofs);
      }

      /// Returns number of DOFs in the mapping.
      unsigned int size()
      {
        return dofMap.size();
      }

      /// Returns the raw data of the mapping.
      DofMap& getMap()
      {
        return dofMap;
      }

      const FiniteElemSpace* getFeSpace() const
      {
        return feSpace;
      }

      const Mesh* getMesh() const
      {
        return mesh;
      }

      DofComm& getDofComm()
      {
        FUNCNAME("ComponentDofMap::getDofComm()");

        TEST_EXIT_DBG(dofComm)("No DOF communicator object defined!\n");

        return *dofComm;
      }

      DofMap::iterator begin()
      {
        return dofMap.begin();
      }

      DofMap::iterator end()
      {
        return dofMap.end();
      }

      /// Recomputes the mapping.
      void update();

      /// Sets the FE space this mapping corresponds to.
      void setFeSpace(const FiniteElemSpace* fe)
      {
        feSpace = fe;
      }

      void setMesh(Mesh* mesh_)
      {
        mesh = mesh_;
      }

      /// Informs the mapping whether a global index must be computed.
      void setGlobalMapping(bool b)
      {
        globalMapping = b;
      }

      /// Sets the DOF communicator.
      void setDofComm(DofComm& dc)
      {
        dofComm = &dc;
      }

      void setMpiComm(MPI::Intracomm& m)
      {
        mpiComm = m;
      }

    private:
      /// Computes a global mapping from the local one.
      void computeGlobalMapping();

      /// Computes the global indices of all DOFs in the mapping that are not owned
      /// by the rank.
      void computeNonLocalIndices();

    private:
      /// DOF communicator for all DOFs on interior boundaries.
      DofComm* dofComm;

      MPI::Intracomm mpiComm;

      /// The FE space this mapping belongs to. This is used only the get the
      /// correct DOF communicator in \ref dofComm.
      const FiniteElemSpace* feSpace;

      const Mesh* mesh;

      /// Mapping data from DOF indices to local and global indices.
      DofMap dofMap;

      /// Set of all DOFs that are in mapping but are not owned by the rank.
      boost::container::flat_set<DegreeOfFreedom> nonRankDofs;

      /// If true, a global index mapping will be computed for all DOFs.
      bool globalMapping;

    public:
      ///
      int nRankDofs, nLocalDofs, nOverallDofs, rStartDofs;
    };


    /// Abstract iterator interface to access containrs containing values of
    /// type \ref ComponentDofMap.
    class ComponentIterator
    {
    public:
      virtual ComponentDofMap& operator*() = 0;

      virtual ComponentDofMap* operator->() = 0;

      virtual bool end() = 0;

      virtual void next() = 0;

      virtual void reset() = 0;
    };


    /// Abstract interface to acces DOF mapping data for each component. Allows
    /// to hide specific implementations, which allow, e.g., to efficiently map
    /// all components having the same FE space to the same DOF mapping.
    class ComponentDataInterface
    {
    public:
      virtual ~ComponentDataInterface() {};
      /// Access via component number
      virtual ComponentDofMap& operator[](int compNumber) = 0;

      /// Acess via FE space pointer
      virtual ComponentDofMap& operator[](const FiniteElemSpace* feSpace) = 0;

      /// Checks whether the DOF mapping is defined for a specific
      /// component number.
      virtual bool isDefinedFor(int compNumber) const = 0;

      /// Returns iterator which iterates over all DOF mappings.
      virtual ComponentIterator& getIteratorData() = 0;

      /// Returns iterator which iterates over the DOF mappings of all
      /// components. If the data is defined for each FE space and more than
      /// one commponent is defined on the same FE space, the iterative will
      /// also iterative multple times over the same DOF mapping object.
      virtual ComponentIterator& getIteratorComponent() = 0;

      /** \brief
       * Initialization of the object.
       *
       * \param[in]   componentSpaces   Set of the FE spaces for each component.
       * \param[in]   feSpaces          Set of all different FE spaces.
       * \param[in]   globalMapping     Mapping is parallel (true) or only
       *                                local (false).
       */
      virtual void init(std::vector<const FiniteElemSpace*>& componentSpaces,
                        std::vector<const FiniteElemSpace*>& feSpaces,
                        bool globalMapping) = 0;

    protected:
      /// The FE spaces for all components.
      std::vector<const FiniteElemSpace*> componentSpaces;

      /// The set of all FE spaces. It uniquly contains all different FE spaces
      /// from \ref feSpaces.
      std::vector<const FiniteElemSpace*> feSpaces;
    };


    /// This class concretizes the interface class \ref ComponentDataInterface. A
    /// DOF mapping is implemented for each component.
    class ComponentData : public ComponentDataInterface
    {
    public:
      ComponentData()
        : iter(this)
      {}

      /// Returns DOF mapping for a given component number.
      ComponentDofMap& operator[](int compNumber)
      {
        TEST_EXIT_DBG(componentData.count(compNumber))
        ("No data for component %d!\n", compNumber);

        return componentData.find(compNumber)->second;
      }

      /// Just to implement the corresponding virtual function in \ref
      /// ComponentDataInterface. Of course it does not work as we have data for
      /// each component. Thus there may be different mappings for the same
      /// FE space.
      ComponentDofMap& operator[](const FiniteElemSpace* feSpace)
      {
        ERROR_EXIT("FE Space access is not possible for component wise defined DOF mappings\n");
        return componentData.find(0)->second;
      }

      /// Return data iterator.
      ComponentIterator& getIteratorData()
      {
        iter.reset();
        return iter;
      }

      /// Return component iterator.
      ComponentIterator& getIteratorComponent()
      {
        iter.reset();
        return iter;
      }

      /// Checks whether the DOF mapping is defined for a specific
      /// component number.
      bool isDefinedFor(int compNumber) const
      {
        return (static_cast<unsigned int>(compNumber) < componentData.size());
      }

      /// Initialization
      void init(std::vector<const FiniteElemSpace*>& f0,
                std::vector<const FiniteElemSpace*>& f1,
                bool globalMapping);

    protected:

      /// Iterator class to iterate over all parallel DOF mappings.
      class Iterator : public ComponentIterator
      {
      public:
        Iterator(ComponentData* d)
          : data(d),
            componentCounter(-1)
        {}

        ComponentDofMap& operator*()
        {
          return (*data)[componentCounter];
        }

        ComponentDofMap* operator->()
        {
          return &((*data)[componentCounter]);
        }


        bool end()
        {
          return (it == data->componentSpaces.end());
        }

        void next()
        {
          ++it;
          ++componentCounter;

          if (it == data->componentSpaces.end())
            componentCounter = -1;
        }

        void reset()
        {
          it = data->componentSpaces.begin();
          componentCounter = 0;
        }

      protected:
        /// Pointer to data class over which the iterator must iterate.
        ComponentData* data;

        /// Internal iterator of the internal data from \ref ComponentData.
        std::vector<const FiniteElemSpace*>::iterator it;

        /// Component number of current iteration.
        int componentCounter;
      };


      /// Data mapping from component numbers to DOF mapping objects.
      std::map<unsigned int, ComponentDofMap> componentData;

      /// Iterator object.
      Iterator iter;



      friend class Iterator;
    };



    /// This class concretizes the interface class \ref ComponentDataInterface. A
    /// DOF mapping is implemented for each finite element space. Thus, different
    /// components sharing the same FE space are handled by the same DOF mapping.
    class FeSpaceData : public ComponentDataInterface
    {
    public:
      FeSpaceData()
        : iterData(this),
          iterComponent(this)
      {}

      /// Returns DOF mapping for a given component number.
      ComponentDofMap& operator[](int compNumber)
      {
        const FiniteElemSpace* feSpace = componentSpaces[compNumber];
        return componentData.find(feSpace)->second;
      }

      /// Returns DOF mapping for a given FE space.
      ComponentDofMap& operator[](const FiniteElemSpace* feSpace)
      {
        TEST_EXIT_DBG(componentData.count(feSpace))("No data for FE space!\n");;

        return componentData.find(feSpace)->second;
      }

      /// Checks whether the DOF mapping is defined for a specific
      /// component number.
      bool isDefinedFor(int compNumber) const
      {
        const FiniteElemSpace* feSpace = componentSpaces[compNumber];
        return static_cast<bool>(componentData.count(feSpace));
      }

      /// Return data iterator.
      ComponentIterator& getIteratorData()
      {
        iterData.reset();
        return iterData;
      }

      /// Return component iterator.
      ComponentIterator& getIteratorComponent()
      {
        iterComponent.reset();
        return iterComponent;
      }

      /// Initialization
      void init(std::vector<const FiniteElemSpace*>& f0,
                std::vector<const FiniteElemSpace*>& f1,
                bool globalMapping);


    protected:

      /// Iterator class to iterate over all parallel DOF mappings.
      class IteratorData : public ComponentIterator
      {
      public:
        IteratorData(FeSpaceData* d)
          : data(d)
        {}

        ComponentDofMap& operator*()
        {
          return (*data)[*it];
        }

        ComponentDofMap* operator->()
        {
          return &((*data)[*it]);
        }

        bool end()
        {
          return (it == data->feSpaces.end());
        }

        void next()
        {
          ++it;
        }

        void reset()
        {
          it = data->feSpaces.begin();
        }

      protected:
        FeSpaceData* data;

        std::vector<const FiniteElemSpace*>::iterator it;
      };


      /// Iterator class to iterate over all component DOF mappings.
      class IteratorComponent : public ComponentIterator
      {
      public:
        IteratorComponent(FeSpaceData* d)
          : data(d)
        {}

        ComponentDofMap& operator*()
        {
          return (*data)[*it];
        }

        ComponentDofMap* operator->()
        {
          return &((*data)[*it]);
        }

        bool end()
        {
          return (it == data->componentSpaces.end());
        }

        void next()
        {
          ++it;
        }

        void reset()
        {
          it = data->componentSpaces.begin();
        }

      protected:
        FeSpaceData* data;

        std::vector<const FiniteElemSpace*>::iterator it;
      };


      std::map<const FiniteElemSpace*, ComponentDofMap> componentData;

      IteratorData iterData;

      IteratorComponent iterComponent;

      friend class IteratorData;

      friend class IteratorComponent;
    };



    /// Used to specify whether a parallel DOF mapping is defined for each
    /// specific component or for each FE space.
    enum DofMappingMode
    {
      COMPONENT_WISE,
      FESPACE_WISE
    };

    /**
     * Implements the mapping from sets of distributed DOF indices to local and
     * global indices. The mapping works for a given set of FE spaces. Furthermore,
     * this class may compute the matrix indices of the set of DOF indices.
     */
    class ParallelDofMapping
    {
    public:
      /** \brief
       * Constructur for parallel DOF mapping.
       *
       * \param[in]  mode            Defines if DOF mapping is defined either per
       *                             component or per FE space.
       * \param[in]  matIndexGlobal  If true, the mat index is defined on global
       *                             DOF indices, otherwise on local ones.
       */
      ParallelDofMapping(DofMappingMode mode,
                         bool matIndexFromGlobal = false);

      ~ParallelDofMapping()
      {
        if (data)
          delete data;
        data = NULL;
      }

      /** \brief
       * Initialize the parallel DOF mapping.
       *
       * \param[in]  componentSpaces    The FE spaces of all components of the
       *                                PDE to be solved.
       * \param[in]  feSpaces           Unique list of FE spaces. Thus, two
       *                                arbitrary elements of this list are always
       *                                different.
       * \param[in]  globalMapping      If true, at least one rank's mapping con-
       *                                taines DOFs that are not owend by the rank.
       */
      void init(std::vector<const FiniteElemSpace*>& componentSpaces,
                std::vector<const FiniteElemSpace*>& feSpaces,
                bool globalMapping = true);

      /// In the case of having only one FE space, this init function can be used.
      void init(const FiniteElemSpace* feSpace,
                bool globalMapping = true)
      {
        std::vector<const FiniteElemSpace*> feSpaces;
        feSpaces.push_back(feSpace);
        init(feSpaces, feSpaces, globalMapping);
      }

      void setMpiComm(MPI::Intracomm& m);

      inline MPI::Intracomm& getMpiComm()
      {
        return mpiComm;
      }

      /// Clear all data.
      void clear(Mesh* mesh = NULL);

      /// Set the DOF communicator objects that are required to exchange information
      /// about DOFs that are on interior boundaries.

      void setDofComms(std::map<Mesh*, MultiLevelDofComm>& dofComms, int level);

      /// Returns the DOF communicator.
      DofComm& getDofComm(const FiniteElemSpace* feSpace)
      {
        return (*data)[feSpace].getDofComm();
      }

      inline bool isMatIndexFromGlobal()
      {
        return needMatIndexFromGlobal;
      }

      /// Access the DOF mapping for a given component number.
      inline ComponentDofMap& operator[](int compNumber)
      {
        return (*data)[compNumber];
      }

      /// Access the DOF mapping for a given FE space
      inline ComponentDofMap& operator[](const FiniteElemSpace* feSpace)
      {
        return (*data)[feSpace];
      }

      /// Checks whether the DOF mapping is defined for a specific
      /// component number.
      inline bool isDefinedFor(int compNumber) const
      {
        return data->isDefinedFor(compNumber);
      }

      /// Returns the number of solution components the mapping is defined on.
      inline int getNumberOfComponents() const
      {
        return static_cast<int>(componentSpaces.size());
      }

      /// Returns \ref nRankDofs, thus the number of DOFs owned by the rank.
      inline int getRankDofs() const
      {
        TEST_EXIT_DBG(nRankDofs >= 0)("Should not happen!\n");

        return nRankDofs;
      }

      inline int getRankDofs(int component) const
      {
        int nDofs = (*data)[component].nRankDofs;
        TEST_EXIT_DBG(nDofs >= 0)("Should not happen!\n");

        return nDofs;
      }

      /// Returns \ref nLocalDofs, thus the number of DOFs in ranks subdomain.
      inline int getLocalDofs() const
      {
        TEST_EXIT_DBG(nLocalDofs >= 0)("Should not happen!\n");

        return nLocalDofs;
      }

      inline int getLocalDofs(int component) const
      {
        int nDofs = (*data)[component].nLocalDofs;
        TEST_EXIT_DBG(nDofs >= 0)("Should not happen!\n");

        return nDofs;
      }

      /// Returns \ref nOverallDofs, thus the number of all DOFs in this mapping.
      inline int getOverallDofs() const
      {
        TEST_EXIT_DBG(nOverallDofs >= 0)("Should not happen!\n");

        return nOverallDofs;
      }

      inline int getOverallDofs(int component) const
      {
        int nDofs = (*data)[component].nOverallDofs;
        TEST_EXIT_DBG(nDofs >= 0)("Should not happen!\n");

        return nDofs;
      }

      /// Returns \ref rStartDofs, thus the smallest global index of a DOF that is
      /// owned by the rank.
      inline int getStartDofs() const
      {
        TEST_EXIT_DBG(rStartDofs >= 0)("Should not happen!\n");

        return rStartDofs;
      }

      inline int getStartDofs(int component) const
      {
        int nDofs = (*data)[component].rStartDofs;
        TEST_EXIT_DBG(nDofs >= 0)("Should not happen!\n");

        return nDofs;
      }

      /// Update the mapping.
      void update(Mesh* mesh = NULL);

      /// Updates only the DOF to matrix index mapping
      void updateMatIndex();

      inline DofToMatIndex::MapType& getMatData(int component)
      {
        return dofToMatIndex.getData(component);
      }

      /// Returns the global matrix index of a given DOF for a given
      /// component number.
      inline int getMatIndex(int ithComponent, DegreeOfFreedom d) const
      {
        return dofToMatIndex.get(ithComponent, d);
      }

      /// Returns the component number and local/global DOF index for a given
      /// matrix row index. Should be used for debugging only!
      inline void getReverseMatIndex(int index, int& component, int& dofIndex) const
      {
        dofToMatIndex.getReverse(index, component, dofIndex);
      }

      /// Returns the local matrix index of a given DOF for a given
      /// component number.
      inline int getLocalMatIndex(int ithComponent, DegreeOfFreedom d) const
      {
        return dofToMatIndex.get(ithComponent, d) - rStartDofs;
      }

      /// Returns the set of unique FE spaces.
      std::vector<const FiniteElemSpace*> getFeSpaces(Mesh* mesh = NULL)
      {
        if (!mesh)
          return feSpaces;

        std::vector<const FiniteElemSpace*> thisMeshSpaces;
        for (size_t i = 0; i < feSpaces.size(); i++)
          if (feSpaces[i]->getMesh() == mesh)
            thisMeshSpaces.push_back(feSpaces[i]);
        return thisMeshSpaces;
      }

      bool hasFeSpace(const FiniteElemSpace* feSpace)
      {
        return (std::find(feSpaces.begin(), feSpaces.end(), feSpace)
                != feSpaces.end());
      }

      // Writes all data of this object to an output stream.
      void serialize(std::ostream& out)
      {
        ERROR_EXIT("MUST BE IMPLEMENTED!\n");
      }

      // Reads the object data from an input stream.
      void deserialize(std::istream& in)
      {
        ERROR_EXIT("MUST BE IMPLEMENTED!\n");
      }

      /// Compute local and global matrix indices.
      void computeMatIndex(bool globalIndex);

#ifndef HAVE_PARALLEL_MTL4
      void createIndexSet(IS& is,
                          int firstComponent,
                          int nComponents);
#endif

      /// Prints out some information about the mapping. May be used during
      /// debugging or parallel solver creation.
      void printInfo();

    protected:
      /// Compute \ref nRankDofs.
      int computeRankDofs();

      /// Compute \ref nLocalDofs.
      int computeLocalDofs();

      /// Compute \ref nOverallDofs.
      int computeOverallDofs();

      /// Compute \ref rStartDofs.
      int computeStartDofs();

    private:
      MPI::Intracomm mpiComm;

      /// Is true if there are DOFs in at least one subdomain that are not owned
      /// by the rank. If the value is false, each rank contains only DOFs that
      /// are also owned by this rank.
      bool globalMapping;

      /// If matrix indices should be computed, this variable defines if the
      /// mapping from DOF indices to matrix row indices is defined on local
      /// or global DOF indices. If true, the mapping is to specify and to use
      /// on global ones, otherwise on local DOF indices.
      /// In most scenarios the mapping stored on local DOF indices is what we
      /// want to have. Only when periodic boundary conditions are used together
      /// with some global matrix approach, the matrix indices must be stored
      /// also for DOFs that are not part of the local subdomain. Thus, the
      /// mapping will be stored on global DOF indices.
      bool needMatIndexFromGlobal;

      /// Maps from components to DOF mappings.
      ComponentDataInterface* data;

      /// Number of DOFs owned by rank.
      int nRankDofs;

      /// Number of DOFs in rank's subdomain.
      int nLocalDofs;

      /// Number of global DOFs (this value is thus the same on all ranks).
      int nOverallDofs;

      /// Smallest global index of a DOF owned by the rank.
      int rStartDofs;

      /// Mapping from global DOF indices to global matrix indices under
      /// consideration of possibly multiple components.
      DofToMatIndex dofToMatIndex;

      /// FE spaces of all components.
      std::vector<const FiniteElemSpace*> componentSpaces;

      /// Set of unique FE spaces.
      std::vector<const FiniteElemSpace*> feSpaces;

      /// Defines the DOF mapping either. The mapping may be defined either for
      /// FE spaces or for component numbers.
      DofMappingMode mode;
    };
  }
}

#endif
