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



/** \file ElementObjectDatabase.h */

#ifndef AMDIS_ELEMENT_OBJECT_DATABASE_H
#define AMDIS_ELEMENT_OBJECT_DATABASE_H

#include <map>
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/container/flat_map.hpp>

//#include "AMDiS_fwd.h"
#include "Containers.h"
#include "Global.h"
#include "Boundary.h"
//#include "Serializer.h"
#include "FiniteElemSpace.h"
#include "Mesh.h"
#include "parallel/ParallelTypes.h"

namespace AMDiS
{
  namespace Parallel
  {

    using boost::container::flat_map;

    /// Just to templatize the typedef.
    template<typename T>
    struct PerBoundMap
    {
      typedef std::map<std::pair<T, T>, BoundaryType> type;
      typedef typename type::iterator iterator;
    };


#pragma pack(push)
#pragma pack(1)
    /// Defines one element object. This may be either a vertex, edge or face.
    struct ElementObjectData
    {
      ElementObjectData(int a = -1, int b = 0)
        : elIndex(a),
          ithObject(b)
      {}

      /// Index of the element this object is part of.
      int elIndex;

      /// Index of the object within the element.
      char ithObject;

      /// Write this element object to disk.
      void serialize(std::ostream& out) const;

      /// Read this element object from disk.
      void deserialize(std::istream& in);

      /// Compare this element object with another one.
      bool operator==(ElementObjectData& cmp) const
      {
        return (elIndex == cmp.elIndex && ithObject == cmp.ithObject);
      }

      /// Define a strict order on element objects.
      bool operator<(const ElementObjectData& rhs) const
      {
        return (elIndex < rhs.elIndex ||
                (elIndex == rhs.elIndex && ithObject < rhs.ithObject));
      }
    };
#pragma pack(pop)


    /** \brief
     * This class is a database of element objects. An element object is either a
     * vertex, edge or the face of a specific element. This database is used to
     * store all objects of all elements of a mesh. The information is stored in a
     * way that makes it possible to identify all elements, which have a given
     * vertex, edge or face in common. If is is known which element is owned by
     * which rank in parallel computations, it is thus possible to get all interior
     * boundaries on object level. This is required, because two elements may share
     * a common vertex without beging neighbours in the definition of AMDiS.
     */
    class ElementObjectDatabase
    {
    public:
      ElementObjectDatabase()
        : feSpace(NULL),
          macroMesh(NULL),
          iterGeoPos(CENTER),
          macroElementRankMap(NULL),
          levelData(NULL)
      {}

      void init(std::vector<Mesh*>&,
                std::map<Mesh*, std::vector<const FiniteElemSpace*>>&);

      Mesh* getMesh()
      {
        return macroMesh;
      }

      /*
       * \param[in]  macroElementRankMap   Maps to each macro element of the mesh
       *                                   the rank that owns this macro element.
       */
      void create(std::map<int, int>& macroElementRankMap,
                  MeshLevelData& levelData);

      void createMacroElementInfo(std::map<Mesh*, std::vector<MacroElement*>>& mel);

      /// Create for a filled object database the membership information for all
      /// element objects. An object is owned by a rank, if the rank has the
      /// heighest rank number of all ranks where the object is part of.
      void updateRankData();

      /// All data from the database is dropped.
      void clear();

      /// Resets iteration;
      bool resetIterator()
      {
        iterGeoPos = CENTER;
        return true; // ???
      }

      /** \brief
       * Iterates over all elements for one geometrical index, i.e., over all
       * vertices, edges or faces in the mesh. The function returns true, if the
       * result is valid. Otherwise the iterator is at the end position.
       *
       * \param[in]  pos   Must be either VERTEX, EDGE or FACE and defines the
       *                   elements that should be traversed.
       */
      inline bool iterate(GeoIndex pos)
      {
        // CENTER marks the variable "iterGeoPos" to be in an undefined state. I.e.,
        // there is no iteration that is actually running.

        if (iterGeoPos == CENTER)
        {
          iterGeoPos = pos;
          switch (iterGeoPos)
          {
          case VERTEX:
            vertexIter = vertexInRank.begin();
            break;
          case EDGE:
            edgeIter = edgeInRank.begin();
            break;
          case FACE:
            faceIter = faceInRank.begin();
            break;
          default:
            ERROR_EXIT("Not GeoIndex %d!\n", iterGeoPos);
          }
        }
        else
        {
          switch (iterGeoPos)
          {
          case VERTEX:
            ++vertexIter;
            break;
          case EDGE:
            ++edgeIter;
            break;
          case FACE:
            ++faceIter;
            break;
          default:
            ERROR_EXIT("Not GeoIndex %d!\n", iterGeoPos);
          }
        }

        switch (iterGeoPos)
        {
        case VERTEX:
          if (vertexIter == vertexInRank.end())
          {
            iterGeoPos = CENTER;
            return false;
          }
          break;
        case EDGE:
          if (edgeIter == edgeInRank.end())
          {
            iterGeoPos = CENTER;
            return false;
          }
          break;
        case FACE:
          if (faceIter == faceInRank.end())
          {
            iterGeoPos = CENTER;
            return false;
          }
          break;
        default:
          ERROR_EXIT("Should not happen!\n");
        }

        return true;
      }


      /// Returns the data of the current iterator position.
      flat_map<int, ElementObjectData>& getIterateData()
      {
        switch (iterGeoPos)
        {
        case VERTEX:
          return vertexIter->second;
          break;
        case EDGE:
          return edgeIter->second;
          break;
        case FACE:
          return faceIter->second;
          break;
        default:
          ERROR_EXIT("Should not happen!\n");

          // Will never be reached, just to avoid compiler warnings.
          return faceIter->second;
        }
      }

      /// Returns the rank owner of the current iterator position.
      int getIterateOwner(int level);

      /// Returns the owner of a macro element vertex.
      int getOwner(DegreeOfFreedom vertex, int level);

      /// Returns the owner of a macro element edge.
      int getOwner(DofEdge edge, int level);

      /// Returns the owner of a macro element face.
      int getOwner(DofFace face, int level);

      /// Returns the rank owner of the current iterator position.
      int getIterateMaxLevel();

      /// Checks if a given vertex DOF is in a given rank.
      int isInRank(DegreeOfFreedom vertex, int rank)
      {
        return (vertexInRank[vertex].count(rank));
      }

      /// Checks if a given edge is in a given rank.
      int isInRank(DofEdge edge, int rank)
      {
        return (edgeInRank[edge].count(rank));
      }

      /// Checks if a given face is in a given rank.
      int isInRank(DofFace face, int rank)
      {
        return (faceInRank[face].count(rank));
      }


      /// Returns a vector with all macro elements that have a given vertex DOF
      /// in common.
      std::vector<ElementObjectData>& getElements(DegreeOfFreedom vertex)
      {
        return vertexElements[vertex];
      }

      /// Returns a vector with all macro elements that have a given edge in common.
      std::vector<ElementObjectData>& getElements(DofEdge edge)
      {
        return edgeElements[edge];
      }

      /// Returns a vector with all macro elements that have a given face in common.
      std::vector<ElementObjectData>& getElements(DofFace face)
      {
        return faceElements[face];
      }


      /// Returns a vector with all macro elements that have a given vertex DOF
      /// in common.
      std::vector<ElementObjectData>& getElementsVertex(int elIndex, int ithVertex)
      {
        ElementObjectData elObj(elIndex, ithVertex);
        DegreeOfFreedom vertex = vertexLocalMap[elObj];
        return vertexElements[vertex];
      }

      /// Returns a vector with all macro elements that have a given edge in common.
      std::vector<ElementObjectData>& getElementsEdge(int elIndex, int ithEdge)
      {
        ElementObjectData elObj(elIndex, ithEdge);
        DofEdge edge = edgeLocalMap[elObj];
        return edgeElements[edge];
      }

      /// Returns a vector with all macro elements that have a given face in common.
      std::vector<ElementObjectData>& getElementsFace(int elIndex, int ithFace)
      {
        ElementObjectData elObj(elIndex, ithFace);
        DofFace face = faceLocalMap[elObj];
        return faceElements[face];
      }

      /// Returns, for a given vertex, a map that maps from rank numbers to
      /// element data objects, which identify on the rank the element which
      /// contains this vertex. If more than one element in one subdomain contains
      /// the vertex, the element with the highest element index is given. If the
      /// vertex is not contained in a rank's subdomain, it will not be considered
      /// in this mapping.
      flat_map<int, ElementObjectData>& getElementsInRank(DegreeOfFreedom vertex)
      {
        return vertexInRank[vertex];
      }

      /// Returns, for a given edge, a map that maps from rank numbers to
      /// element data objects, which identify on the rank the element which
      /// contains this edge. If more than one element in one subdomain contains
      /// the edge, the element with the highest element index is given. If the
      /// edge is not contained in a rank's subdomain, it will not be considered
      /// in this mapping.
      flat_map<int, ElementObjectData>& getElementsInRank(DofEdge edge)
      {
        return edgeInRank[edge];
      }

      /// Returns, for a given face, a map that maps from rank numbers to
      /// element data objects, which identify on the rank the element which
      /// contains this face. If the face is not contained in a rank's subdomain,
      /// it will not be considered in this mapping.
      flat_map<int, ElementObjectData>& getElementsInRank(DofFace face)
      {
        return faceInRank[face];
      }

      /// Get degree of a DOF, thus the number of ranks which contain it.
      inline int getDegree(DegreeOfFreedom dof)
      {
        return vertexInRank[dof].size();
      }

      /// Get degree of an edge, thus the number of ranks which contain it.
      inline int getDegree(DofEdge edge)
      {
        return edgeInRank[edge].size();
      }

      /// Get degree of a face, thus the number of ranks which contain it.
      inline int getDegree(DofFace face)
      {
        return faceInRank[face].size();
      }

      /// Returns to an element object data the appropriate vertex DOF.
      DegreeOfFreedom getVertexLocalMap(ElementObjectData& data)
      {
        TEST_EXIT_DBG(vertexLocalMap.count(data))("Should not happen!\n");

        return vertexLocalMap[data];
      }

      /// Returns to an element object data the appropriate edge.
      DofEdge getEdgeLocalMap(ElementObjectData& data)
      {
        TEST_EXIT_DBG(edgeLocalMap.count(data))("Should not happen!\n");

        return edgeLocalMap[data];
      }

      /// Returns to an element object data the appropriate face.
      DofFace getFaceLocalMap(ElementObjectData& data)
      {
        TEST_EXIT_DBG(faceLocalMap.count(data))("Should not happen!\n");

        return faceLocalMap[data];
      }

      PerBoundMap<DegreeOfFreedom>::type& getPeriodicVertices()
      {
        return periodicVertices;
      }

      PerBoundMap<DofEdge>::type& getPeriodicEdges()
      {
        return periodicEdges;
      }

      PerBoundMap<DofFace>::type& getPeriodicFaces()
      {
        return periodicFaces;
      }

      inline bool getEdgeReverseMode(ElementObjectData& obj0,
                                     ElementObjectData& obj1)
      {
        if (macroMesh->getDim() == 2)
          return true;

        if (edgeReverseMode.empty())
          return false;

        return static_cast<bool>(edgeReverseMode.count(std::make_pair(obj0, obj1)));
      }

      inline bool getFaceReverseMode(ElementObjectData& obj0,
                                     ElementObjectData& obj1)
      {
        if (faceReverseMode.empty())
          return false;

        return static_cast<bool>(faceReverseMode.count(std::make_pair(obj0, obj1)));
      }

      /// Returns true if there is periodic data.
      bool hasPeriodicData()
      {
        return (periodicVertices.size() != 0);
      }

      /// Returns true if the given boundary type is larger or equal to the smallest
      /// periodic boundary ID in mesh. See \ref smallestPeriodicBcType for more
      /// information.
      bool isValidPeriodicType(BoundaryType t) const
      {
        return (t >= smallestPeriodicBcType);
      }

      inline Element* getElementPtr(int index, Mesh* mesh)
      {
        FUNCNAME_DBG("ElementObjectDatabase::getElementPtr()");
        TEST_EXIT_DBG(macroElIndexMap[index][mesh])
        ("No element pointer in macroElIndex map. Something is wrong.\n");
        return macroElIndexMap[index][mesh];
      }

      inline MacroElIndexMap& getElIndexMap()
      {
        return macroElIndexMap;
      }

      inline int getElementType(int index)
      {
        return macroElIndexTypeMap[index];
      }

      void setData(std::map<int, int>& rankMap,
                   MeshLevelData& ld)
      {
        macroElementRankMap = &rankMap;
        levelData = &ld;
      }

      /// Write the element database to disk.
      void serialize(std::ostream& out);

      /// Read the element database from disk.
      void deserialize(std::istream& in);

      /// Returns the estimated memory usage of an object of this class.
      unsigned long calculateMemoryUsage();

    protected:
      /** \brief
       * Adds an element to the object database. If the element is part of a
       * periodic boundary, all information about subobjects of the element on
       * this boundary are collected.
       *
       * \param[in]  Element    the element.
       */
      void addElement(Element* el);

      void addElementPeriodicBoundary(ElInfo* elInfo);

      /// Adds the i-th DOF vertex of an element to the object database.
      void addVertex(Element* el, int ith);

      /// Adds the i-th edge of an element to the object database.
      void addEdge(Element* el, int ith);

      /// Adds the i-th face of an element to the object database.
      void addFace(Element* el, int ith);

      /// Creates final data of the periodic boundaries. Must be called after all
      /// elements of the mesh are added to the object database. Then this
      /// functions search for indirectly connected vertices in periodic
      /// boundaries. This is only the case, if there are more than one boundary
      /// conditions. Then, e.g., in 2D, all edges of a square are iterectly
      /// connected. In 3D, if the macro mesh is a box, all eight vertex nodes and
      /// always four of the 12 edges are indirectly connected.
      void createPeriodicData();

      /// Creates on all boundaries the reverse mode flag.
      void createReverseModeData();

      BoundaryType getNewBoundaryType();

      BoundaryType provideConnectedPeriodicBoundary(BoundaryType b0,
          BoundaryType b1);

      /// Some auxiliary function to write the element object database to disk.
      void serialize(std::ostream& out, std::vector<ElementObjectData>& elVec);

      /// Some auxiliary function to read the element object database from disk.
      void deserialize(std::istream& in, std::vector<ElementObjectData>& elVec);

      /// Some auxiliary function to write the element object database to disk.
      void serialize(std::ostream& out, flat_map<int, ElementObjectData>& data);

      /// Some auxiliary function to read the element object database from disk.
      void deserialize(std::istream& in, flat_map<int, ElementObjectData>& data);

      int getOwner(std::vector<ElementObjectData>& objData, int level);
    private:
      const FiniteElemSpace* feSpace;

      /// The macro mesh that is used to store all its element information in the database.
      Mesh* macroMesh;

      /// Used to get mesh or element pointers. meshes[0] is always equal to \ref macroMesh
      std::vector<Mesh*> meshes;

      /// Corresponding FE spaces of meshes.
      std::vector<const FiniteElemSpace*> feSpaces;


      /// Maps to each vertex DOF all element objects that represent this vertex.
      flat_map<DegreeOfFreedom, std::vector<ElementObjectData>> vertexElements;

      /// Maps to each edge all element objects that represent this edge.
      flat_map<DofEdge, std::vector<ElementObjectData>> edgeElements;

      /// Maps to each face all element objects that represent this edge.
      flat_map<DofFace, std::vector<ElementObjectData>> faceElements;


      /// Temporary object to speed up creation of \ref vertexElements
      std::map<DegreeOfFreedom, std::vector<ElementObjectData>> tmpVertexElements;

      /// Temporary object to speed up creation of \ref edgeElements
      std::map<DofEdge, std::vector<ElementObjectData>> tmpEdgeElements;

      /// Temporary object to speed up creation of \ref faceElements
      std::map<DofFace, std::vector<ElementObjectData>> tmpFaceElements;



      /// Maps to an element object the corresponding vertex DOF.
      flat_map<ElementObjectData, DegreeOfFreedom> vertexLocalMap;

      /// Maps to an element object the corresponding edge.
      flat_map<ElementObjectData, DofEdge> edgeLocalMap;

      /// Maps to an element object the corresponding face.
      flat_map<ElementObjectData, DofFace> faceLocalMap;


      /// Defines to each vertex DOF a map that maps to each rank number the element
      /// objects that have this vertex DOF in common.
      flat_map<DegreeOfFreedom, flat_map<int, ElementObjectData>> vertexInRank;

      /// Defines to each edge a map that maps to each rank number the element
      /// objects that have this edge in common.
      flat_map<DofEdge, flat_map<int, ElementObjectData>> edgeInRank;

      /// Defines to each face a map that maps to each rank number the element
      /// objects that have this face in common.
      flat_map<DofFace, flat_map<int, ElementObjectData>> faceInRank;


      /// Vertex iterator to iterate over \ref vertexInRank
      flat_map<DegreeOfFreedom, flat_map<int, ElementObjectData>>::iterator vertexIter;

      /// Edge iterator to iterate over \ref edgeInRank
      flat_map<DofEdge, flat_map<int, ElementObjectData>>::iterator edgeIter;

      /// Face iterator to iterate over \ref faceInRank
      flat_map<DofFace, flat_map<int, ElementObjectData>>::iterator faceIter;


      /// Defines the geometrical iteration index of the iterators. I.e., the value
      /// is either VERTEX, EDGE or FACE and the corresponding element objects are
      /// traversed. The value CENTER is used to define a not defined states of the
      /// iterators, i.e., if no iteration is running.
      GeoIndex iterGeoPos;

      std::map<std::pair<BoundaryType, BoundaryType>, BoundaryType> bConnMap;

      /// The following three data structures store periodic DOFs, edges and faces.
      PerBoundMap<DegreeOfFreedom>::type periodicVertices;
      PerBoundMap<DofEdge>::type periodicEdges;
      PerBoundMap<DofFace>::type periodicFaces;

      /// Defines the smallest boudary ID for periodic boundary conditions. This is
      /// required to distinguish between "real" periodic boundaries and periodic
      /// boundary IDs that are set by the parallel algorithm for indirectly
      /// connected boundaries.
      BoundaryType smallestPeriodicBcType;

      /// Stores to each vertex all its periodic associations.
      std::map<DegreeOfFreedom, std::set<BoundaryType>> periodicDofAssoc;

      /// Stores to each edge all its periodic associations.
      std::map<DofEdge, std::set<DofEdge>> periodicEdgeAssoc;

      /// Stores all interior edge boundaries which have reverse mode enabled.
      std::set<std::pair<ElementObjectData, ElementObjectData>> edgeReverseMode;

      /// Stores all interior face boundaries which have reverse mode enabled.
      std::set<std::pair<ElementObjectData, ElementObjectData>> faceReverseMode;

      std::map<int, int>* macroElementRankMap;

      /// Maps to each macro element index a pointer to the corresponding element.
      MacroElIndexMap macroElIndexMap;

      /// Maps to each macro element index the type of this element.
      flat_map<int, int> macroElIndexTypeMap;

      MeshLevelData* levelData;

      friend class ParallelDebug;
    };

  }
}

#endif
