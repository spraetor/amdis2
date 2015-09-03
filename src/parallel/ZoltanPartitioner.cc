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

#ifdef HAVE_ZOLTAN

#include "parallel/ZoltanPartitioner.h"
#include "Traverse.h"
#include "ElInfo.h"
#include "Initfile.h"

using namespace std;

namespace AMDiS
{
  namespace Parallel
  {

    ZoltanPartitioner::ZoltanPartitioner(string name,
                                         MPI::Intracomm* comm)
      : MeshPartitioner(name, comm),
        zoltan(*comm),
        elWeights(NULL)
    {
      /* Read configuration for Zoltan
       * format in initfile:
       * <tag><NameOfZoltanParameter>:<ValueOfParameter>
       * e.g.
       *  zoltan parameter->LB_METHOD: GRAPH
       */
      std::string tag="zoltan parameter->";
      std::string tagInitial="zoltan parameter->";

      Parameters::get("zoltan->parameter tag", tag);
      Parameters::get("zoltan->parameter initial tag", tagInitial);

      Parameters::getParameterMap(tag,paramMap,2);
      Parameters::getParameterMap(tagInitial,paramMapInitial,2);

    }


    bool ZoltanPartitioner::partition(map<int, double>& weights,
                                      PartitionMode mode)
    {
      FUNCNAME("ZoltanPartitioner::partition()");

      if (!boxPartitioning)
      {
        zoltan.Set_Num_Obj_Fn(ZoltanFunctions::getNumObj, this);
        zoltan.Set_Obj_List_Fn(ZoltanFunctions::getObjectList, this);
        zoltan.Set_Num_Geom_Fn(ZoltanFunctions::getNumGeom, this);
        zoltan.Set_Geom_Multi_Fn(ZoltanFunctions::getGeomMulti, this);
        zoltan.Set_Num_Edges_Multi_Fn(ZoltanFunctions::getNumEdgeMulti, this);
        zoltan.Set_Edge_List_Multi_Fn(ZoltanFunctions::getEdgeListMulti, this);
      }
      else
      {
        zoltan.Set_Num_Obj_Fn(ZoltanFunctions::box_getNumObj, this);
        zoltan.Set_Obj_List_Fn(ZoltanFunctions::box_getObjectList, this);
        zoltan.Set_Num_Geom_Fn(ZoltanFunctions::box_getNumGeom, this);
        //zoltan.Set_Geom_Multi_Fn(ZoltanFunctions::box_getGeomMulti, this);
        zoltan.Set_Num_Edges_Multi_Fn(ZoltanFunctions::box_getNumEdgeMulti, this);
        zoltan.Set_Edge_List_Multi_Fn(ZoltanFunctions::box_getEdgeListMulti, this);
      }

      if (mode != INITIAL)
        elWeights = &weights;

      int changes;
      int nGid, nLid;

      int nImportEls;
      ZOLTAN_ID_PTR import_global_ids, import_local_ids;
      int* import_procs;
      int* import_to_part;

      int nExportEls;
      ZOLTAN_ID_PTR export_global_ids, export_local_ids;
      int* export_procs;
      int* export_to_part;

      /*
       *  Set default configuration for Zoltan
       */

      if (mode == INITIAL)
      {
        zoltan.Set_Param("LB_APPROACH", "PARTITION");
        zoltan.Set_Param("LB_METHOD", "GRAPH");
        zoltan.Set_Param("REDUCE_DIMENSIONS", "1");
        zoltan.Set_Param("DEGENERATE_RATIO", "1.1");
        zoltan.Set_Param("RCB_RECTILINEAR_BLOCKS", "1");
        zoltan.Set_Param("AVERAGE_CUTS", "1");
        zoltan.Set_Param("RCB_RECOMPUTE_BOX", "1");
        if (boxPartitioning)
          zoltan.Set_Param("LB_METHOD", "GRAPH");
      }
      else
      {
        zoltan.Set_Param("LB_APPROACH", "REPARTITION");
        zoltan.Set_Param("LB_METHOD", "GRAPH");
        zoltan.Set_Param("REFTREE_INITPATH", "CONNECTED");
        zoltan.Set_Param("REDUCE_DIMENSIONS", "1");
        zoltan.Set_Param("DEGENERATE_RATIO", "1.1");
        zoltan.Set_Param("RCB_RECTILINEAR_BLOCKS", "1");
        zoltan.Set_Param("AVERAGE_CUTS", "1");
        zoltan.Set_Param("RCB_RECOMPUTE_BOX", "1");
      }

      /*
       * Overwrite  default config of zoltan with Values
       * defined in init-file.
       */
      if (mode == INITIAL)
      {
        std::map<string,string>::const_iterator itr;
        for(itr=paramMapInitial.begin(); itr!= paramMapInitial.end(); ++itr)
          if( zoltan.Set_Param((*itr).first, (*itr).second) != ZOLTAN_OK)
          {
            ERROR_EXIT("Wrong parameter for Zoltan in Initfile (paramInitialMap): %s : %s \n", (*itr).first.c_str(), (*itr).second.c_str());
          }
      }
      else
      {
        std::map<string,string>::const_iterator itr;
        for(itr=paramMap.begin(); itr!= paramMap.end(); ++itr)
          if( zoltan.Set_Param((*itr).first, (*itr).second) != ZOLTAN_OK)
          {
            ERROR_EXIT("Wrong parameter for Zoltan in Initfile (paramMap): %s : %s \n", (*itr).first.c_str(), (*itr).second.c_str());
          }
      }



      zoltan.Set_Param("OBJ_WEIGHT_DIM", "1");

#if (DEBUG != 0)
      zoltan.Set_Param("DEBUG_LEVEL", "1");
#else
      zoltan.Set_Param("DEBUG_LEVEL", "0");
#endif


      int err = zoltan.LB_Partition(changes, nGid, nLid,
                                    nImportEls,
                                    import_global_ids,
                                    import_local_ids,
                                    import_procs,
                                    import_to_part,
                                    nExportEls,
                                    export_global_ids,
                                    export_local_ids,
                                    export_procs,
                                    export_to_part);

      recvElements.clear();
      sendElements.clear();


      // Check if new partitioning would lead to a empty partition
      int createEmptyPartition;
      if (!boxPartitioning)
      {
        createEmptyPartition=
          (mesh->getMacroElements().size() == nExportEls && nImportEls == 0);
      }
      else
      {
        int realNrOfExpEl=0;
        for (int i = 0; i < nExportEls; i++)
        {
          int sendElIndex = export_global_ids[i];
          realNrOfExpEl+=boxSplitting[sendElIndex].size();
        }
        createEmptyPartition=
          (mesh->getMacroElements().size() == realNrOfExpEl && nImportEls == 0);
      }

      mpi::globalMax(createEmptyPartition);

      if (createEmptyPartition > 0)
        err = ZOLTAN_FATAL;

      // if a valid partition was produced, collect all elements to be send or received
      // in a list
      if (err == ZOLTAN_OK && changes != 0)
      {
        if (nImportEls > 0)
        {
          for (int i = 0; i < nImportEls; i++)
          {
            int recvElIndex = import_global_ids[i];

            if (!boxPartitioning)
            {
              elementInRank[recvElIndex] = true;
              recvElements[import_procs[i]].push_back(recvElIndex);
            }
            else
            {
              set<int>& box = boxSplitting[recvElIndex];

              for (set<int>::iterator it = box.begin(); it != box.end(); ++it)
              {
                elementInRank[*it] = true;
                recvElements[import_procs[i]].push_back(*it);
              }
            }
          }
        }

        if (nExportEls > 0)
        {
          for (int i = 0; i < nExportEls; i++)
          {
            int sendElIndex = export_global_ids[i];

            if (!boxPartitioning)
            {
              elementInRank[sendElIndex] = false;
              sendElements[export_procs[i]].push_back(sendElIndex);
            }
            else
            {
              set<int>& box = boxSplitting[sendElIndex];

              for (set<int>::iterator it = box.begin(); it != box.end(); ++it)
              {
                elementInRank[*it] = false;
                sendElements[export_procs[i]].push_back(*it);
              }
            }
          }
        }

        createPartitionMap();
      }

      zoltan.LB_Free_Part(&import_global_ids, &import_local_ids,
                          &import_procs, &import_to_part);
      zoltan.LB_Free_Part(&export_global_ids, &export_local_ids,
                          &export_procs, &export_to_part);
      elWeights = NULL;

      if (err != ZOLTAN_OK)
        return false;

      return true;
    }


    void ZoltanPartitioner::createPartitionMap()
    {
      FUNCNAME("ZoltanPartitioner::getPartitionMap()");

      int mpiSize = mpiComm->Get_size();

      vector<int> localElements;
      for (map<int, bool>::iterator it = elementInRank.begin();
           it != elementInRank.end(); ++it)
        if (it->second)
          localElements.push_back(it->first);

      int nLocalElements = localElements.size();
      vector<int> nPartitionElements(mpiSize);
      vector<int> elDist(mpiSize + 1);
      mpiComm->Allgather(&nLocalElements, 1, MPI_INT, &(elDist[1]), 1, MPI_INT);

      elDist[0] = 0;
      nPartitionElements[0] = elDist[1];

      for (int i = 2; i <= mpiSize; i++)
      {
        nPartitionElements[i - 1] = elDist[i];
        elDist[i] += elDist[i - 1];
      }

      int nOverallElements = elDist[mpiSize];

      TEST_EXIT_DBG(nOverallElements == static_cast<int>(elementInRank.size()))
      ("Number of elements differs: %d %d!\n", nOverallElements, elementInRank.size());

      vector<int> partitionElements(nOverallElements);

      // distribute partition elements
      mpiComm->Allgatherv(&(localElements[0]),
                          nLocalElements,
                          MPI_INT,
                          &(partitionElements[0]),
                          &(nPartitionElements[0]),
                          &(elDist[0]),
                          MPI_INT);

      // fill partitionMap
      partitionMap.clear();
      for (int i = 0; i < mpiSize; i++)
        for (int j = 0; j < nPartitionElements[i]; j++)
          partitionMap[partitionElements[elDist[i] + j]] = i;
    }



    // === Zoltan query function for "normal" element partitioning. ===


    int ZoltanFunctions::getNumObj(void* data, int* ierr)
    {
      ZoltanPartitioner* zoltan = (ZoltanPartitioner*)data;
      map<int, bool>& elInRank = zoltan->getElementInRank();

      int nObjects = 0;
      for (map<int, bool>::iterator it = elInRank.begin(); it != elInRank.end(); ++it)
        if (it->second)
          nObjects++;

      *ierr = ZOLTAN_OK;

      return nObjects;
    }


    void ZoltanFunctions::getObjectList(void* data,
                                        int sizeGid,
                                        int sizeLid,
                                        ZOLTAN_ID_PTR globalId,
                                        ZOLTAN_ID_PTR localId,
                                        int wgt_dim,
                                        float* obj_wgts,
                                        int* ierr)
    {
      ZoltanPartitioner* zoltan = (ZoltanPartitioner*)data;
      map<int, bool>& elInRank = zoltan->getElementInRank();
      int localCounter = 0;

      for (map<int, bool>::iterator it = elInRank.begin(); it != elInRank.end(); ++it)
      {
        if (it->second)
        {
          globalId[localCounter] = it->first;
          localId[localCounter] = localCounter;

          if (wgt_dim > 0)
          {
            TEST_EXIT_DBG(wgt_dim == 1)("Should not happen!\n");

            if (zoltan->elWeights)
            {
              TEST_EXIT_DBG(zoltan->elWeights->count(it->first))
              ("No element weight for element index %d\n", it->first);

              obj_wgts[localCounter] = static_cast<float>((*zoltan->elWeights)[it->first]);
            }
            else
            {
              obj_wgts[localCounter] = 1.0;
            }
          }

          localCounter++;
        }
      }

      *ierr = ZOLTAN_OK;
    }


    int ZoltanFunctions::getNumGeom(void* data, int* ierr)
    {
      ZoltanPartitioner* zoltan = (ZoltanPartitioner*)data;
      *ierr = ZOLTAN_OK;

      return zoltan->getMesh()->getGeo(WORLD);
    }


    void ZoltanFunctions::getGeomMulti(void* data,
                                       int num_gid_entries,
                                       int num_lid_entries,
                                       int num_obj,
                                       ZOLTAN_ID_PTR global_ids,
                                       ZOLTAN_ID_PTR local_ids,
                                       int num_dim,
                                       double* geom_vec,
                                       int* ierr)
    {
      FUNCNAME("ZoltanFunctions::getGeomMulti()");

      ZoltanPartitioner* zoltan = (ZoltanPartitioner*)data;
      Mesh* mesh = zoltan->getMesh();

      TEST_EXIT_DBG(num_dim == mesh->getGeo(WORLD))("Should not happen!\n");

      map<int, WorldVector<double>> elIndexToCoords;
      DimVec<double> bary(mesh->getDim(), 1.0 / mesh->getGeo(VERTEX));
      WorldVector<double> coords;
      TraverseStack stack;
      ElInfo* elInfo =
        stack.traverseFirst(mesh, 0, Mesh::CALL_EL_LEVEL | Mesh::FILL_COORDS);
      while (elInfo)
      {
        elInfo->coordToWorld(bary, coords);
        elIndexToCoords[elInfo->getElement()->getIndex()] = coords;

        elInfo = stack.traverseNext(elInfo);
      }

      int c = 0;
      for (int i = 0; i < num_obj; i++)
      {
        int elIndex = global_ids[i];

        TEST_EXIT_DBG(elIndexToCoords.count(elIndex))("Should not happen!\n");

        for (int j = 0; j < num_dim; j++)
          geom_vec[c++] = elIndexToCoords[elIndex][j];
      }

      *ierr = ZOLTAN_OK;
    }


    void ZoltanFunctions::getNumEdgeMulti(void* data,
                                          int num_gid_entries,
                                          int num_lid_entries,
                                          int num_obj,
                                          ZOLTAN_ID_PTR global_ids,
                                          ZOLTAN_ID_PTR local_ids,
                                          int* num_edges,
                                          int* ierr)
    {
      FUNCNAME("ZoltanFunctions::getNumEdgeMulti()");

      ZoltanPartitioner* zoltan = (ZoltanPartitioner*)data;

      for (int i = 0; i < num_obj; i++)
      {
        TEST_EXIT_DBG(zoltan->elNeighbours.count(global_ids[i]))("Should not happen!\n");

        num_edges[i] = zoltan->elNeighbours[global_ids[i]].size();
      }

      *ierr = ZOLTAN_OK;
    }


    void ZoltanFunctions::getEdgeListMulti(void* data,
                                           int num_gid_entries,
                                           int num_lid_entries,
                                           int num_obj,
                                           ZOLTAN_ID_PTR global_ids,
                                           ZOLTAN_ID_PTR local_ids,
                                           int* num_edges,
                                           ZOLTAN_ID_PTR nbor_global_id,
                                           int* nbor_procs,
                                           int wgt_dim,
                                           float* ewgts,
                                           int* ierr)
    {
      FUNCNAME("ZoltanFunctions::getEdgeListMulti()");

      ZoltanPartitioner* zoltan = (ZoltanPartitioner*)data;

      int c = 0;
      for (int i = 0; i < num_obj; i++)
      {
        TEST_EXIT_DBG(static_cast<int>(zoltan->elNeighbours[global_ids[i]].size()) ==
                      num_edges[i])
        ("Wrong sizes for global el index %d: %d %d\n",
         global_ids[i],
         zoltan->elNeighbours[global_ids[i]].size(),
         num_edges[i]);

        for (int j = 0; j < num_edges[i]; j++)
        {
          int neighIndex = zoltan->elNeighbours[global_ids[i]][j];
          nbor_global_id[c] = neighIndex;
          nbor_procs[c] = zoltan->partitionMap[neighIndex];
          c++;
        }
      }


      TEST_EXIT_DBG(wgt_dim == 0)("Edge weights are not handled yet!\n");

      *ierr = ZOLTAN_OK;
    }



    // === Zoltan query functions for box partitioning. ===

    int ZoltanFunctions::box_getNumObj(void* data, int* ierr)
    {
      FUNCNAME("ZoltanFunctions::box_getNumObj()");

      ZoltanPartitioner* zoltan = (ZoltanPartitioner*)data;
      map<int, bool>& elInRank = zoltan->getElementInRank();

      int nObjects = 0;
      for (map<int, bool>::iterator it = elInRank.begin(); it != elInRank.end(); ++it)
        if (it->second)
          nObjects++;

      TEST_EXIT_DBG(nObjects % 6 == 0)("Should not happen!\n");

      nObjects = nObjects / 6;

      *ierr = ZOLTAN_OK;

      return nObjects;
    }


    void ZoltanFunctions::box_getObjectList(void* data,
                                            int sizeGid,
                                            int sizeLid,
                                            ZOLTAN_ID_PTR globalId,
                                            ZOLTAN_ID_PTR localId,
                                            int wgt_dim,
                                            float* obj_wgts,
                                            int* ierr)
    {
      ZoltanPartitioner* zoltan = (ZoltanPartitioner*)data;
      map<int, bool>& elInRank = zoltan->getElementInRank();
      set<int> boxes;
      map<int, double> boxWeight;

      for (map<int, bool>::iterator it = elInRank.begin(); it != elInRank.end(); ++it)
        if (it->second)
        {
          int boxId = zoltan->elInBox[it->first];
          boxes.insert(boxId);

          if (zoltan->elWeights)
          {
            TEST_EXIT_DBG(zoltan->elWeights->count(it->first))
            ("No element weight for element index %d\n", it->first);

            boxWeight[boxId] += (*zoltan->elWeights)[it->first];
          }
        }

      int localCounter = 0;
      for (set<int>::iterator it = boxes.begin(); it != boxes.end(); ++it)
      {
        globalId[localCounter] = *it;
        localId[localCounter] = localCounter;

        if (wgt_dim > 0)
        {
          TEST_EXIT_DBG(wgt_dim == 1)("Should not happen!\n");

          if (zoltan->elWeights)
            obj_wgts[localCounter] = static_cast<float>(boxWeight[*it]);
          else
            obj_wgts[localCounter] = 1.0;
        }

        localCounter++;
      }

      *ierr = ZOLTAN_OK;
    }


    int ZoltanFunctions::box_getNumGeom(void* data, int* ierr)
    {
      ZoltanPartitioner* zoltan = (ZoltanPartitioner*)data;
      *ierr = ZOLTAN_OK;

      return zoltan->getMesh()->getGeo(WORLD);
    }


    void ZoltanFunctions::box_getGeomMulti(void* data,
                                           int num_gid_entries,
                                           int num_lid_entries,
                                           int num_obj,
                                           ZOLTAN_ID_PTR global_ids,
                                           ZOLTAN_ID_PTR local_ids,
                                           int num_dim,
                                           double* geom_vec,
                                           int* ierr)
    {
      FUNCNAME("ZoltanFunctions::box_getGeomMulti()");

      ERROR_EXIT("Not yet supported!\n");
    }


    void ZoltanFunctions::box_getNumEdgeMulti(void* data,
        int num_gid_entries,
        int num_lid_entries,
        int num_obj,
        ZOLTAN_ID_PTR global_ids,
        ZOLTAN_ID_PTR local_ids,
        int* num_edges,
        int* ierr)
    {
      FUNCNAME("ZoltanFunctions::box_getNumEdgeMulti()");

      ZoltanPartitioner* zoltan = (ZoltanPartitioner*)data;

      for (int i = 0; i < num_obj; i++)
      {
        TEST_EXIT_DBG(zoltan->boxNeighbours.count(global_ids[i]))("Should not happen!\n");

        num_edges[i] = zoltan->boxNeighbours[global_ids[i]].size();
      }

      *ierr = ZOLTAN_OK;
    }


    void ZoltanFunctions::box_getEdgeListMulti(void* data,
        int num_gid_entries,
        int num_lid_entries,
        int num_obj,
        ZOLTAN_ID_PTR global_ids,
        ZOLTAN_ID_PTR local_ids,
        int* num_edges,
        ZOLTAN_ID_PTR nbor_global_id,
        int* nbor_procs,
        int wgt_dim,
        float* ewgts,
        int* ierr)
    {
      FUNCNAME("ZoltanFunctions::box_getEdgeListMulti()");

      ZoltanPartitioner* zoltan = (ZoltanPartitioner*)data;

      int c = 0;
      for (int i = 0; i < num_obj; i++)
      {
        TEST_EXIT_DBG(static_cast<int>(zoltan->boxNeighbours[global_ids[i]].size()) ==
                      num_edges[i])
        ("Wrong sizes for global el index %d: %d %d\n",
         global_ids[i],
         zoltan->boxNeighbours[global_ids[i]].size(),
         num_edges[i]);

        set<int>& neigh = zoltan->boxNeighbours[global_ids[i]];
        for (set<int>::iterator it = neigh.begin(); it != neigh.end(); ++it)
        {
          nbor_global_id[c] = *it;
          nbor_procs[c] = zoltan->partitionMap[*(zoltan->boxSplitting[*it].begin())];
          c++;
        }
      }

      TEST_EXIT_DBG(wgt_dim == 0)("Edge weights are not handled yet!\n");

      *ierr = ZOLTAN_OK;
    }


  }
}

#endif
