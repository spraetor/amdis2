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



/** \file ZoltanPartitioner.h */

#ifndef AMDIS_ZOLTAN_PARTITIONER_H
#define AMDIS_ZOLTAN_PARTITIONER_H

#ifdef HAVE_ZOLTAN

#include <map>
#include <set>
#include <string>

#include <zoltan_cpp.h>
#include "AMDiS_fwd.h"
#include "parallel/MeshPartitioner.h"


namespace AMDiS { namespace Parallel {

  class ZoltanPartitioner : public MeshPartitioner
  {
  public:
    ZoltanPartitioner(std::string name, MPI::Intracomm *comm);

    ~ZoltanPartitioner() {}
      
    /// \ref MeshPartitioner::partition
    bool partition(std::map<int, double> &weights, PartitionMode mode = INITIAL);

    void createPartitionMap(std::map<int, int> &pMap)
    {
      pMap = partitionMap;
    }

  protected:
    void createPartitionMap();

  protected:
    /// Data structure used by the Zoltan library.
    Zoltan zoltan;

    /// Pointer to a map that weights all macro elements in rank's partition.
    std::map<int, double> *elWeights;

    /// parmeter map to configure zoltan
    std::map<std::string,std::string> paramMapInitial;
    std::map<std::string,std::string> paramMap;

    friend class ZoltanFunctions;
  };


  class ZoltanFunctions
  {
  public:
    // === Zoltan query function for "normal" element partitioning. ===

    static int getNumObj(void *data, int *ierr);   

    static void getObjectList(void *data, 
			      int sizeGid, 
			      int sizeLid,
			      ZOLTAN_ID_PTR globalId, 
			      ZOLTAN_ID_PTR localId,
			      int wgt_dim, 
			      float *obj_wgts, 
			      int *ierr);

    static int getNumGeom(void *data, int *ierr);

    static void getGeomMulti(void *data, 
			     int num_gid_entries, 
			     int num_lid_entries, 
			     int num_obj, 
			     ZOLTAN_ID_PTR global_ids, 
			     ZOLTAN_ID_PTR local_ids, 
			     int num_dim, 
			     double *geom_vec, 
			     int *ierr); 

    static void getNumEdgeMulti(void *data, 
				int num_gid_entries, 
				int num_lid_entries, 
				int num_obj, 
				ZOLTAN_ID_PTR global_ids, 
				ZOLTAN_ID_PTR local_ids, 
				int *num_edges, 
				int *ierr);

    static void getEdgeListMulti(void *data, 
				 int num_gid_entries, 
				 int num_lid_entries, 
				 int num_obj, 
				 ZOLTAN_ID_PTR global_ids, 
				 ZOLTAN_ID_PTR local_ids, 
				 int *num_edges, 
				 ZOLTAN_ID_PTR nbor_global_id, 
				 int *nbor_procs, 
				 int wgt_dim, 
				 float *ewgts, 
				 int *ierr);


    // === Zoltan query functions for box partitioning. ===

    static int box_getNumObj(void *data, int *ierr);   

    static void box_getObjectList(void *data, 
				  int sizeGid, 
				  int sizeLid,
				  ZOLTAN_ID_PTR globalId, 
				  ZOLTAN_ID_PTR localId,
				  int wgt_dim, 
				  float *obj_wgts, 
				  int *ierr);

    static int box_getNumGeom(void *data, int *ierr);

    static void box_getGeomMulti(void *data, 
				 int num_gid_entries, 
				 int num_lid_entries, 
				 int num_obj, 
				 ZOLTAN_ID_PTR global_ids, 
				 ZOLTAN_ID_PTR local_ids, 
				 int num_dim, 
				 double *geom_vec, 
				 int *ierr); 

    static void box_getNumEdgeMulti(void *data, 
				    int num_gid_entries, 
				    int num_lid_entries, 
				    int num_obj, 
				    ZOLTAN_ID_PTR global_ids, 
				    ZOLTAN_ID_PTR local_ids, 
				    int *num_edges, 
				    int *ierr);

    static void box_getEdgeListMulti(void *data, 
				     int num_gid_entries, 
				     int num_lid_entries, 
				     int num_obj, 
				     ZOLTAN_ID_PTR global_ids, 
				     ZOLTAN_ID_PTR local_ids, 
				     int *num_edges, 
				     ZOLTAN_ID_PTR nbor_global_id, 
				     int *nbor_procs, 
				     int wgt_dim, 
				     float *ewgts, 
				     int *ierr);

  };

} }

#endif // HAVE_ZOLTAN

#endif
