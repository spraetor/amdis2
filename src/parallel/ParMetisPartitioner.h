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



/** \file ParMetisPartitioner.h */

#ifndef AMDIS_PARMETIS_PARTITIONER_H
#define AMDIS_PARMETIS_PARTITIONER_H

#include <map>
#include <set>
#include <parmetis.h>
#include <mpi.h>

#include "AMDiS_fwd.h"
#include "Global.h"
#include "parallel/MeshPartitioner.h"
#include "parallel/ParallelDofMapping.h"

namespace AMDiS { namespace Parallel {

  class ParMetisGraph;

  class ParMetisMesh
  {
  public:
    ParMetisMesh(Mesh *mesh, MPI::Intracomm *comm, 
		 std::map<int, bool>& elementInRank,
		 DofMap *mapLocalGlobal);

    ~ParMetisMesh();

    inline void setParMetisIndex(int amdisIndex, int parMetisIndex) 
    {
      elem_a2p[amdisIndex] = parMetisIndex + 1;
    }

    inline int getParMetisIndex(int amdisIndex) 
    {
      int result = elem_a2p[amdisIndex];
      TEST_EXIT(result > 0)("invalid index\n");
      return result - 1;
    }

    inline void setAMDiSIndex(int parMetisIndex, int amdisIndex) 
    {
      elem_p2a[parMetisIndex] = amdisIndex;
    }

    inline int getAMDiSIndex(int parMetisIndex) 
    {
      return elem_p2a[parMetisIndex];
    }

    inline int *getAMDiSIndices() 
    {
      return elem_p2a;
    }

    inline int *getElementPtr() 
    { 
      return eptr; 
    }

    inline int *getElementInd() 
    { 
      return eind; 
    }

    inline int *getElementDist() 
    { 
      return elmdist; 
    }

    inline int getDim() 
    { 
      return dim; 
    }

    inline float *getXYZ() 
    { 
      return xyz; 
    }

    inline int getNumElements() 
    { 
      return nElements; 
    }

  protected:
    int *eptr;

    int *eind;

    /* \brief
     * Array that specifies the distribution of the mesh elements.
     *
     * elmdist[0] = 0;
     * elmdist[1] = number of elements of proc 0;
     * elmdist[2] = elmdist[1] + number of elements of proc 1;
     *    ...
     */
    int *elmdist;

    int dim;

    float *xyz;

    int nElements;

    std::map<int, int> elem_a2p;

    int *elem_p2a;

    /// The MPI communicator that should be used for mesh partition.
    MPI::Intracomm *mpiComm;
  };


  class ParMetisGraph
  {
  public:
    ParMetisGraph(ParMetisMesh *parMetisMesh,
		  MPI::Intracomm *comm,
		  int ncommonnodes = -1);

    ~ParMetisGraph();

    inline int *getXAdj() 
    { 
      return xadj; 
    }

    inline int *getAdjncy() 
    { 
      return adjncy; 
    }

    void print();

  protected:
    ParMetisMesh *parMetisMesh;

    int *xadj;

    int *adjncy;
  };


  class ParMetisPartitioner : public MeshPartitioner
  {
  public:
    ParMetisPartitioner(std::string name, MPI::Intracomm *comm)
      : MeshPartitioner(name, comm),
        parMetisMesh(NULL),
	itr(1000000.0)
    {}

    ~ParMetisPartitioner();

    /// \ref MeshPartitioner::partition
    bool partition(std::map<int, double> &elemWeights,
		   PartitionMode mode = INITIAL);

    void createPartitionMap(std::map<int, int>& partitionMap);

    void setItr(double value)
    {
      itr = value;
    }

  protected:
    // Returns true, if the mesh partitioning could be distributed. If this is not the
    // case, i.e., because there are empty partitions, the function returns false.
    bool distributePartitioning(int *part);

  protected:
    ParMetisMesh *parMetisMesh;

    double itr;
  };
} }

#endif
