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



/** \file SimplePartitioner.h */

#ifndef AMDIS_SIMPLE_PARTITIONER_H
#define AMDIS_SIMPLE_PARTITIONER_H

#include "AMDiS_fwd.h"
#include "Global.h"
#include "parallel/MeshPartitioner.h"

namespace AMDiS { namespace Parallel {

  /**
   * The "Simple partitioner" does not change the initial partitioning which is more
   * or less a random assignment of elements to ranks. This partitioner may be useful
   * for either debugging purposes or if the number of macro elements is equal to the
   * number of processes. In this case, neither ParMetis nor Zoltan will be able to
   * compute a valid partition. But a random one-to-one partition is the best possible
   * in this case.
   */
  class SimplePartitioner : public MeshPartitioner
  {
  public:
    SimplePartitioner(std::string name, MPI::Intracomm *comm)
      : MeshPartitioner(name, comm)
    {}

    ~SimplePartitioner() {}

    /// \ref MeshPartitioner::partition
    bool partition(std::map<int, double> &elemWeights, PartitionMode mode = INITIAL)
    {
      return true;
    }

    void createPartitionMap(std::map<int, int>& pMap)
    {
      pMap = partitionMap;
    }
  };
} }

#endif
