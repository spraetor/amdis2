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



/** \file CheckerPartitioner.h */

#ifndef AMDIS_CHECKER_PARTITIONER_H
#define AMDIS_CHECKER_PARTITIONER_H

#include "AMDiS_fwd.h"
#include "Global.h"
#include "Initfile.h"
#include "parallel/MeshPartitioner.h"

namespace AMDiS { namespace Parallel {

  class CheckerPartitioner : public MeshPartitioner
  {
  public:
    CheckerPartitioner(std::string name, MPI::Intracomm *comm);

    ~CheckerPartitioner() {}

    bool createInitialPartitioning();

    void createTetrahedronStripes();

    /// \ref MeshPartitioner::partition
    bool partition(std::map<int, double> &elemWeights, PartitionMode mode = INITIAL)
    {
      return true;
    }

    void createPartitionMap(std::map<int, int>& pMap)
    {
      pMap = partitionMap;
    }

  protected:
    int mpiRank, mpiSize;

    /// 0: standard mode: each node gets one box
    /// 1: x-stripes: each node gets one x-stripe of boxes
    /// 2: y-stripes: each node gets one y-stripe of boxes
    /// 3: z-stripes: each node gets one y-stripe of boxes
    /// 4: tetrahedron-stripes: alias Hieram mode :)
    int mode;

    /// Only used in mode 4.
    std::map<int, int> elStripeInRank;

    bool multilevel;
  };
} }

#endif
