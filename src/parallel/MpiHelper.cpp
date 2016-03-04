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


#include <mpi.h>
#include "MpiHelper.h"

namespace AMDiS
{
  namespace Parallel
  {

    namespace mpi
    {

      void globalAdd(MPI::Intracomm& mpiComm, double& value)
      {
        double valCopy = value;
        mpiComm.Allreduce(&valCopy, &value, 1, MPI_DOUBLE, MPI_SUM);
      }

      void globalAdd(MPI::Intracomm& mpiComm, int& value)
      {
        int valCopy = value;
        mpiComm.Allreduce(&valCopy, &value, 1, MPI_INT, MPI_SUM);
      }

      void globalMin(double& value)
      {
        double valCopy = value;
        MPI::COMM_WORLD.Allreduce(&valCopy, &value, 1, MPI_DOUBLE, MPI_MIN);
      }

      void globalMin(int& value)
      {
        int valCopy = value;
        MPI::COMM_WORLD.Allreduce(&valCopy, &value, 1, MPI_INT, MPI_MIN);
      }

      void globalMax(double& value)
      {
        double valCopy = value;
        MPI::COMM_WORLD.Allreduce(&valCopy, &value, 1, MPI_DOUBLE, MPI_MAX);
      }

      void globalMax(int& value)
      {
        int valCopy = value;
        MPI::COMM_WORLD.Allreduce(&valCopy, &value, 1, MPI_INT, MPI_MAX);
      }

      void startRand()
      {
        srand(time(0) * (MPI::COMM_WORLD.Get_rank() + 1));
      }
    }
  }
}

