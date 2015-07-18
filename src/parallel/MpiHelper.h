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



/** \file MpiHelper.h */

#ifndef AMDIS_MPIHELPER_H
#define AMDIS_MPIHELPER_H

#include <mpi.h>
#include <time.h>
#include <stdlib.h>
#include "Global.h"
#include "traits/category.hpp"

namespace AMDiS { namespace Parallel {

  namespace mpi 
  {
  
  // TODO: replace globalAdd, global* by boost::mpi::all_reduce using functors
  
    void globalAdd(MPI::Intracomm &mpiComm, double &value);

    inline void globalAdd(double &value)
    {
      globalAdd(MPI::COMM_WORLD, value);
    }
    
    void globalAdd(MPI::Intracomm &mpiComm, int &value);

    inline void globalAdd(int &value)
    {
      globalAdd(MPI::COMM_WORLD, value);
    }

    template<typename VectorType>
    void globalAdd_dispatch(VectorType &value, tag::vector)
    {
      typedef typename traits::category<VectorType>::size_type size_type;
      for (size_type i = 0; i < num_rows(value); i++)
	globalAdd(value[i]);
    }
    
    template<typename T, typename Tag>
    void globalAdd_dispatch(T &value, Tag t)
    {
      WARNING("Unknown type for globalAdd. Can not sum up the values of all processors!\n");
    }
    
    template<typename T>
    void globalAdd(T &value) 
    {
      globalAdd_dispatch(value, typename traits::category<T>::tag());
    }

    void globalMin(double &value);

    void globalMin(int &value);

    template<typename T>
    void globalMin(T &value) 
    {
      WARNING("Unknown type for globalMin. Can not determine minimal value of all processors!\n");
    }
    
    void globalMax(double &value);

    void globalMax(int &value);

    template<typename T>
    void globalMax(T &value) 
    {
      WARNING("Unknown type for globalMax. Can not determine maximal value of all processors!\n");
    }
    
    void startRand();

    /** \brief
     * In many situations a rank computes a number of local DOFs. Then all
     * ranks want to know the number of global DOFs and the starting 
     * displacment number of the DOF numbering in each rank.
     *
     * \param[in]   mpiComm        The MPI communicator.
     * \param[in]   nRankDofs      The number of local DOFs.
     * \param[out]  rStartDofs     Displacment of the DOF numbering. On rank n
     *                             this is the sum of all local DOF numbers in
     *                             ranks 0 to n - 1.
     * \param[out]  nOverallDofs   Global sum of nRankDofs. Is equal on all
     *                             ranks.
     */
    inline void getDofNumbering(MPI::Intracomm& mpiComm,
				int nRankDofs, 
				int& rStartDofs, 
				int& nOverallDofs)
    {
      rStartDofs = 0;
      nOverallDofs = 0;
      mpiComm.Scan(&nRankDofs, &rStartDofs, 1, MPI_INT, MPI_SUM);
      rStartDofs -= nRankDofs;
      mpiComm.Allreduce(&nRankDofs, &nOverallDofs, 1, MPI_INT, MPI_SUM);
    }
  }

} }

#endif
