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



/** \file PeriodicBC.h */

#ifndef AMDIS_PERIODICBC_H
#define AMDIS_PERIODICBC_H

#include <map>
#include <vector>
#include "AMDiS_fwd.h"
#include "BoundaryCondition.h"
#include "FixVec.h"

namespace AMDiS {

  template<typename T>
  class DimVecLess 
  {
  public:
    bool operator()(const DimVec<T> &v1, const DimVec<T> &v2) const
    {
      int  size = v1.getSize();
      for (int i = 0; i < size; i++) {
	if (v1[i] < v2[i]) 
	  return true;
	if (v1[i] > v2[i]) 
	  return false;
      }
      return false;
    }
  };

  class PeriodicDOFMapping
  {
  public:
    PeriodicDOFMapping(const BasisFunction *basFcts);

    ~PeriodicDOFMapping();

  public:
    static PeriodicDOFMapping* providePeriodicDOFMapping(const BasisFunction *basFcts);

    const DegreeOfFreedom *getDOFPermutation(const DimVec<int> &vertexPermutation);

  protected:
    /// Basis functions the mapping object is defined on.
    const BasisFunction *basFcts;

    std::map<DimVec<int>, DegreeOfFreedom*, DimVecLess<int> > dofPermutation;

    /// Maps from coordinates (in barycentric coordinates) of a local basis
    /// function to the local index of that basis function.
    std::map<DimVec<double>, int, DimVecLess<double> > indexOfCoords;

    /// Global array that provids for each existing basis functions a unique
    /// mapping object.
    static std::vector<PeriodicDOFMapping*> mappings;
  };


  /**
   * \ingroup Assembler
   *
   * \brief
   * Sub class of BoundaryCondition. Implements Periodic boundary conditions.
   */
  class PeriodicBC : public BoundaryCondition
  {
  public:
    /// Constructor.
    PeriodicBC(BoundaryType type, const FiniteElemSpace *rowFeSpace);

    ~PeriodicBC();

    void initMatrix(DOFMatrix* matrix);

    void fillBoundaryCondition(DOFMatrix *matrix,
			       ElInfo *elInfo,
			       const DegreeOfFreedom *dofIndices,
			       const BoundaryType *localBound,
			       int nBasFcts);

    void exitMatrix(DOFMatrix* matrix);

    void exitVector(DOFVectorBase<double>* vector);

    /// We are defining periodic boundary conditions, so return always true here.
    bool isPeriodic()
    {
      return true;
    }

  protected:
    VertexVector *associated;

    PeriodicDOFMapping *periodicDOFMapping;

    DOFMatrix *masterMatrix;
  };

}

#endif
