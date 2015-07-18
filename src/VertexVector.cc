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

#include <numeric>

#include "VertexVector.h"
#include "DOFAdmin.h"
#include "DOFIterator.h"

namespace AMDiS {

  VertexVector::VertexVector(const DOFAdmin *a, std::string n)
    : DOFVectorDOF(),
      name(n),
      feSpace(NULL),
      admin(a)
  {
    const_cast<DOFAdmin*>(admin)->addDOFIndexed(this);
    const_cast<DOFAdmin*>(admin)->addDOFContainer(this);
  }


  VertexVector::~VertexVector()
  {
    const_cast<DOFAdmin*>(admin)->removeDOFIndexed(this);
    const_cast<DOFAdmin*>(admin)->removeDOFContainer(this);
  }


  void VertexVector::set(DegreeOfFreedom val)
  {
    DOFIteratorBase it(const_cast<DOFAdmin*>(admin), USED_DOFS);
    for (it.reset(); !it.end(); ++it)
      if (!it.isDofFree()) 
	operator[](it.getDOFIndex()) = val;
  } 

  
  void VertexVector::resize(int size) 
  {
    int oldSize = static_cast<int>(vec.size());
    vec.resize(size);
    for (int i = oldSize; i < size; i++)
      vec[i] = i;
  }
  

  void VertexVector::compressDofContainer(int size, std::vector<DegreeOfFreedom> &newDOF) 
  {
    DOFContainer::compressDofContainer(size, newDOF);
    int totalSize = getAdmin()->getSize();
    for (int i = size; i < totalSize; i++)
      vec[i] = i;
  }


  void VertexVector::changeDofIndices(std::map<DegreeOfFreedom, DegreeOfFreedom>& dofIndexMap)
  {
    std::vector<DegreeOfFreedom> tmp(vec.size());
    std::iota(tmp.begin(), tmp.end(), 0);
    
    for (auto const& d : dofIndexMap)
      if (vec[d.first] == -1)
	tmp[d.second] = -1;
      else 
	tmp[d.second] = dofIndexMap[vec[d.first]];

    std::swap(vec, tmp);
  }

}
