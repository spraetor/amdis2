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



/** \file VertexVector.h */

#ifndef AMDIS_VERTEXVECTOR_H
#define AMDIS_VERTEXVECTOR_H

#include "DOFVector.h"

namespace AMDiS {

  class VertexVector : public DOFVectorDOF
  {
  public:
    struct Iterator : public DOFIterator<DegreeOfFreedom> 
    {
      Iterator(VertexVector *c, DOFIteratorType type)
	: DOFIterator<DegreeOfFreedom>(const_cast<DOFAdmin*>(c->getAdmin()), 
				       dynamic_cast<DOFIndexed<DegreeOfFreedom>*>(c), 
				       type)
      { }   
    };

    /// Constructor
    VertexVector(const DOFAdmin *admin, std::string name);

    /// Destructor, calls \ref removeDOFIndexed and \ref removeDOFContainer on \ref admin.
    ~VertexVector();

    void freeDofIndex(DegreeOfFreedom dof)
    {
      FUNCNAME_DBG("VertexVector::freeDofIndex()");
      TEST_EXIT_DBG(dof < static_cast<int>(vec.size()))("Should not happen!\n");

      vec[dof] = dof;
    }

    const DOFAdmin *getAdmin() const
    { 
      return admin; 
    }

    void resize(int size);

    void set(DegreeOfFreedom val);

    void compressDofContainer(int size, std::vector<DegreeOfFreedom> &newDOF);

    void changeDofIndices(std::map<DegreeOfFreedom, DegreeOfFreedom>& dofIndexMap);

  protected:
    const DOFAdmin *admin;
  };
  
} // end namespace AMDiS

#endif
