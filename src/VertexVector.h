/** \file VertexVector.h */

#pragma once

#include "DOFVector.h"

namespace AMDiS 
{

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
