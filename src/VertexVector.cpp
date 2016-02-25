#include "VertexVector.hpp"

#include <utility> // std::swap

#include "DOFAdmin.hpp"
#include "DOFIterator.hpp"

namespace AMDiS
{
  VertexVector::VertexVector(DOFAdmin const* a, std::string n)
    : DOFVectorDOF(),
      admin(a)
  {
    name = n;
    feSpace = NULL;
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


  void VertexVector::compressDofContainer(int size, std::vector<DegreeOfFreedom>& newDOF)
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

} // end namespace AMDiS
