#include "SystemVector.h"

namespace AMDiS 
{

  int SystemVector::getUsedSize() const 
  {
    int totalSize = 0;
    for (size_t i = 0; i < vectors.size(); i++)
      totalSize += vectors[i]->getUsedSize();
    return totalSize;
  }

  
  double& SystemVector::operator[](DegreeOfFreedom index) 
  {
    DegreeOfFreedom localIndex = index;
    DegreeOfFreedom vectorIndex = 0;

    while (localIndex >= vectors[vectorIndex]->getUsedSize())
      localIndex -= vectors[vectorIndex++]->getUsedSize();

    return (*(vectors[vectorIndex]))[localIndex];
  }

  
  double SystemVector::operator[](DegreeOfFreedom index) const 
  {
    DegreeOfFreedom localIndex = index;
    DegreeOfFreedom vectorIndex = 0;

    while (localIndex >= vectors[vectorIndex]->getUsedSize())
      localIndex -= vectors[vectorIndex++]->getUsedSize();

    return (*(vectors[vectorIndex]))[localIndex];
  }

  
  void SystemVector::set(double value) 
  {
    for (size_t i = 0; i < vectors.size(); i++)
      vectors[i]->set(value);
  }

  
  void SystemVector::setCoarsenOperation(RefineCoarsenOperation op) 
  { 
    for (size_t i = 0; i < vectors.size(); i++)
      vectors[i]->setCoarsenOperation(op); 
  }
  

  void SystemVector::setRefineOperation(RefineCoarsenOperation op) 
  { 
    for (size_t i = 0; i < vectors.size(); i++)
      vectors[i]->setRefineOperation(op); 
  }

  
  SystemVector& SystemVector::operator=(double value) 
  {
    for (size_t i = 0; i < vectors.size(); i++)
      (*(vectors[i])) = value;
    return *this;
  }

  
  SystemVector& SystemVector::operator=(const SystemVector& rhs) 
  {
    TEST_EXIT_DBG(rhs.vectors.size() == vectors.size())("Invalied sizes!\n");
    for (size_t i = 0; i < vectors.size(); i++)
      (*(vectors[i])) = (*(rhs.getDOFVector(i)));

    return *this;
  }

  
  void SystemVector::copy(const SystemVector& rhs) 
  {
    TEST_EXIT_DBG(getSize() == rhs.getSize())("Invalid sizes!\n");
    for (size_t i = 0; i < vectors.size(); i++)
      vectors[i]->copy(*(const_cast<SystemVector&>(rhs).getDOFVector(i)));
  }

  
  void SystemVector::interpol(std::vector<std::function<double(WorldVector<double>)> >& f) 
  {
    for (size_t i = 0; i < vectors.size(); i++)
      vectors[i]->interpol(f[i]);
  }

  
  void SystemVector::interpol(SystemVector *v, double factor) 
  {
    for (int i = 0; i < v->getSize(); i++)
      vectors[i]->interpol(v->getDOFVector(i), factor);
  }

  
  int SystemVector::calcMemoryUsage() const
  {
    int result = 0;
    for (size_t i = 0; i < vectors.size(); i++)
      result += vectors[i]->calcMemoryUsage();
    result += sizeof(SystemVector);

    return result;
  }
  
} // end namespace AMDiS
