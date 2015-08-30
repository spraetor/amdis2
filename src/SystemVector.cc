#include "SystemVector.h"
#include "MatrixVector.h"
#include "DOFVector.h"
#include "DOFVectorOperations.h"

namespace AMDiS
{
  SystemVector::SystemVector(std::string name_,
                             std::vector<const FiniteElemSpace*> feSpace_,
                             int size,
                             bool createVec_)
    : name(name_),
      componentSpaces(feSpace_),
      vectors(size),
      createVec(createVec_)
  {
    if (createVec_)
      for (int i = 0; i < size; i++)
        vectors[i] = new DOFVector<double>(componentSpaces[i], "tmp");
  }


  SystemVector::SystemVector(const SystemVector& rhs)
    : name(rhs.getName()),
      componentSpaces(rhs.getFeSpaces()),
      vectors(rhs.getSize())
  {
    for (size_t i = 0; i < vectors.size(); i++)
      vectors[i] = new DOFVector<double>(*rhs.getDOFVector(i));
  }


  SystemVector::~SystemVector()
  {
    if (createVec)
    {
      for (size_t i = 0; i < vectors.size(); i++)
        delete vectors[i];
    }
  }


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


  SystemVector& SystemVector::operator=(SystemVector const& rhs)
  {
    TEST_EXIT_DBG(rhs.vectors.size() == vectors.size())("Invalied sizes!\n");
    for (size_t i = 0; i < vectors.size(); i++)
      (*(vectors[i])) = (*(rhs.getDOFVector(i)));

    return *this;
  }


  void SystemVector::copy(SystemVector const& rhs)
  {
    TEST_EXIT_DBG(getSize() == rhs.getSize())("Invalid sizes!\n");
    for (size_t i = 0; i < vectors.size(); i++)
      vectors[i]->copy(*(const_cast<SystemVector&>(rhs).getDOFVector(i)));
  }


  void SystemVector::interpol(std::vector<std::function<double(WorldVector<double>)>>& f)
  {
    for (size_t i = 0; i < vectors.size(); i++)
      vectors[i]->interpol(f[i]);
  }


  void SystemVector::interpol(SystemVector* v, double factor)
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


  /* ----- OPERATORS WITH SYSTEM-VECTORS ------------------------------------ */


  SystemVector& operator*=(SystemVector& x, double d)
  {
    for (int i = 0; i < x.getSize(); i++)
      *(x.getDOFVector(i)) *= d;
    return x;
  }


  SystemVector operator*(SystemVector x, double d)
  {
    for (int i = 0; i < x.getSize(); i++)
      (*(x.getDOFVector(i))) *= d;
    return x;
  }


  SystemVector operator*(double d, SystemVector x)
  {
    for (int i = 0; i < x.getSize(); i++)
      (*(x.getDOFVector(i))) *= d;
    return x;
  }


  double operator*(SystemVector const& x, SystemVector const& y)
  {
    TEST_EXIT_DBG(x.getSize() == y.getSize())("invalid size\n");
    double result = 0.0;
    for (int i = 0; i < x.getSize(); i++)
      result += (*x.getDOFVector(i)) * (*y.getDOFVector(i));
    return result;
  }


  SystemVector& operator+=(SystemVector& x, SystemVector const& y)
  {
    TEST_EXIT_DBG(x.getSize() == y.getSize())("invalid size\n");
    for (int i = 0; i < x.getSize(); i++)
      (*(x.getDOFVector(i))) += (*(y.getDOFVector(i)));
    return x;
  }


  SystemVector& operator-=(SystemVector& x, SystemVector const& y)
  {
    TEST_EXIT_DBG(x.getSize() == y.getSize())("invalid size\n");
    for (int i = 0; i < x.getSize(); i++)
      (*(x.getDOFVector(i))) -= (*(y.getDOFVector(i)));
    return x;
  }


  SystemVector operator+(SystemVector x, SystemVector const& y)
  {
    TEST_EXIT_DBG(x.getSize() == y.getSize())("invalid size\n");
    for (int i = 0; i < x.getSize(); i++)
      (*(x.getDOFVector(i))) += (*(y.getDOFVector(i)));
    return x;
  }


  double norm(SystemVector const* x)
  {
    double result = 0.0;
    for (int i = 0; i < x->getSize(); i++)
      result += x->getDOFVector(i)->squareNrm2();
    return std::sqrt(result);
  }


  double L2Norm(SystemVector const* x)
  {
    double result = 0.0;
    for (int i = 0; i < x->getSize(); i++)
      result += x->getDOFVector(i)->L2NormSquare();
    return std::sqrt(result);
  }


  double H1Norm(SystemVector const* x)
  {
    double result = 0.0;
    for (int i = 0; i < x->getSize(); i++)
      result += x->getDOFVector(i)->H1NormSquare();
    return std::sqrt(result);
  }

} // end namespace AMDiS
