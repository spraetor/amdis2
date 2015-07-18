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



/** \file SystemVector.h */

#ifndef AMDIS_SYSTEMVECTOR_H
#define AMDIS_SYSTEMVECTOR_H

#include "MatrixVector.h"
#include "DOFVector.h"
#include "CreatorInterface.h"

namespace AMDiS {

  /// A system vector is a vector of dof vectors used for vector valued problems.
  class SystemVector
  {
  public:
    /// Constructor.
    SystemVector(std::string name_,
		 std::vector<const FiniteElemSpace*> feSpace_, 
		 int size,
		 bool createVec_ = false)
      : name(name_),
	componentSpaces(feSpace_),
	vectors(size),
	createVec(createVec_)
    {
      if (createVec_)
	for (int i = 0; i < size; i++)
	  vectors[i] = new DOFVector<double>(componentSpaces[i], "tmp");
    }

    /// Copy Constructor.
    SystemVector(const SystemVector& rhs)
      : name(rhs.getName()),
	componentSpaces(rhs.getFeSpaces()),
	vectors(rhs.getSize())
    {
      for (size_t i = 0; i < vectors.size(); i++)
	vectors[i] = new DOFVector<double>(*rhs.getDOFVector(i));
    }

    /// Destructor, deletes all DOFVectors
    ~SystemVector() 
    {
      if (createVec) {
	for (size_t i = 0; i < vectors.size(); i++)
	  delete vectors[i];
      }
    }

    /// Sets \ref vectors[index] = vec.
    void setDOFVector(int index, DOFVector<double> *vec) 
    {
      TEST_EXIT_DBG(index < getSize())
	("Invalid index %d!\n", index);
      vectors[index] = vec;
    }

    /// Returns \ref vectors[index].
    DOFVector<double> *getDOFVector(int index) 
    {
      TEST_EXIT_DBG(index < getSize())
	("Invalid index %d!\n", index);
      return vectors[index];
    }

    /// Returns \ref vectors[index].
    const DOFVector<double> *getDOFVector(int index) const 
    {
      TEST_EXIT_DBG(index < getSize())
	("Invalid index %d!\n", index);
      return vectors[index];
    }

    std::string getName() const
    {
      return name;
    }

    /// Returns sum of used vector sizes.
    int getUsedSize() const;
    
    /// Returns number of contained vectors.
    int getSize() const 
    {
      return static_cast<int>(vectors.size());
    }

    /// Returns the fe space for a given component.
    const FiniteElemSpace *getFeSpace(int i) const 
    { 
      return componentSpaces[i]; 
    }

    /// Returns the fe spaces for all components.
    std::vector<const FiniteElemSpace*> getFeSpaces() const 
    {
      return componentSpaces;
    }

    /// Here the system vector is interpreted as one large vector. The given
    /// is used as a global index which indicates a local vector number and
    /// a local index on this vector. The operator returns this local vector
    /// at the local index.
    double& operator[](int index);

    /// For const access.
    double operator[](int index) const;

    /// Sets all entries in all vectors to value.
    void set(double value);
    
    /// Sets all entries in all vectors to value.
    SystemVector& operator=(double value);

    /// Copy assignment function
    void copy(const SystemVector& rhs);
    
    /// Copy-Assignement operator.
    SystemVector& operator=(const SystemVector& rhs);

    /// Set the coarsen operation for all DOFVectors
    void setCoarsenOperation(RefineCoarsenOperation op);

    /// Set the refine operation for all DOFVectors
    void setRefineOperation(RefineCoarsenOperation op);

    void interpol(std::vector<std::function<double(WorldVector<double>)> > *f);

    void interpol(SystemVector *v, double factor);

    int calcMemoryUsage() const;

  protected:
    /// Name of the system vector
    std::string name;

    /// Finite element space.
    std::vector<const FiniteElemSpace*> componentSpaces;

    /// Local dof vectors.
    std::vector<DOFVector<double>*> vectors;
    
    bool createVec;
  };


  /// multiplication with scalar
  inline const SystemVector& operator*=(SystemVector& x, double d) 
  {
    for (int i = 0; i < x.getSize(); i++)
      *(x.getDOFVector(i)) *= d;
    return x;
  }

  /// scalar product
  inline double operator*(SystemVector& x, SystemVector& y) 
  {
    TEST_EXIT_DBG(x.getSize() == y.getSize())("invalid size\n");
    double result = 0.0;
    for (int i = 0; i < x.getSize(); i++)
      result += (*x.getDOFVector(i)) * (*y.getDOFVector(i));
    return result;
  }

  /// addition of two system vectors
  inline const SystemVector& operator+=(SystemVector& x, const SystemVector& y) 
  {
    TEST_EXIT_DBG(x.getSize() == y.getSize())("invalid size\n");
    for (int i = 0; i < x.getSize(); i++)
      (*(x.getDOFVector(i))) += (*(y.getDOFVector(i)));
    return x;
  }

  /// subtraction of two system vectors.
  inline const SystemVector& operator-=(SystemVector& x, SystemVector& y) 
  {
    TEST_EXIT_DBG(x.getSize() == y.getSize())("invalid size\n");
    for (int i = 0; i < x.getSize(); i++)
      (*(x.getDOFVector(i))) -= (*(y.getDOFVector(i)));
    return x;
  }

  /// multiplication with a scalar
  inline SystemVector operator*(SystemVector& x, double d) 
  {
    SystemVector result = x;
    for (int i = 0; i < x.getSize(); i++)
      (*(result.getDOFVector(i))) *= d;
    return result;
  }

  /// multiplication with a scalar
  inline SystemVector operator*(double d, SystemVector& x) 
  {
    SystemVector result = x;
    for (int i = 0; i < x.getSize(); i++)
      (*(result.getDOFVector(i))) *= d;
    return result;
  }

  /// addition of two system vectors
  inline SystemVector operator+(const SystemVector& x, const SystemVector& y)
  {
    TEST_EXIT_DBG(x.getSize() == y.getSize())("invalid size\n");
    SystemVector result = x;
    for (int i = 0; i < x.getSize(); i++)
      (*(result.getDOFVector(i))) += (*(y.getDOFVector(i)));
    return result;
  }

  /// Calls SystemVector::set(). Used for solving.
  inline void set(SystemVector& x, double value) 
  {
    x.set(value);
  } 

  /// Calls SystemVector::set(). Used for solving.
  inline void setValue(SystemVector& x, double value) 
  {
    x.set(value);
  }

  /// Norm of system vector.
  inline double norm(SystemVector* x) 
  {
    double result = 0.0;
    for (int i = 0; i < x.getSize(); i++)
      result += x->getDOFVector(i)->squareNrm2();
    return std::sqrt(result);
  }

  /// L2 norm of system vector.
  inline double L2Norm(SystemVector* x) 
  {
    double result = 0.0;
    for (int i = 0; i < x.getSize(); i++)
      result += x->getDOFVector(i)->L2NormSquare();
    return std::sqrt(result);
  }

  /// H1 norm of system vector.
  inline double H1Norm(SystemVector* x) 
  {
    double result = 0.0;
    for (int i = 0; i < x.getSize(); i++)
      result += x->getDOFVector(i)->H1NormSquare();
    return std::sqrt(result);
  }

  /// Returns SystemVector::getUsedSize().
  inline int size(SystemVector* vec) 
  {
    return vec->getUsedSize();
  }

} // end namespace AMDiS

#endif // AMDIS_SYSTEMVECTOR_H
