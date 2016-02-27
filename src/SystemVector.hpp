#pragma once

#include <string>
#include <vector>

#include "AMDiS_fwd.hpp"
#include "AMDiS_base.hpp"
#include "Log.hpp"
#include "MatrixVector_fwd.hpp"

namespace AMDiS
{

  /// A system vector is a vector of dof vectors used for vector valued problems.
  class SystemVector
  {
  public:
    /// Constructor.
    SystemVector(std::string name_,
                 std::vector<FiniteElemSpace const*> feSpace_,
                 int size,
                 bool createVec_ = false);

    /// Copy Constructor.
    SystemVector(SystemVector const& rhs);

    /// Destructor, deletes all DOFVectors
    ~SystemVector();

    /// Sets \ref vectors[index] = vec.
    void setDOFVector(int index, DOFVector<double>* vec)
    {
      TEST_EXIT_DBG(index < getSize())
      ("Invalid index %d!\n", index);
      vectors[index] = vec;
    }

    /// Returns \ref vectors[index].
    DOFVector<double>* getDOFVector(int index)
    {
      TEST_EXIT_DBG(index < getSize())
      ("Invalid index %d!\n", index);
      return vectors[index];
    }

    /// Returns \ref vectors[index].
    DOFVector<double> const* getDOFVector(int index) const
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
    FiniteElemSpace const* getFeSpace(int i) const
    {
      return componentSpaces[i];
    }

    /// Returns the fe spaces for all components.
    std::vector<FiniteElemSpace const*> getFeSpaces() const
    {
      return componentSpaces;
    }

    /// Here the system vector is interpreted as one large vector. The given
    /// is used as a global index which indicates a local vector number and
    /// a local index on this vector. The operator returns this local vector
    /// at the local index.
    double& operator[](DegreeOfFreedom index);

    /// For const access.
    double operator[](DegreeOfFreedom index) const;

    /// Sets all entries in all vectors to value.
    void set(double value);

    /// Sets all entries in all vectors to value.
    SystemVector& operator=(double value);

    /// Copy assignment function
    void copy(SystemVector const& rhs);

    /// Copy-Assignement operator.
    SystemVector& operator=(SystemVector const& rhs);

    /// Set the coarsen operation for all DOFVectors
    void setCoarsenOperation(RefineCoarsenOperation op);

    /// Set the refine operation for all DOFVectors
    void setRefineOperation(RefineCoarsenOperation op);

    void interpol(std::vector<std::function<double(WorldVector<double>)>>& f);

    void interpol(SystemVector* v, double factor);

    size_t calcMemoryUsage() const;

  protected:
    /// Name of the system vector
    std::string name;

    /// Finite element space.
    std::vector<FiniteElemSpace const*> componentSpaces;

    /// Local dof vectors.
    std::vector<DOFVector<double>*> vectors;

    bool createVec;
  };


  /* ----- OPERATORS WITH SYSTEM-VECTORS ------------------------------------ */


  /// multiplication with scalar
  SystemVector& operator*=(SystemVector& x, double d);

  /// multiplication with a scalar
  inline SystemVector operator*(SystemVector x, double d) 
  { 
    return x *= d; 
  }

  /// multiplication with a scalar
  inline SystemVector operator*(double d, SystemVector x) 
  { 
    return x *= d; 
  }

  /// scalar product
  double operator*(SystemVector const& x, SystemVector const& y);

  /// addition of two system vectors
  SystemVector& operator+=(SystemVector& x, SystemVector const& y);

  /// subtraction of two system vectors.
  SystemVector& operator-=(SystemVector& x, SystemVector const& y);

  /// addition of two system vectors
  inline SystemVector operator+(SystemVector x, SystemVector const& y) 
  { 
    return x += y; 
  }

  /// addition of two system vectors
  inline SystemVector operator-(SystemVector x, SystemVector const& y) 
  { 
    return x -= y; 
  }


  /// Norm of system vector.
  double norm(SystemVector const* x);

  /// L2 norm of system vector.
  double L2Norm(SystemVector const* x);

  /// H1 norm of system vector.
  double H1Norm(SystemVector const* x);


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

  /// Returns SystemVector::getUsedSize().
  inline size_t size(SystemVector const& vec)
  {
    return static_cast<size_t>(vec.getUsedSize());
  }

  /// Returns SystemVector::getUsedSize().
  inline size_t size(SystemVector const* vec)
  {
    return size(*vec);
  }

} // end namespace AMDiS
