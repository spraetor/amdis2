/** \file DOFVector.h */

#pragma once

#include <vector>
#include <memory>
#include <map>

#include "AMDiS_fwd.hpp"
#include "BasisFunction.hpp"
#include "Boundary.hpp"
#include "CreatorInterface.hpp"
#include "DOFContainer.hpp"
#include "DOFIndexed.hpp"
#include "DOFIterator.hpp"
#include "DOFMatrix.hpp"
#include "DOFVectorBase.hpp"
#include "FiniteElemSpace.hpp"
#include "FixVec.hpp"
#include "Flag.hpp"
#include "Global.hpp"
#include "RCNeighbourList.hpp"
#include "SurfaceQuadrature.hpp"
#include "Traits.hpp"

namespace AMDiS
{
  /** \ingroup DOFAdministration
   * \brief
   * The DOFs described above are just integers that can be used as indices into
   * vectors and matrices. During refinement and coarsening of the mesh, the
   * number of used DOFs, the meaning of one integer index, and even the total
   * range of DOFs change. To be able to handle these changes automatically for
   * all vectors, which are indexed by the DOFs, special data structures are
   * used which contain such vector data. Lists of these structures are kept in
   * DOFAdmin, so that all vectors in the lists can be resized together with the
   * range of DOFs. During refinement and coarsening of elements, values can be
   * interpolated automatically to new DOFs, and restricted from old DOFs.
   */
  template <class T>
  class DOFVector : public DOFVectorBase<T>
  {
  public:
    using Super = DOFVectorBase<T>;
    using value_type = Value_t<Super>;
    using size_type  = Size_t<Super>;
    using reference  = typename Super::reference;
    using const_reference = typename Super::const_reference;

  public:
    /** \ingroup DOFAdministration
     * \brief
     * Enables the access of DOFVector<T>::Iterator. Alias for DOFIterator<T>
     */
    class Iterator : public DOFIterator<T>
    {
    public:
      Iterator(DOFIndexed<T>* c, DOFIteratorType type)
        : DOFIterator<T>(c, type)
      {}

      Iterator(DOFAdmin* admin, DOFIndexed<T>* c, DOFIteratorType type)
        : DOFIterator<T>(admin, c, type)
      {}
    };

    class Creator : public CreatorInterface<DOFVector<T>>
    {
    public:
      Creator(FiniteElemSpace* feSpace_)
        : feSpace(feSpace_)
      {}

      DOFVector<T>* create()
      {
        return new DOFVector<T>(feSpace, "");
      }

      void free(DOFVector<T>* vec)
      {
        delete vec;
      }

    private:
      FiniteElemSpace* feSpace;
    };

  public:
    /// Empty constructor. No initialization!
    DOFVector()
      : DOFVectorBase<T>()
    {}

    /// Constructs a DOFVector with name n belonging to FiniteElemSpace f
    DOFVector(const FiniteElemSpace* f, std::string n, bool addToSynch = false);

    /// Initialization.
    void init(const FiniteElemSpace* f, std::string n, bool addToSynch = false);

    /// Copy Constructor
    DOFVector(const DOFVector& rhs) : DOFVectorBase<T>()
    {
      *this = rhs;
      this->name = rhs.name + "copy";
      if (this->feSpace && this->feSpace->getAdmin())
        (dynamic_cast<DOFAdmin*>(this->feSpace->getAdmin()))->addDOFIndexed(this);
    }

    /// Destructor
    ~DOFVector();

    /// Returns iterator to the begin of \ref vec
    typename std::vector<T>::iterator begin()
    {
      return vec.begin();
    }

    /// Returns iterator to the end of \ref vec
    typename std::vector<T>::iterator end()
    {
      return vec.end();
    }

    /// Returns const_iterator to the begin of \ref vec
    typename std::vector<T>::const_iterator begin() const
    {
      return vec.begin();
    }

    /// Returns const_iterator to the end of \ref vec
    typename std::vector<T>::const_iterator end() const
    {
      return vec.end();
    }

    /// Used by DOFAdmin to compress this DOFVector. Implementation of
    /// \ref DOFIndexedBase::compress()
    virtual void compressDOFIndexed(int first, int last,
                                    std::vector<DegreeOfFreedom>& newDof);

    /// Restriction after coarsening. Implemented for DOFVector<double>
    void coarseRestrict(RCNeighbourList&, int) {}

    /// Interpolation after refinement. Implemented for DOFVector<double>
    void refineInterpol(RCNeighbourList&, int) {}

    /// Returns \ref vec
    std::vector<T>& getVector()
    {
      return vec;
    }

    /// Returns size of \ref vec
    int getSize() const
    {
      return int(vec.size());
    }

    /// Returns used size of the vector.
    int getUsedSize() const
    {
      return this->feSpace->getAdmin()->getUsedSize();
    }

    /// Resizes \ref vec to n
    void resize(int n)
    {
      FUNCNAME_DBG("DOFVector<T>::resize()");
      TEST_EXIT_DBG(n >= 0)("Can't resize DOFVector to negative size\n");
      vec.resize(n);
    }

    /// Resizes \ref vec to n and inits new values with init
    void resize(int n, T init)
    {
      FUNCNAME_DBG("DOFVector<T>::resize()");
      TEST_EXIT_DBG(n >= 0)("Can't resize DOFVector to negative size\n");
      vec.resize(n, init);
    }

    /// Returns \ref vec[i]
    const_reference operator[](DegreeOfFreedom i) const
    {
      FUNCNAME_DBG("DOFVector<T>::operator[]");
      TEST_EXIT_DBG(i >= 0 && i < static_cast<int>(vec.size()))
      ("Illegal vector index %d.\n", i);
      return vec[i];
    }

    /// Returns \ref vec[i]
    reference operator[](DegreeOfFreedom i)
    {
      FUNCNAME_DBG("DOFVector<T>::operator[]");

      TEST_EXIT_DBG(i >= 0 && i < static_cast<int>(vec.size()))
      ("Illegal vector index %d.\n", i);

      return vec[i];
    }

    /// Calculates Integral of this DOFVector
    T Int(Quadrature* q = NULL) const
    {
      return Int(-1, q);
    }

    /** \brief
     * Calculates Integral of this DOFVector.
     *
     * \param[in]  meshlevel  Then mesh level on which the integral should be
     *                        calculated. Usually, only -1 for the leaf level
     *                        makes sence, but other values can be used, e.g.,
     *                        for debugging.
     * \param[in]  q          Quadrature object. If not specified, the function
     *                        creates a new quadrature object.
     */
    T Int(int meshLevel, Quadrature* q = NULL) const;


    /// Calculates Integral of this DOFVector over parts of the domain
    /// boundary, indicated by boundaryType. Implemented for DOFVector<double>
    double IntOnBoundary(BoundaryType /*boundary*/, Quadrature* /*q*/ = NULL) const
    {
      FUNCNAME("DOFVector::IntOnBoundary())");
      TEST_EXIT(false)("Please implement your integration\n");
      return 0.0;
    }

    /// Calculates Integral of this DOFVector times normal vector over parts
    /// of the domain boundary, indicated by boundaryType. Implemented for
    /// DOFVector<WorldVector<double> >
    double IntOnBoundaryNormal(BoundaryType /*boundary*/, Quadrature* /*q*/ = NULL) const
    {
      FUNCNAME("DOFVector::IntOnBoundaryNormal())");
      TEST_EXIT(false)("Please implement your integration\n");
      return 0.0;
    }

    /// Calculates L1 norm of this DOFVector
    double L1Norm(Quadrature* q = NULL) const;

    /// Calculates L2 norm of this DOFVector
    double L2Norm(Quadrature* /*q*/ = NULL) const
    {
      return std::sqrt(L2NormSquare());
    }

    /// Calculates square of L2 norm of this DOFVector
    double L2NormSquare(Quadrature* q = NULL) const;

    /// Calculates H1 norm of this DOFVector
    double H1Norm(Quadrature* /*q*/ = NULL) const
    {
      return std::sqrt(H1NormSquare());
    }

    /// Calculates square of H1 norm of this DOFVector
    double H1NormSquare(Quadrature* q = NULL) const;

    /// Calculates euclidian norm of this DOFVector
    double nrm2() const;

    /// Returns square of the euclidian norm.
    double squareNrm2() const;

    /// Calculates l2 norm of this DOFVector
    double l2norm() const
    {
      return nrm2();
    }

    /// Calculates the absolute sum of this DOFVector
    T asum() const;

    /// Calculates the l1 norm of this DOFVector
    double l1norm() const
    {
      return asum();
    }

    /// Calculates doublewell of this DOFVector
    double DoubleWell(Quadrature* q = NULL) const;

    /// Calculates the sum of this DOFVector
    T sum() const;

    /// Sets \ref vec[i] = val, i=0 , ... , size
    void set(T val);

    /// Assignment operator for setting complete vector to a certain value d
    DOFVector<T>& operator=(T scal)
    {
      set(scal);
      return *this;
    }

    /// Assignment operator between two vectors
    DOFVector<T>& operator=(DOFVector<T> const& y);

    // point wise multiplication
    DOFVector<T>& operator*=(DOFVector<T> const& y);

    // multiplication with scalar
    DOFVector<T>& operator*=(T scal);

    // addition
    DOFVector<T>& operator+=(DOFVector<T> const& y);

    // subtraction
    DOFVector<T>& operator-=(DOFVector<T> const& y);


    /// vec[i] = v.vec[i]
    void copy(const DOFVector<T>& v);

    /// Returns minimum of DOFVector.
    T min() const;

    /// Returns maximum of DOFVector.
    T max() const;

    /// Returns absolute maximum of DOFVector.
    T absMax() const;

    /// Returns the average value of the DOFVector.
    T average() const;

    ///
    size_t calcMemoryUsage() const;

    /// Computes the coefficients of the interpolant of the function fct and
    /// stores these in the DOFVector
    // implementation in Expressions.h
    template <class Term>
    inline void interpol(Term const& term);

    // implementation in Expressions.h
    inline void interpol(DOFVector<T>* v, double factor = 1.0);


    /// Eval DOFVector at given point p. If oldElInfo != NULL the search for
    /// the element, where p is inside, starts from oldElInfo. implemented for:
    /// double, WorldVector< double >
    T evalAtPoint(WorldVector<double> const& /*p*/,
                  ElInfo* /*oldElInfo*/ = NULL) const
    {
      FUNCNAME("DOFVector::evalAtPoint())");
      TEST_EXIT(false)("Please implement your evaluation\n");
    }

    T operator()(WorldVector<double> const& p) const
    {
      return evalAtPoint(p, NULL);
    }

    /// Determine the DegreeOfFreedom that has coords with minimal euclidean
    /// distance to WorldVector p. return true if DOF is found, and false
    /// otherwise.
    bool getDofIdxAtPoint(WorldVector<double>& p,
                          DegreeOfFreedom& idx,
                          ElInfo* oldElInfo = NULL,
                          bool useOldElInfo = false) const;


    DOFVector<Gradient_t<T>>* getGradient(DOFVector<Gradient_t<T>>* grad) const;

    WorldVector<DOFVector<T>*>* getGradient(WorldVector<DOFVector<T>*>* grad) const;

    DOFVector<Gradient_t<T>>* getRecoveryGradient(DOFVector<Gradient_t<T>>* grad) const;

  protected:
    /// Data container
    std::vector<T> vec;
  };

  template<>
  double DOFVector<double>::IntOnBoundary(
    BoundaryType boundaryType, Quadrature* q) const;

  template<>
  double DOFVector<WorldVector<double>>::IntOnBoundaryNormal(
                                       BoundaryType boundaryType, Quadrature* q) const;

  template<>
  double DOFVector<double>::evalAtPoint(WorldVector<double> const& p,
                                        ElInfo* oldElInfo) const;

  template<>
  WorldVector<double> DOFVector<WorldVector<double>>::evalAtPoint(WorldVector<double> const& p,
      ElInfo* oldElInfo) const;

  template<>
  void DOFVector<double>::refineInterpol(RCNeighbourList&, int);

  template<>
  void DOFVector<double>::coarseRestrict(RCNeighbourList&, int);



  /** \ingroup DOFAdministration
   * \brief
   * A DOFVector that stores DOF indices.
   */
  class DOFVectorDOF : public DOFVector<DegreeOfFreedom>,
    public DOFContainer
  {
  public:
    /// Calls constructor of DOFVector<DegreeOfFreedom> and registers itself
    /// as DOFContainer at DOFAdmin
    DOFVectorDOF(const FiniteElemSpace* feSpace_, std::string name_)
      : DOFVector<DegreeOfFreedom>(feSpace_, name_)
    {
      feSpace->getAdmin()->addDOFContainer(this);
    }

    /// Deregisters itself at DOFAdmin.
    ~DOFVectorDOF()
    {
      feSpace->getAdmin()->removeDOFContainer(this);
    }

    /// Implements DOFContainer::operator[]() by calling
    /// DOFVector<DegreeOfFreedom>::operator[]()
    DegreeOfFreedom& operator[](DegreeOfFreedom i)
    {
      return DOFVector<DegreeOfFreedom>::operator[](i);
    }

    const DegreeOfFreedom& operator[](DegreeOfFreedom i) const
    {
      return DOFVector<DegreeOfFreedom>::operator[](i);
    }

    /// Implements DOFIndexedBase::getSize()
    int getSize() const
    {
      return DOFVector<DegreeOfFreedom>::getSize();
    }

    /// Implements DOFIndexedBase::resize()
    void resize(int size)
    {
      DOFVector<DegreeOfFreedom>::resize(size);
    }

    void freeDOFContent(DegreeOfFreedom dof);

  protected:
    DOFVectorDOF();
  };

  // ===========================================================================

  template <class T>
  inline void checkFeSpace(const FiniteElemSpace* feSpace, const std::vector<T>& vec)
  {
    FUNCNAME_DBG("checkFeSpace()");
    TEST_EXIT_DBG(feSpace)("feSpace is NULL\n");
    TEST_EXIT_DBG(feSpace->getAdmin())("admin is NULL\n");
    TEST_EXIT_DBG(static_cast<int>(vec.size()) >= feSpace->getAdmin()->getUsedSize())
    ("size = %d too small: admin->sizeUsed = %d\n", vec.size(),
     feSpace->getAdmin()->getUsedSize());
  }

  WorldVector<DOFVector<double>*>* transform(DOFVector<WorldVector<double>>* vec,
      WorldVector<DOFVector<double>*>* result);

  template <class T>
  std::vector<DOFVector<double>*>* transform(DOFVector<Gradient_t<T>>* vec,
      std::vector<DOFVector<double>*>* res);



  template <class Value>
  inline size_t size(DOFVector<Value> const& v)
  {
    return v.getUsedSize();
  }

  template <class Value>
  inline size_t num_rows(DOFVector<Value> const& v)
  {
    return v.getUsedSize();
  }

  template <class Value>
  inline size_t num_cols(DOFVector<Value> const& /*v*/)
  {
    return 1;
  }


} // end namespace AMDiS


#include "DOFVectorBase.hh"
#include "DOFVector.hh"
