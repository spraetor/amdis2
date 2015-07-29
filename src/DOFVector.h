/** \file DOFVector.h */

#pragma once
 
#include <vector> 
#include <memory> 
#include <map> 

#include "AMDiS_fwd.h"
#include "FixVec.h"
#include "Global.h" 
#include "Flag.h" 
#include "RCNeighbourList.h" 
#include "DOFIterator.h"
#include "DOFIndexed.h"
#include "DOFContainer.h"
#include "Boundary.h"
#include "CreatorInterface.h"
#include "DOFMatrix.h" 
#include "BasisFunction.h"
#include "FiniteElemSpace.h"
#include "SurfaceQuadrature.h"
#include "Traits.h"
// #include "Expressions.h"

namespace AMDiS 
{
  template <class T> 
  class DOFVectorBase : public DOFIndexed<T>
  {
  public:

    DOFVectorBase() 
      : feSpace(NULL),
	elementVector(3),
        boundaryManager(NULL),
        nBasFcts(0)
    {}
    
    DOFVectorBase(const FiniteElemSpace *f, std::string n);

    virtual ~DOFVectorBase();

    /// For the given element, this function returns an array of all DOFs of 
    /// this DOFVector that are defined on this element.
    virtual void getLocalVector(const Element *el, 
				DenseVector<T>& localVec) const;

    /// Evaluates the DOF vector at a set of quadrature points defined on the 
    /// given element.
    void getVecAtQPs(const ElInfo *elInfo, 
		     const Quadrature *quad,
		     const FastQuadrature *quadFast,
		     DenseVector<T>& vecAtQPs) const;

    /// Evaluates the gradient of a DOF vector at a set of quadrature points
    /// defined on the given element.
    void getGrdAtQPs( const ElInfo *elInfo,
		      const Quadrature *quad,
		      const FastQuadrature *quadFast,
		      DenseVector<Gradient_t<T>> &grdAtQPs) const;

    /// Evaluates the comp'th component of the derivative of a DOF vector at a
    /// set of quadrature points defined on the given element.
    void getDerivativeAtQPs( const ElInfo *elInfo,
		      const Quadrature *quad,
		      const FastQuadrature *quadFast,
		      int comp,
		      DenseVector<T> &derivativeAtQPs) const;

    /// Evaluates the jacobian of a DOF vector at a set of quadrature points
    ///  defined on the given element.
    void getD2AtQPs(const ElInfo *elInfo,
		    const Quadrature *quad,
		    const FastQuadrature *quadFast,
		    DenseVector<typename D2Type<T>::type> &d2AtQPs) const;

    /// Returns the FE space the DOF vector is defined on.
    inline const FiniteElemSpace* getFeSpace() const 
    {
      return feSpace;
    }

    /// Assembles the element vector for the given ellement and adds the
    /// element matrix to the current DOF vector.
    void assemble(T factor, ElInfo *elInfo,			    
            		  const BoundaryType *bound, 
            		  Operator *op = NULL);

    void addElementVector(T sign,
                  			  const ElementVector& elVec, 
                  			  const BoundaryType *bound,
                  			  ElInfo *elInfo,
                  			  bool add = true); 

    ///
    void assembleOperator(Operator &op);
    
    /// That function must be called after the matrix assembling has been
    /// finished. This makes it possible to start some cleanup or matrix
    /// data compressing procedures.
    void finishAssembling();
 
    inline void addOperator(Operator* op, 
			    double *factor = NULL,
			    double *estFactor = NULL) 
    {
      operators.push_back(op);
      operatorFactor.push_back(factor);
      operatorEstFactor.push_back(estFactor);
    }

    inline std::vector<double*>::iterator getOperatorFactorBegin() 
    {
      return operatorFactor.begin();
    }

    inline std::vector<double*>::iterator getOperatorFactorEnd() 
    {
      return operatorFactor.end();
    }

    inline std::vector<double*>::iterator getOperatorEstFactorBegin() 
    {
      return operatorEstFactor.begin();
    }

    inline std::vector<double*>::iterator getOperatorEstFactorEnd() 
    {
      return operatorEstFactor.end();
    }


    inline std::vector<Operator*>::iterator getOperatorsBegin() 
    {
      return operators.begin();
    }

    inline std::vector<Operator*>::iterator getOperatorsEnd() 
    {
      return operators.end();
    }

    Flag getAssembleFlag();

    /// Evaluates \f[ u_h(x(\lambda)) = \sum_{i=0}^{m-1} vec[ind[i]] * 
    /// \varphi^i(\lambda) \f] where \f$ \varphi^i \f$ is the i-th basis
    /// function, \f$ x(\lambda) \f$ are the world coordinates of lambda
    /// and \f$ m \f$ is the number of basis functions
    T evalUh(const DimVec<double>& lambda, DegreeOfFreedom* ind);

    inline std::vector<Operator*>& getOperators() 
    { 
      return operators; 
    }

    inline std::vector<double*>& getOperatorFactor() 
    { 
      return operatorFactor; 
    }

    inline std::vector<double*>& getOperatorEstFactor() 
    { 
      return operatorEstFactor; 
    }

    /// Returns \ref name
    inline std::string getName() const 
    { 
      return name; 
    } 

    inline void setName(std::string n)
    {
      name = n;
    }

    inline BoundaryManager* getBoundaryManager() const 
    { 
      return boundaryManager; 
    }

    inline void setBoundaryManager(BoundaryManager *bm) 
    {
      boundaryManager = bm;
    }

    inline void setDirichletDofValue(DegreeOfFreedom dof,
				     T value)
    {
      dirichletDofValues[dof] = value;
    }

    std::map<DegreeOfFreedom, T>& getDirichletValues()
    {
      return dirichletDofValues;
    }

  protected:
    ///
    const FiniteElemSpace *feSpace;

    ///
    std::string name;

    ///
    ElementVector elementVector;

    ///
    std::vector<Operator*> operators;

    ///
    std::vector<double*> operatorFactor;

    ///
    std::vector<double*> operatorEstFactor;

    ///
    BoundaryManager *boundaryManager;

    /// Number of basis functions of the used finite element space.
    int nBasFcts;

    /// Dimension of the mesh this DOFVectorBase belongs to
    int dim;

    std::map<DegreeOfFreedom, T> dirichletDofValues;
  };



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
    typedef typename DOFVectorBase<T>::value_type	value_type;
    typedef typename DOFVectorBase<T>::size_type	size_type;
    typedef typename DOFVectorBase<T>::reference	reference;
    typedef typename DOFVectorBase<T>::const_reference	const_reference;
  
  public:
    /** \ingroup DOFAdministration
     * \brief
     * Enables the access of DOFVector<T>::Iterator. Alias for DOFIterator<T>
     */
    class Iterator : public DOFIterator<T> {
    public:
      Iterator(DOFIndexed<T> *c, DOFIteratorType type)
	: DOFIterator<T>(c, type)
      {}

      Iterator(DOFAdmin *admin, DOFIndexed<T> *c, DOFIteratorType type)
	: DOFIterator<T>(admin, c, type)
      {}
    };

    class Creator : public CreatorInterface<DOFVector<T> > {
    public:
      Creator(FiniteElemSpace *feSpace_) 
        : feSpace(feSpace_) 
      {}

      DOFVector<T> *create() 
      { 
	return new DOFVector<T>(feSpace, ""); 
      }

      void free(DOFVector<T> *vec) 
      { 
	delete vec; 
      }

    private:
      FiniteElemSpace *feSpace;
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
    virtual ~DOFVector();

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
				    std::vector<DegreeOfFreedom> &newDof);
    
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
      return vec.size();
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
    double IntOnBoundary(BoundaryType boundary, Quadrature* q = NULL) const
    {
      FUNCNAME("DOFVector::IntOnBoundary())");
      TEST_EXIT(false)("Please implement your integration\n");
      return 0.0;
    }

    /// Calculates Integral of this DOFVector times normal vector over parts 
    /// of the domain boundary, indicated by boundaryType. Implemented for 
    /// DOFVector<WorldVector<double> >
    double IntOnBoundaryNormal(BoundaryType boundary, Quadrature* q = NULL) const
    {
      FUNCNAME("DOFVector::IntOnBoundaryNormal())");
      TEST_EXIT(false)("Please implement your integration\n");
      return 0.0;
    }

    /// Calculates L1 norm of this DOFVector
    double L1Norm(Quadrature* q = NULL) const;
 
    /// Calculates L2 norm of this DOFVector
    double L2Norm(Quadrature* q = NULL) const 
    {
      return std::sqrt(L2NormSquare());
    }

    /// Calculates square of L2 norm of this DOFVector
    double L2NormSquare(Quadrature* q = NULL) const;

    /// Calculates H1 norm of this DOFVector
    double H1Norm(Quadrature* q = NULL) const 
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
    DOFVector<T>& operator=(T d) 
    {
      set(d); 
      return *this;
    } 

    /// Assignment operator between two vectors
    DOFVector<T>& operator=(const DOFVector<T>&); 

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
    template <class Expr>
    inline void interpol(Expr expr);
    
    // implementation in Expressions.h
    inline void interpol(std::function<T(WorldVector<double>)> f);

    // implementation in Expressions.h
    inline void interpol(DOFVector<T> *v, double factor = 1.0);


    /// Eval DOFVector at given point p. If oldElInfo != NULL the search for 
    /// the element, where p is inside, starts from oldElInfo. implemented for:
    /// double, WorldVector< double >
    T evalAtPoint(WorldVector<double> &p, 
		  ElInfo *oldElInfo = NULL) const 
    {
      FUNCNAME("DOFVector::evalAtPoint())");
      TEST_EXIT(false)("Please implement your evaluation\n");
    }

    /// Determine the DegreeOfFreedom that has coords with minimal euclidean 
    /// distance to WorldVector p. return true if DOF is found, and false 
    /// otherwise.
    bool getDofIdxAtPoint(WorldVector<double> &p, 
				DegreeOfFreedom &idx, 
				ElInfo *oldElInfo = NULL, 
				bool useOldElInfo = false) const;


    DOFVector<Gradient_t<T>>* getGradient(DOFVector<Gradient_t<T>> *grad) const;

    WorldVector<DOFVector<T>*> *getGradient(WorldVector<DOFVector<T>*> *grad) const;

    DOFVector<Gradient_t<T>>* getRecoveryGradient(DOFVector<Gradient_t<T>> *grad) const;

  protected: 

    /// Data container
    std::vector<T> vec; 
  }; 


  template<>
  double DOFVector<double>::IntOnBoundary(
    BoundaryType boundaryType, Quadrature* q) const;

  template<>
  double DOFVector<WorldVector<double> >::IntOnBoundaryNormal(
    BoundaryType boundaryType, Quadrature* q) const;

  template<>
  double DOFVector<double>::evalAtPoint(WorldVector<double> &p, 
					 ElInfo *oldElInfo) const;

  template<>
  WorldVector<double> DOFVector<WorldVector<double> >::evalAtPoint(WorldVector<double> &p, 
								   ElInfo *oldElInfo) const;

  template<>
  void DOFVector<double>::refineInterpol(RCNeighbourList&, int);

  template<>
  void DOFVector<double>::coarseRestrict(RCNeighbourList&, int); 

  inline double min(const DOFVector<double>& v) 
  {
    return v.min();
  } 

  inline double max(const DOFVector<double>& v) 
  {
    return v.max();
  }


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

  template<typename T>
  double norm(DOFVector<T> *vec) 
  {
    return vec->nrm2();
  }

  template<typename T>
  double L2Norm(DOFVector<T> *vec) 
  {
    return vec->L2Norm();
  }

  template<typename T>
  double H1Norm(DOFVector<T> *vec) 
  {
    return vec->H1Norm();
  }

  template<typename T>
  void print(DOFVector<T> *vec) 
  {
    vec->print();
  }
  
  
  // point wise multiplication
  template<typename T>
  const DOFVector<T>& operator*=(DOFVector<T>& x, const DOFVector<T>& y);

  // multiplication with scalar
  template<typename T>
  const DOFVector<T>& operator*=(DOFVector<T>& x, T scal);

  // scalar product
  template<typename T>
  T operator*(DOFVector<T>& x, DOFVector<T>& y);

  // addition
  template<typename T>
  const DOFVector<T>& operator+=(DOFVector<T>& x, const DOFVector<T>& y);

  // subtraction
  template<typename T>
  const DOFVector<T>& operator-=(DOFVector<T>& x, const DOFVector<T>& y);

  template<typename T>
  const DOFVector<T>& operator*(const DOFVector<T>& v, double d);

  template<typename T>
  const DOFVector<T>& operator*(double d, const DOFVector<T>& v);

  template<typename T>
  const DOFVector<T>& operator+(const DOFVector<T>&v1 , const DOFVector<T>& v2);


  template<typename T>
  inline void set(DOFVector<T>& vec, T d) 
  {
    vec.set(d);
  }

  template<typename T>
  inline void setValue(DOFVector<T>& vec, T d) 
  {
    vec.set(d);
  }

  template<typename T>
  inline void checkFeSpace(const FiniteElemSpace* feSpace, const std::vector<T>& vec)
  {
    FUNCNAME_DBG("checkFeSpace()");
    TEST_EXIT_DBG(feSpace)("feSpace is NULL\n");
    TEST_EXIT_DBG(feSpace->getAdmin())("admin is NULL\n");
    TEST_EXIT_DBG(static_cast<int>(vec.size()) >= feSpace->getAdmin()->getUsedSize())
      ("size = %d too small: admin->sizeUsed = %d\n", vec.size(),
       feSpace->getAdmin()->getUsedSize());
  }

  WorldVector<DOFVector<double>*> *transform(DOFVector<WorldVector<double> > *vec,
					     WorldVector<DOFVector<double>*> *result);
  
  template<typename T>
  std::vector<DOFVector<double>*> *transform(DOFVector<Gradient_t<T>> *vec,
					     std::vector<DOFVector<double>*> *res);
  
} // end namespace AMDiS

#include "DOFVector.hh"
