/** \file BasisFunction.h */

#pragma once

#include <string>

#include "AMDiS_fwd.h"
#include "CreatorInterface.h"
#include "Boundary.h"
#include "FixVec.h"

namespace AMDiS 
{
  /// Function interface for evaluating basis functions.
  struct BasFctType
  {
    virtual ~BasFctType() {}
    virtual double operator()(const DimVec<double>&) const = 0;
  };


  /// Function interface for evaluating gradients of basis functions. 
  struct GrdBasFctType
  {
    virtual ~GrdBasFctType() {}
    virtual void operator()(const DimVec<double>&, 
                  			    DenseVector<double>&) const = 0;
  };
  

  /// Function interface for evaluating second derivative of basis functions.
  struct D2BasFctType
  {
    virtual ~D2BasFctType() {}
    virtual void operator()(const DimVec<double>&, DimMat<double>&) const = 0;
  };
  
  
  typedef BasFctType *BFptr;
  typedef GrdBasFctType *GBFptr;
  typedef D2BasFctType *DBFptr;

  /** \ingroup FEMSpace
   * \brief
   * Base class for finite element basis functions. In order to build up a
   * finite element space, we have to specify a set of local basis functions.
   * Together with the correspondig DOF administration and the underlying mesh,
   * the finite element space is given. 
   * This class holds the local basis functions and their derivatives of the
   * reference element. They are evaluated at barycentric coordinates, so they
   * can be used on every element of the mesh.  
   */
  class BasisFunction
  {  
  protected:
    /// Creates a BasisFunction object of given dim and degree 
    BasisFunction(std::string name, int dim, int degree);

    /// destructor
    virtual ~BasisFunction();

  public:
    /// compares two BasisFunction objects.
    virtual bool operator==(const BasisFunction& a) const 
    {
      return a.getName() == name;
    }

    /// Returns !(*this == b)
    bool operator!=(const BasisFunction& b) const 
    {
      return !(operator == (b));
    }

    /// Used by \ref getDOFIndices and \ref getVec
    virtual int* orderOfPositionIndices(const Element* el, GeoIndex position, 
					                              int positionIndex) const = 0;

    /** \brief
     * The second argument 'bound' has to be a pointer to a vector which has 
     * to be filled. Its length is \ref nBasFcts (the number of basis functions
     * in the used finite element space). After calling this function, the i-th 
     * entry of the array is the boundary type of the i-th basis function of this
     * element.
     * 
     * This function needs boundary information within the ElInfo object; thus, 
     * all routines using this function on the elements need the FILL_BOUND 
     * flag during mesh traversal;
     */
    virtual void getBound(const ElInfo*, BoundaryType*) const {}

    /// Returns \ref degree of BasisFunction
    int getDegree() const 
    { 
      return degree; 
    }

    /// Returns \ref dim of BasisFunction
    int getDim() const 
    { 
      return dim; 
    }

    /// Returns \ref nBasFcts which is the number of local basis functions
    int getNumber() const 
    { 
      return nBasFcts; 
    }

    /// Returns \ref name of BasisFunction
    std::string getName() const 
    { 
      return name; 
    }

    /// Returns \ref nDOF[i]
    int getNumberOfDofs(int i) const;

    /// Returns \ref nDOF
    DimVec<int>* getNumberOfDofs() const 
    { 
      return nDOF; 
    }

    /// Initialisation of the \ref nDOF vector. Must be implemented by sub classes
    virtual void setNDOF() = 0;

    /// Returns the barycentric coordinates of the i-th basis function.
    virtual DimVec<double> *getCoords(int i) const = 0;

    /** \brief
     * Fills a vector with interpolation coefficients of the
     * function f; if indices is a pointer to NULL, the coefficient for all 
     * basis functions are calculated and the i-th entry in the vector is the 
     * coefficient of the i-th basis function; if indices is non NULL, only the 
     * coefficients for a subset of the local basis functions has to be 
     * calculated; n is the number of those basis functions, indices[0], . . . 
     * , indices[n-1] are the local indices of the basis functions where the
     * coefficients have to be calculated, and the i-th entry in the return 
     * vector is then the coefficient of the indices[i]-th basis function; coeff 
     * may be a pointer to a vector which has to be filled.
     * such a function usually needs vertex coordinate information; thus, all 
     * routines using this function on the elements need the FILL COORDS flag 
     * during mesh traversal.
     * Must be implemented by sub classes.
     */
    virtual void interpol(const ElInfo *el_info, 
                  			  int n, 
                  			  const int *indices, 
                  			  std::function<double(WorldVector<double>)> f,
                  			  DenseVector<double> &coeff) const = 0;

    /// WorldVector<double> valued interpol function.
    virtual void interpol(const ElInfo *el_info, 
                  			  int no, 
                  			  const int *b_no,
                  			  std::function<WorldVector<double>(WorldVector<double>)> f, 
                  			  DenseVector<WorldVector<double> >& coeff) const = 0;

    /// Returns the i-th local basis function
    inline BasFctType *getPhi(int i) const 
    { 
      return (*phi)[i]; 
    }

    /// Returns the gradient of the i-th local basis function
    inline GrdBasFctType *getGrdPhi(int i) const 
    { 
      return (*grdPhi)[i]; 
    }

    /// Returns the second derivative of the i-th local basis function
    inline D2BasFctType *getD2Phi(int i) const 
    { 
      return (*d2Phi)[i]; 
    }

    /** \brief
     * Approximates the L2 scalar products of a given function with all basis 
     * functions by numerical quadrature and adds the corresponding values to a 
     * DOF vector;
     * f is a pointer for the evaluation of the given function in world 
     * coordinates x and returns the value of that function at x; if f is a NULL
     *  pointer, nothing is done;
     * fh is the DOF vector where at the i-th entry the approximation of the L2 
     * scalar product of the given function with the i-th global basis function 
     * of fh->feSpace is stored;
     * quad is the quadrature for the approximation of the integral on each leaf 
     * element of fh->feSpace->mesh; if quad is a NULL pointer, a default 
     * quadrature which is exact of degree 2*fh->feSpace->basFcts->degree-2 is 
     * used.
     * The integrals are approximated by looping over all leaf elements, 
     * computing the approximations to the element contributions and adding 
     * these values to the vector fh by add element vec().
     * The vector fh is not initialized with 0.0; only the new contributions are 
     * added
     */
    // TODO: add template function and move virtual function to private section
    virtual void l2ScpFctBas(Quadrature*,
                  			     std::function<double(WorldVector<double>)> /*f*/,
                  			     DOFVector<double>* /*fh*/)
    {}

    /// WorldVector<double> valued l2ScpFctBas function
    virtual void l2ScpFctBas(Quadrature* ,
                  			     std::function<WorldVector<double>(WorldVector<double>)> /*f*/,
                  			     DOFVector<WorldVector<double> >* /*fh*/) 
    {}


    /// Interpolates a DOFIndexed<double> after refinement
    virtual void refineInter(DOFIndexed<double> *, RCNeighbourList*, int)
    {}

    /// Interpolates a DOFIndexed<double> after coarsening
    virtual void coarseInter(DOFIndexed<double> *, RCNeighbourList*, int)
    {}

    /// Restricts a DOFIndexed<double> after coarsening
    virtual void coarseRestr(DOFIndexed<double> *, RCNeighbourList*, int)
    {}

    /// Interpolates a DOFVector<WorldVector<double> > after refinement
    virtual void refineInter(DOFVector<WorldVector<double> >*, RCNeighbourList*, int)
    {}

    /// Interpolates a DOFVector<WorldVector<double> > after coarsening
    virtual void coarseInter(DOFVector<WorldVector<double> >*, RCNeighbourList*, int)
    {}

    /// Restricts a DOFVector<WorldVector<double> > after coarsening
    virtual void coarseRestr(DOFVector<WorldVector<double> >*, RCNeighbourList*, int)
    {}

    /// Returns local dof indices of the element for the given fe space.
    virtual void getLocalIndices(const Element *el,
                        				 const DOFAdmin *admin,
                        				 std::vector<DegreeOfFreedom> &indices) const
    {}

    ///
    virtual void getLocalDofPtrVec(const Element *el, 
                        				   const DOFAdmin *admin,
                        				   std::vector<const DegreeOfFreedom*>& vec) const
    {}


    /// Evaluates elements value at barycentric coordinates lambda with local 
    /// coefficient vector uh.
    template <class T>
    T evalUh(const DimVec<double>& lambda, const DenseVector<T>& uh) const;


    /** \brief
     * Evaluates the gradient at barycentric coordinates lambda. Lambda is the
     * Jacobian of the barycentric coordinates. uh is the local coefficient
     * vector. If val is not NULL the result will be stored 
     * there, otherwise a pointer to a static local variable is returned which 
     * will be overwritten after the next call.
     */
    template <class T> Gradient_t<T>& 
    evalGrdUh(const DimVec<double>& lambda,
      	      const DimVec<WorldVector<double> >& Lambda,
      	      const DenseVector<T>& uh,
      	      Gradient_t<T>& val) const;


    /** \brief
     * Evaluates the second derivative at barycentric coordinates lambda. 
     * Lambda is the Jacobian of the barycentric coordinates. uh is the local 
     * coefficient vector. If val is not NULL the result will be stored 
     * there, otherwise a pointer to a static local variable is returned which 
     * will be overwritten after the next call.
     */
    const WorldMatrix<double>& 
    evalD2Uh(const DimVec<double>& lambda,
      	     const DimVec<WorldVector<double> >& Lambda,
      	     const ElementVector& uh,
      	     WorldMatrix<double>* val) const;

    /**
    * override this method, if the base of your finite element space is not
    * nodal
    */
    virtual bool isNodal() const = 0;

  protected:
    /// Textual description
    std::string name;     

    /// Number of basisfunctions on one Element                 
    int nBasFcts;

    /// Maximal degree of the basis functions                 
    int degree;

    /// Dimension of the basis functions                  
    int dim;

    /// Dimension of the world.
    int dow;

    /// Number of DOFs at the different positions                  
    DimVec<int> *nDOF;

    /// Vector of the local functions
    std::vector<BasFctType*> *phi;

    /// Vector of gradients
    std::vector<GrdBasFctType*> *grdPhi;

    /// Vector of second derivatives
    std::vector<D2BasFctType*> *d2Phi;
  };
  
  /** 
   * \brief
   * Interface for creators of concrete BasisFunctions. 
   */
  class BasisFunctionCreator : public CreatorInterface<BasisFunction>
  { 
  public:
    virtual ~BasisFunctionCreator() {}

    /// Sets \ref problem
    void setDim(int dim_) 
    { 
      dim = dim_; 
    }

  protected:
    /// dimension of the mesh
    int dim;
  };

} // end namespace AMDiS

#include "BasisFunction.hh"
