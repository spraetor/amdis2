/** \file OperatorTerm.h */

#pragma once

#include <set>

#include "AMDiS_fwd.h"
#include "SubAssembler.h"

namespace AMDiS 
{

  /** 
   * \ingroup Assembler
   * 
   * \brief
   * Base class for ZeroOrderTerm, FirstOrderTerm and SecondOrderTerm. 
   * OperatorTerms are the building blocks of an Operator. Each OperatorTerm
   * has its properties which are regarded, when constructing 
   * an Assembler for the corresponding Operator.
   */
  class OperatorTerm
  {
  public:
    /// Constructs an OperatorTerm with initially no properties.
    /// degree_ is used to determine the degree of the needed quadrature
    /// for the assemblage.  
    OperatorTerm(int deg) 
      : properties(0), 
      	degree(deg),
      	dimOfWorld(Global::getGeo(WORLD)),
      	bOne(-1)
    {}

    /// Destructor.
    virtual ~OperatorTerm() {}

    /// Virtual method. It's called by SubAssembler::initElement() for
    /// each OperatorTerm belonging to this SubAssembler. E.g., vectors
    /// and coordinates at quadrature points can be calculated here.
    void initElement(const ElInfo* elInfo, SubAssembler* subAssembler, 
		                 Quadrature *quad = NULL) 
    {
      initImpl(elInfo, subAssembler, quad);
    }

    /// Returs \auxFeSpaces, the list of all aux fe spaces the operator makes 
    /// use off.
    std::set<const FiniteElemSpace*>& getAuxFeSpaces()
    {
      return auxFeSpaces;
    }

    /// Specifies whether the matrix of the term is symmetric
    void setSymmetric(bool symm);

    /// Returns true, if the term is piecewise constant, returns false otherwise.
    bool isPWConst() const
    { 
      return (degree == 0);
    }

    /// Returns true, if the term has a symmetric matrix, returns false otherwise.
    bool isSymmetric();

    /// Returns \ref degree.
    int getDegree() const
    { 
      return degree; 
    }

    /// Sets one component of the b vector to be one. See \ref bOne.
    void setB(int b)
    {
      bOne = b;
    }

    /// Evaluation of the OperatorTerm at all quadrature points.
    void eval(int nPoints,
      	      const mtl::dense_vector<double>& uhAtQP,
      	      const mtl::dense_vector<WorldVector<double> >& grdUhAtQP,
      	      const mtl::dense_vector<WorldMatrix<double> >& D2UhAtQP,
      	      mtl::dense_vector<double>& result,
      	      double factor) const
    {
      evalImpl(nPoints, uhAtQP, grdUhAtQP, D2UhAtQP, result, factor);
    }
    
  private:
    // default behavior: init nothing
    virtual void initImpl(const ElInfo*, SubAssembler*, Quadrature*) { }
    
    // must be implemented by derived class
    virtual void evalImpl(int nPoints,
                  			  const mtl::dense_vector<double>& uhAtQP,
                  			  const mtl::dense_vector<WorldVector<double> >& grdUhAtQP,
                  			  const mtl::dense_vector<WorldMatrix<double> >& D2UhAtQP,
                  			  mtl::dense_vector<double>& result,
                  			  double factor) const = 0;
    
  protected:
    /// Stores the properties of this OperatorTerm
    Flag properties;

    /// Polynomial degree of the term. Used to detemine the degree of the quadrature.
    int degree;

    /// Stores the dimension of the world.
    int dimOfWorld;

    /// List off all fe spaces, the operator term makes use off.
    std::set<const FiniteElemSpace*> auxFeSpaces;

    /// Pointer to the Operator this OperatorTerm belongs to.
    Operator* operat;

    /// In many cases, the vector b in the evaluation \f$ \Lambda \cdot b\f$ has
    /// zeros in all components expect one that is set to one. Using the function
    /// \ref lb is then unnecessary time consuming. Instead, this variable
    /// defines the component of the vector b to be one. The function \ref lb_one
    /// is used if this variable is not -1.
    int bOne;

    /// Flag for piecewise constant terms
    static const Flag PW_CONST;

    /// Flag for symmetric terms
    static const Flag SYMMETRIC;

    friend class SubAssembler;
    friend class ZeroOrderAssembler;
    friend class FirstOrderAssembler;
    friend class SecondOrderAssembler;
    friend class Operator;
  };
  
  
  // forward declarations
  class ZeroOrderTerm;
  class FirstOrderTerm;
  class SecondOrderTerm;
  
  /// helper class to adopt the correct OperatorTerm based on the term order
  template <int Order>
  struct GetTerm 
  {
    typedef if_then_else< Order == 0, ZeroOrderTerm, 
	          if_then_else< Order == 1, FirstOrderTerm, 
	          if_then_else< Order == 2, SecondOrderTerm,
				                              OperatorTerm > > > type;
  };

  
  /// basic interface for OperatorTerms based on expressions
  template <class Expr, int Order = -1>
  struct GenericOperatorTerm : public GetTerm<Order>::type
  {
    typedef typename GetTerm<Order>::type Super;
    
    /// Expression term stored as copy
    Expr expr;
    
    /// constructor
    /// adds all feSpaces provided by the expression term to auxFeSpaces liste
    GenericOperatorTerm(const Expr& expr_)
      : Super(expr_.getDegree()), expr(expr_) 
    {
      expr.insertFeSpaces(this->auxFeSpaces);
#ifndef NDEBUG
      test_auxFeSpaces(this->auxFeSpaces);
#endif
    }

  private:
    /// \brief Implements OperatorTerm::initImpl().
    /// calls init() for \ref expr
    virtual void initImpl(const ElInfo* elInfo,
                  			  SubAssembler* subAssembler,
                  			  Quadrature* quad) override
    {
      expr.initElement(this, elInfo, subAssembler, quad, NULL);
    }
    
    /// test for only one mesh allowed in expressions
    template <class FeSpaceList>
    void test_auxFeSpaces(FeSpaceList const& auxFeSpaces)
    {
      typedef typename FeSpaceList::const_iterator fe_iter;
      if (auxFeSpaces.size() > 0) {
      	Mesh* mesh0 = (*auxFeSpaces.begin())->getMesh();
      	for (fe_iter it = auxFeSpaces.begin(); it != auxFeSpaces.end(); it++) {
      	  if ((*it)->getMesh() != mesh0) {
      	    ERROR_EXIT("Only one mesh allowed in expression.\n");
      	  }
      	}
      }
    }
  };
  

  template <class Expr>
  struct GenericOperatorTerm<Expr, -1> : public GenericOperatorTerm<Expr, -2>
  {
    typedef GenericOperatorTerm<Expr, -2> Super;
    GenericOperatorTerm(const Expr& expr_) : Super(expr_) { }
    
  private:
    // Implements OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
                  			  const mtl::dense_vector<double>& uhAtQP,
                  			  const mtl::dense_vector<WorldVector<double> >& grdUhAtQP,
                  			  const mtl::dense_vector<WorldMatrix<double> >& D2UhAtQP,
                  			  mtl::dense_vector<double>& result,
                  			  double factor) const override {};
  };
  
} // end namespace AMDiS
