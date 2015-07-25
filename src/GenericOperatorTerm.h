#pragma once

#include "OperatorTerm.h"
#include "ZeroOrderTerm.h"
#include "FirstOrderTerm.h"
#include "SeconOrderTerm.h"

namespace AMDiS
{
  
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
    typedef typename GetTerm<Order>::type super;
    
    /// Expression term stored as copy
    Expr expr;
    
    /// constructor
    /// adds all feSpaces provided by the expression term to auxFeSpaces liste
    GenericOperatorTerm(const Expr& expr_)
      : super(term_.getDegree()), expr(expr_) 
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
      expr.init(elInfo, subAssembler, quad, NULL);
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
    typedef GenericOperatorTerm<Expr, -2> super;
    GenericOperatorTerm(const Expr& expr_) : super(expr_) { }
    
  private:
    // Implements OperatorTerm::eval().
    virtual void evalImpl(int nPoints,
			  const mtl::dense_vector<double>& uhAtQP,
			  const mtl::dense_vector<WorldVector<double> >& grdUhAtQP,
			  const mtl::dense_vector<WorldMatrix<double> >& D2UhAtQP,
			  mtl::dense_vector<double>& result,
			  double factor) override {};
  };
  
} // end namespace AMDiS
