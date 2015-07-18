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



/** \file Expressions.h */

#ifndef AMDIS_EXPRESSIONS_BASE_H
#define AMDIS_EXPRESSIONS_BASE_H

#include "AMDiS_fwd.h"
#include "OperatorTerm.h"
#include "Functors.h"
#include "MatrixVectorOperations.h"

#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include "expressions/LazyOperatorTerm.h"
#include "expressions/expressions.h"


/** \brief Expressions provide an easy way of automated generation of
 * 'arbitrary' operator-terms out of some elementary operations, by using a 
 * recursive definition of the term. All necessary data will be initialized 
 * when an expression as part of the term uses this data.
 * Since no virtual functions, like in the AbstractFunction classes, are used
 * the overhead of a vtable is removed.
 * 
 * Usage:
 * addZOT(Operator, Term) 
 *  ... add a zeroOrderTerm to Operator, i.e. (Term(x) * u, v)
 * 
 * addFOT(Operator, Term, FirstOrderType) 
 *  ... add a firstOrderTerm to Operator, (if Term::value_type = double) 
 *      i.e. (Term(x) * 1 * grad(u), v), rsp. (Term(x) * 1 * u, grad(v))
 * addFOT(Operator, Term, FirstOrderType)
 *  ... add a firstOrderTerm to Operator, (if Term::value_type = WorldVector)
 *      i.e. (Term(x) * b * grad(u), v), rsp. (Term(x) * u, grad(v))
 * addFOT<I>(Operator, Term, FirstOrderType)
 *  ... add a firstOrderTerm to Operator, 
 *      i.e. (Term(x) * e_I * grad(u), v), rsp. (Term(x) * e_I * u, grad(v))
 * 
 * addSOT(Operator, Term)
 *  ... add a secondOrderTerm to Operator, i.e. (Term(x) * grad(u), grad(v))
 * addSOT<I,J>(Operator, Term)
 *  ... add a secondOrderTerm to Operator, i.e. (E_IJ * Term(x) * grad(u), grad(v))
 * 
 * where Operator is eather a pointer or reference, FirstOrderType in {GRD_PHI, GRD_PSI}
 * and Term a componation of elementary terms by + - * /
 *  - constant(value) / value  ... a constant value
 *  - valueOf(DOFVector)       ... values of a DOFVector at QP
 *  - gradientOf(DOFVector)    ... gradient of a DOFVector at QP
 *  - derivative<I>(DOFVector) ... I'th partial derivative
 *  - X()                      ... coordinate at quadrature points
 *  - pow<I>(Term)             ... I'th power of a term
 *  - sqrt(Term)               ... square root of a term
 *  - Exp(Term)                ... exponential function of a term
 *  - function_<F>(Term)       ... evaluates F()(Term(iq))
 *  - function_(F f, Term)     ... evaluates f(Term(iq))
 * 
 * 
 * with F a functor that implements
 *   typedef (...) result_type;
 *   int getDegree(int d0);
 *   result_type operator()(const T0& v0) const;
 * 
 * respective
 *   int getDegree(int d0, int d1);
 *   result_type operator()(const T0& v0, const T1& v1) const;
 * 
 * respective
 *   int getDegree(int d0, int d1, int d2);
 *   result_type operator()(const T0& v0, const T1& v1, const T2& v2) const;
 * 
 * where the d0, d1, d2 give the polynomial degrees of the v0, v1, v2 terms.
 * */


namespace AMDiS {

  /// Integrate an expression over a domain.
  /** If the expression does not contain any DOFVector, a mesh must be given as second argument */
  template <class Expr>
  inline typename boost::enable_if<traits::is_expr<Expr>, typename Expr::value_type>::type
  integrate(Expr expr, Mesh* mesh_ = NULL);

  
  /* ----- ACCUMULATION OF AND EXPRESSION OVER DOFS ------------------------- */

  
  /// Accumulate the values of an expression at the Degrees of freedom
  template <class Expr, class Functor>
  inline typename boost::enable_if<traits::is_expr<Expr>, typename Expr::value_type>::type
  accumulate(Expr expr, Functor f, typename Expr::value_type value0);

  
  /// Maximum of an expression at DOFs, using the \ref accumulate function.
  template <class Expr>
  inline typename boost::enable_if<traits::is_expr<Expr>, typename Expr::value_type>::type
  max(Expr expr)
  {
    typename Expr::value_type value0 = std::numeric_limits<typename Expr::value_type>::min();
    value0 = accumulate(expr, functors::max<typename Expr::value_type>(), value0);
    
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalMax(value0);
#endif
    return value0;
  }
  

  /// Minimum of an expression at DOFs, using the \ref accumulate function.
  template <class Expr>
  inline typename boost::enable_if<traits::is_expr<Expr>, typename Expr::value_type>::type
  min(Expr expr)
  {
    typename Expr::value_type value0 = std::numeric_limits<typename Expr::value_type>::max();
    value0 = accumulate(expr, functors::min<typename Expr::value_type>(), value0);
    
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalMin(value0);
#endif
    return value0;
  }
  

  /// Maximum of  absolute values of an expression at DOFs, using the \ref accumulate function.
  template <class Expr>
  inline typename boost::enable_if<traits::is_expr<Expr>, typename Expr::value_type>::type
  abs_max(Expr expr)
  {
    typename Expr::value_type value0 = 0;
    value0 = accumulate(expr, functors::abs_max<typename Expr::value_type>(), value0);
    
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalMax(value0);
#endif
    return value0;
  }

  
  /// Minimum of  absolute values of an expression at DOFs, using the \ref accumulate function.
  template <class Expr>
  inline typename boost::enable_if<traits::is_expr<Expr>, typename Expr::value_type>::type
  abs_min(Expr expr)
  {
    typename Expr::value_type value0 = std::numeric_limits<typename Expr::value_type>::max();
    value0 = accumulate(expr, functors::abs_min<typename Expr::value_type>(), value0);
    
#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalMin(value0);
#endif
    return value0;
  }
  
  
  /* ------ ASSIGNMENT OF AN EXRESION TO A DOFVECTOR ------------------------ */

  
  /// Assign an expression to a DOFVector
  template <class T, class Expr>
  inline typename boost::enable_if<
    typename and_<traits::is_expr<Expr>, 
		  traits::is_convertible<typename Expr::value_type, T>
		  >::type
    >::type
  transformDOF(Expr expr, DOFVector<T>* result);

  
  /// \brief Assign an expression to a DOFVector 
  /// (using multi-mesh if \ref expr and \ref result vector are on different meshes)
  template <class T, class Expr>
  inline typename boost::enable_if<
    typename and_<traits::is_expr<Expr>, 
		  traits::is_convertible<typename Expr::value_type, T>
		  >::type
    >::type
  transformDOF_mm(Expr expr, DOFVector<T>* result);

  
  /// Assign an expression to a DOFVector
  template <class T, class Expr>
  typename boost::enable_if<
    typename and_<traits::is_expr<Expr>, 
		  traits::is_convertible<typename Expr::value_type, T>
		  >::type,
    DOFVector<T>& >::type
  operator<<(DOFVector<T>& result, const Expr& expr)
  {
    transformDOF(expr, &result);
    return result;
  }

  
  /* ----- CONVERT EXPRESSION TO STRING ------------------------------------- */

  
  /// Print an expression to an output stream
  template <class Expr>
  typename boost::enable_if<typename traits::is_expr<Expr>::type, 
			    std::ostream& >::type
  operator<<(std::ostream& result, const Expr& expr)
  {
    result << expr.str();
    return result;
  }

} // end namespace AMDiS

#include "Expressions.hh"

#endif // AMDIS_EXPRESSIONS_BASE_H