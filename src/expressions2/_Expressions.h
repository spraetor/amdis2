/** \file Expressions.h */

#pragma once

#include <ostream>

#include <AMDiS_fwd.h>
#include <OperatorTerm.h>
#include <MatrixVectorOperations.h>
#include <Log.h>

#include "MathOperations.hpp"
#include "TermConcepts.hpp"

namespace AMDiS
{
  /// Integrate an termession over a domain.
  /** If the termession does not contain any DOFVector, a mesh must be given as second argument */
  template <class M>
  inline Value_t<M> integrate(BaseTerm<M> const& term, Mesh* mesh_opt = NULL);


  /* ----- ACCUMULATION OF AND EXPRESSION OVER DOFS ------------------------- */
  

  namespace detail
  {
    /// Accumulate the values of an termession at the Degrees of freedom
    template <class M, class F>
    inline Value_t<M> accumulate(BaseTerm<M> const& term, F f, Value_t<M> value0);
    
  } // end namespace detail

  /// Maximum of an termession at DOFs, using the \ref accumulate function.
  template <class T>
  inline Value_t<T>
  max(BaseTerm<T> const& term)
  {
    Value_t<T> value0 = std::numeric_limits<Value_t<T>>::min();
    value0 = detail::accumulate(term.sub(), functors::max<Value_t<T>>(), value0);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalMax(value0);
#endif
    return value0;
  }


  /// Minimum of an termession at DOFs, using the \ref accumulate function.
  template <class T>
  inline Value_t<T>
  min(BaseTerm<T> const& term)
  {
    Value_t<T> value0 = std::numeric_limits<Value_t<T>>::max();
    value0 = detail::accumulate(term.sub(), functors::min<Value_t<T>>(), value0);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalMin(value0);
#endif
    return value0;
  }


  /// Maximum of  absolute values of an termession at DOFs, using the \ref accumulate function.
  template <class T>
  inline Value_t<T>
  abs_max(BaseTerm<T> const& term)
  {
    Value_t<T> value0 = 0;
    value0 = detail::accumulate(term.sub(), functors::abs_max<Value_t<T>>(), value0);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalMax(value0);
#endif
    return value0;
  }


  /// Minimum of  absolute values of an termession at DOFs, using the \ref accumulate function.
  template <class T>
  inline Value_t<T>
  abs_min(BaseTerm<T> const& term)
  {
    Value_t<T> value0 = std::numeric_limits<Value_t<T>>::max();
    value0 = detail::accumulate(term.sub(), functors::abs_min<Value_t<T>>(), value0);

#ifdef HAVE_PARALLEL_DOMAIN_AMDIS
    Parallel::mpi::globalMin(value0);
#endif
    return value0;
  }
  
  /* ------ ASSIGNMENT OF AN EXRESION TO A DOFVECTOR ------------------------ */

  namespace detail
  {
    /// Assign an term to a DOFVector
    template <class T, class M>
    inline void transformDOF(BaseTerm<M> const& term, DOFVector<T>* result);
    
  } // end namespace detail

  /// Assign an term to a DOFVector
  template <class T, class M>
  inline DOFVector<T>&
  operator<<(DOFVector<T>& result, BaseTerm<M> const& term)
  {
    detail::transformDOF(term.sub(), &result);
    return result;
  }


  /// Assign a coords-functor to a DOFVector
  template <class T, class F>
  inline requires::CoordsFunctor<DOFVector<T>&, F>
  operator<<(DOFVector<T>& result, F&& f)
  {
    detail::transformDOF(toTerm(std::forward<F>(f)), &result);
    return result;
  }


  /// Assign a scalar value to a DOFVector
  template <class T, class S>
  inline Requires_t<traits::is_arithmetic<S>, DOFVector<T>> &
  operator<<(DOFVector<T>& result, S scalar)
  {
    detail::transformDOF(toTerm(scalar), &result);
    return result;
  }

  /* ----- CONVERT EXPRESSION TO STRING ------------------------------------- */


  /// Print an termession to an output stream
  template <class Term>
  inline requires::Term<std::ostream&, Term>
  operator<<(std::ostream& result, const Term& term)
  {
    result << term.str();
    return result;
  }


  /* ----- IMPLEMENTATION OF DOFVECTOR::INTERPOL METHODS --------------------- */


  template <class T>
  template <class Term>
  inline void DOFVector<T>::interpol(Term&& term)
  {
    detail::transformDOF(toTerm(std::forward<Term>(term)), this);
  }

  template <class T>
  void DOFVector<T>::interpol(DOFVector<T>* v, double factor)
  {
    this->interpol(factor * valueOf(v));
  }


} // end namespace AMDiS

#include "_Expressions.hh"
