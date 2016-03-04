#pragma once

// std c++ headers
#include <ostream>

// AMDiS includes
#include "AMDiS_fwd.hpp"
#include "Log.hpp"
#include "MatrixVectorOperations.hpp"
#include "OperatorTerm.hpp"

#include "expressions/MathOperations.hpp"
#include "expressions/TermConcepts.hpp"

namespace AMDiS
{
  /// Integrate an term over a domain.
  /** If the term does not contain any DOFVector, a mesh must be given as second argument */
  template <class M>
  inline Value_t<M> integrate(BaseTerm<M> const& term, Mesh* mesh_opt = NULL);

  
  /// Integrate an term over the boundary of a domain.
  /** If the term does not contain any DOFVector, a mesh must be given as second 
   *  argument to the boundary_wrapper. */
  template <class M>
  inline Value_t<M> integrate(BaseTerm<M> const& term, BoundaryWrapper b);
  

  // ________ ACCUMULATION OF AND EXPRESSION OVER DOFS _________________________
  

  namespace detail
  {
    /// Accumulate the values of an term at the Degrees of freedom
    template <class M, class F>
    inline Value_t<M> accumulate(BaseTerm<M> const& term, F f, Value_t<M> value0);
    
  } // end namespace detail

  
  /// Maximum of an term at DOFs, using the \ref accumulate function.
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


  /// Minimum of an term at DOFs, using the \ref accumulate function.
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


  /// Maximum of  absolute values of an term at DOFs, using the \ref accumulate function.
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


  /// Minimum of  absolute values of an term at DOFs, using the \ref accumulate function.
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
  
  
  // _________ ASSIGNMENT OF AN EXRESION TO A DOFVECTOR ________________________

  
  namespace detail
  {
    /// Assign an term to a DOFVector
    template <class T, class M>
    inline void transformDOF(BaseTerm<M> const& term, DOFVector<T>* result);
    
  } // end namespace detail

  /// Assign an term to a DOFVector
  template <class T, class M>
  DOFVector<T>& operator<<(DOFVector<T>& result, BaseTerm<M> const& term)
  {
    detail::transformDOF(term.sub(), &result);
    return result;
  }


  /// Assign a coords-functor to a DOFVector
  template <class T, class F>
    requires::CoordsFunctor<DOFVector<T>&, F>
  operator<<(DOFVector<T>& result, F&& f)
  {
    detail::transformDOF(toTerm(std::forward<F>(f)), &result);
    return result;
  }


  /// Assign a scalar value to a DOFVector
  template <class T, class S>
    Requires_t<concepts::Arithmetic<S>, DOFVector<T>> &
  operator<<(DOFVector<T>& result, S scalar)
  {
    detail::transformDOF(toTerm(scalar), &result);
    return result;
  }

  
  // _________ CONVERT EXPRESSION TO STRING ____________________________________


  /// Print an term to an output stream
  template <class M>
  inline std::ostream&
  operator<<(std::ostream& result, BaseTerm<M> const& term)
  {
    result << term.sub().str();
    return result;
  }


  // _________ IMPLEMENTATION OF DOFVECTOR::INTERPOL METHODS ___________________


  template <class T>
  template <class Term>
  void DOFVector<T>::interpol(Term const& term)
  {
    detail::transformDOF(toTerm(term), this);
  }

  template <class T>
  void DOFVector<T>::interpol(DOFVector<T>* v, double factor)
  {
    this->interpol(factor * valueOf(v));
  }


} // end namespace AMDiS

#include "Expressions.hh"
