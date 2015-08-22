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
  /* ------ ASSIGNMENT OF AN EXRESION TO A DOFVECTOR ------------------------ */

  
  /// Assign an expression to a DOFVector
  template <class T, class Term>
  inline requires::Term< void, Term >
  transformDOF(Term term, DOFVector<T>* result);

  
  /// Assign an expression to a DOFVector
  template <class T, class Term>
  inline requires::Term< DOFVector<T>&, Term >
  operator<<(DOFVector<T>& result, Term&& term)
  {
    transformDOF(std::forward<Term>(term), &result);
    return result;
  }

  
  /// Assign an expression to a DOFVector
  template <class T, class F>
  inline requires::CoordsFunctor< DOFVector<T>&, F >
  operator<<(DOFVector<T>& result, F&& f)
  {
    transformDOF(ToTerm<F>::get(std::forward<F>(f)), &result);
    return result;
  }
  
  /* ----- CONVERT EXPRESSION TO STRING ------------------------------------- */

  
  /// Print an expression to an output stream
  template <class Term>
  inline requires::Term< std::ostream&, Term >
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
    transformDOF(toTerm(std::forward<Term>(term)), this);
  }

  template <class T>
  void DOFVector<T>::interpol(DOFVector<T> *v, double factor)
  {
    this->interpol(factor * valueOf(v));
  }


} // end namespace AMDiS

#include "_Expressions.hh"
