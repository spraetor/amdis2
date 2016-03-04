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



/** \file AbstractFunction.h */

#ifndef AMDIS_ABSTRACTFUNCTION_H
#define AMDIS_ABSTRACTFUNCTION_H

#include "Global.h"

#include <boost/preprocessor/arithmetic/inc.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

namespace AMDiS
{

  /**
   * \ingroup Common
   *
   * \brief
   * Interface for binary functions.
   */
  template<typename ReturnType,
           typename ArgumentType1,
           typename ArgumentType2>
  class BinaryAbstractFunction
  {
  public:
    /// Constructor.
    BinaryAbstractFunction(int degree = 0) :
      degree_(degree)
    {}

    virtual ~BinaryAbstractFunction() {}

    /// Returns \ref degree_.
    inline int getDegree() const
    {
      return degree_;
    }

    /// Deligates the evaluation to overriden method f.
    virtual ReturnType operator()(const ArgumentType1& x,
                                  const ArgumentType2& y) const = 0;

  protected:
    int degree_;
  };

  /**
   * \ingroup Common
   *
   * \brief
   * Interface for tertiary functions.
   */
  template<typename ReturnType,
           typename ArgumentType1,
           typename ArgumentType2,
           typename ArgumentType3>
  class TertiaryAbstractFunction
  {
  public:
    /// Constructor.
    TertiaryAbstractFunction(int degree = 0) :
      degree_(degree)
    {}

    virtual ~TertiaryAbstractFunction() {}

    /// Returns \ref degree_.
    inline int getDegree() const
    {
      return degree_;
    }

    /// function evaluation.
    virtual ReturnType operator()(const ArgumentType1& x,
                                  const ArgumentType2& y,
                                  const ArgumentType3& z) const = 0;

  protected:
    int degree_;
  };

  /**
   * \ingroup Common
   *
   * \brief
   * Interface for quart functions.
   */
  template<typename ReturnType,
           typename ArgumentType1,
           typename ArgumentType2,
           typename ArgumentType3,
           typename ArgumentType4>
  class QuartAbstractFunction
  {
  public:
    /// Constructor.
    QuartAbstractFunction(int degree = 0) :
      degree_(degree)
    {}

    virtual ~QuartAbstractFunction() {}

    /// Returns \ref degree_.
    inline int getDegree() const
    {
      return degree_;
    }

    /// function evaluation.
    virtual ReturnType operator()(const ArgumentType1& x,
                                  const ArgumentType2& y,
                                  const ArgumentType3& z,
                                  const ArgumentType4& u) const = 0;

  protected:
    int degree_;
  };

#if HAS_VARIADIC_TEMPLATES

  // C++11 implementation of AbstractFunction with arbitrary nr of arguments using variadic templates
  /**
   * \ingroup Common
   *
   * \brief
   * An AbstractFunction object represents a function
   * f : ArgumentTypes... -> ReturnType.
   *
   * AbstractFunction is a pure virtual class interface class.
   * To create your own function you have to derive AbstractFunction and
   * overload operator().
   */
  template<typename ReturnType, typename... Args>
  class AbstractFunction
  {
  public:
    AbstractFunction(int degree = 0)
      : degree_(degree)
    {}

    virtual ~AbstractFunction() {}

    /// Returns \ref degree_.
    inline int getDegree() const
    {
      return degree_;
    }

    /// function evaluation.
    virtual ReturnType operator()(Args const& ... args) const = 0;

  protected:
    int degree_;
  };

#else

  template<typename ReturnType, typename ArgumentType>
  class AbstractFunction
  {
  public:
    /// Constructor.
    AbstractFunction(int degree = 0) :
      degree_(degree)
    {}

    virtual ~AbstractFunction() {}

    /// Returns \ref degree_.
    inline int getDegree() const
    {
      return degree_;
    }

    /// Deligates the evaluation to overriden method f.
    virtual ReturnType operator()(const ArgumentType& x) const = 0;

  protected:
    int degree_;
  };

  ///////////////////////////////////////////////////////////////
  // test of AbstractFunction with arbitrary number of arguments

#define ABSTRACT_FUNCTION_MACRO(z, n, _) \
  template< typename ReturnType, \
            BOOST_PP_ENUM_PARAMS_Z(z, n, typename T) \
            > struct AbstractFunction ## n { \
    AbstractFunction ## n (int degree = 0) : \
      degree_(degree) \
    {} \
    virtual ~AbstractFunction ## n () {} \
    inline int getDegree() const \
    { \
      return degree_; \
    } \
    virtual ReturnType operator()(BOOST_PP_ENUM_BINARY_PARAMS_Z(z, n, const T, &t)) const = 0; \
    \
protected: \
    int degree_; \
  }; \


#if 0 // description of the macro
  template<typename ReturnType, BOOST_PP_ENUM_PARAMS(N, typename T)>// expands to typename T0, typename T1, typename T2...
  class CONCAT_STR(AbstractFunction,N)
  {
  public:
    CONCAT_STR(AbstractFunction,N)(int degree = 0) :
      degree_(degree)
    {
    }

    virtual ~CONCAT_STR(AbstractFunction,N)() {}

    /// Returns \ref degree_.
    inline int getDegree() const
    {
      return degree_;
    }

    /// function evaluation.
    virtual ReturnType operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, const T, &t)) const = 0;

  protected:
    int degree_;
  };
#endif

  BOOST_PP_REPEAT_FROM_TO(1, 11, ABSTRACT_FUNCTION_MACRO, nil)

#endif
}

#endif // AMDIS_ABSTRACTFUNCTION_H
