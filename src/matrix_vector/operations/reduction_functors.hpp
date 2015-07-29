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



/** \file reduction_functors.h */

#pragma once

#include <cmath>

#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/vector/reduction_functors.hpp>
#include <boost/numeric/mtl/operation/dot.hpp>
#include <boost/integer_traits.hpp>

#include "traits/mult_type.hpp"

#include "operations/functors.hpp"
#include "operations/assign.hpp"

namespace AMDiS {

  namespace functors
  {

    // unary reduction functors: import from mtl
    using MTL_VEC::infinity_norm_functor;
    using MTL_VEC::sum_functor;
    using MTL_VEC::product_functor;
    using MTL_VEC::max_functor;
    using MTL_VEC::min_functor;
    
    /// Reduction functor to calculate the ||v||_1
    /** Defines \ref init, \ref update and \ref post_reduction
     *  so that ||v||_1 := |v_0| + |v_1| + ...:
     *  init: result = 0, 
     *  update: result += |v_i|, 
     *  post_reduction: result = result
     **/
    template <class A>
    struct one_norm_functor
    {
      typedef A result_type;
      
      template <typename Value>
      static inline void init(Value& value)
      {
	using math::zero;
	value= zero(value);
      }

      template <typename Value, typename Element>
      static inline void update(Value& value, const Element& x)
      {    
	using std::abs;
	value+= abs(x);
      }

      template <typename Value>
      static inline void finish(Value& value, const Value& value2)
      {
	value+= value2;
      }

      template <typename Value>
      static inline Value post_reduction(const Value& value)
      {
	return value;
      }
    };


    /// Reduction functor to calculate the |v|_2
    /** Defines \ref init, \ref update and \ref post_reduction
     *  so that ||v||_2 := sqrt(|v_0|^2 + |v_1|^2 + ...):
     *  init: result = 0, 
     *  update: result += |v_i|^2, 
     *  post_reduction: result = sqrt(result)
     **/
    template <class A>
    struct two_norm_functor
    {
      typedef A result_type;
      
      template <typename Value>
      static inline void init(Value& value)
      {
	using math::zero;
	value= zero(value);
      }

      template <typename Value, typename Element>
      static inline void update(Value& value, const Element& x)
      {    
	using mtl::squared_abs;
	value+= squared_abs(x);
      }

      template <typename Value>
      static inline void finish(Value& value, const Value& value2)
      {
	value+= value2;
      }

      // After reduction compute square root
      template <typename Value>
      static inline Value post_reduction(const Value& value)
      {
	using std::sqrt;
	return sqrt(value);
      }
    };

    /// Same as \ref two_norm_functor without the root at the end:
    /**
     *  post_reduction: result = result
     **/
    template <class A>
    struct unary_dot_functor
	: two_norm_functor<A>
    {
      template <typename Value>
      static inline Value post_reduction(const Value& value)
      {
	return value;
      }
    };
    
    /// \cond HIDDEN_SYMBOLS
    template <class A, class B, class ConjOp>
    struct dot_functor_aux
    {
      typedef typename if_< boost::is_same<ConjOp, MTL_VEC::detail::with_conj>,
			    typename mtl::sfunctor::conj<A>::result_type,
			    A >::type A_;
			    
      typedef typename traits::mult_type<A_, B>::type result_type;
      
      template <typename Value>
      static inline void init(Value& value)
      {
	using math::zero;
	value= zero(value);
      }

      template <typename Value, typename Element1, typename Element2>
      static inline void update(Value& value, const Element1& x, const Element2& y)
      {    
	value+= ConjOp()(x) * y;
      }

      template <typename Value>
      static inline void finish(Value& value, const Value& value2, const Value& value3)
      {
	value+= ConjOp()(value2) * value3;
      }

      template <typename Value>
      static inline Value post_reduction(const Value& value)
      {
	return value;
      }
    };
    /// \endcond
    
    
    /// Binary reduction functor (scalar product)
    /** Same as reduction functors, but \ref update has two arguments:
     *  init: result = 0, 
     *  update: result += v_i^H * w_i, 
     *  post_reduction: result =result
     **/
    template <class A, class B>
    struct dot_functor 
	: dot_functor_aux<A,B, MTL_VEC::detail::with_conj> {};
	
    
    /// Binary reduction functor (scalar product)
    /** Same as reduction functors, but \ref update has two arguments:
     *  init: result = 0, 
     *  update: result += v_i^T * w_i, 
     *  post_reduction: result =result
     **/
    template <class A, class B>
    struct dot_real_functor 
	: dot_functor_aux<A,B, MTL_VEC::detail::without_conj> {};

	
    template <class ResultType, class InitAssign, class UpdateAssign, 
	      class PostOp = identity<ResultType>, class FinishAssign = UpdateAssign>
    struct general_unary_reduction_functor
    {
      typedef ResultType result_type;
      
      template <class Value>
      static inline void init(Value& value)
      {
	InitAssign op;
	op(value);
      }

      template <class Value, class Element>
      static inline void update(Value& value, const Element& x)
      {   
	UpdateAssign op;
	op(value, x);
      }

      template <class Value>
      static inline void finish(Value& value, const Value& value2)
      {
	FinishAssign op;
	op(value, value2);
      }

      // After reduction compute square root
      template <class Value>
      static inline Value post_reduction(const Value& value)
      {
	PostOp op;
	return op(value);
      }
    };
	
    // max(v0, v1, v2, v3,...)
    template <class T>
    struct max_reduction_functor
	: general_unary_reduction_functor<T, 
	    AMDiS::assign::min_value<T>, AMDiS::assign::max<T> > {};
	
    // max(|v0|,|v1|,|v2|,...)
    template <class T>
    struct abs_max_reduction_functor
	: general_unary_reduction_functor<T, 
	    AMDiS::assign::ct_value<T, int, 0>, 
	    compose<AMDiS::assign::max<T>, 2, abs<T> > > {};
	
    // min(v0, v1, v2, v3, ...)
    template <class T>
    struct min_reduction_functor
	: general_unary_reduction_functor<T, 
	    AMDiS::assign::max_value<T>, AMDiS::assign::min<T> > {};
	
    // min(|v0|,|v1|,|v2|,...)
    template <class T>
    struct abs_min_reduction_functor
	: general_unary_reduction_functor<T, 
	    AMDiS::assign::max_value<T>, 
	    compose<AMDiS::assign::min<T>, 2, abs<T> > > {};
	
    // v0+v1+v2+v3+...
    template <class T>
    struct sum_reduction_functor
	: general_unary_reduction_functor<T, 
	    AMDiS::assign::ct_value<T, int, 0>, AMDiS::assign::plus<T> > {};
	
    // v0*v1*v2*v3*...
    template <class T>
    struct prod_reduction_functor
	: general_unary_reduction_functor<T, 
	    AMDiS::assign::ct_value<T, int, 1>, AMDiS::assign::multiplies<T> > {};
	
  } // end namespace functors
} // end namespace AMDiS
