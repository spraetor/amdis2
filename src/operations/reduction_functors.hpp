#pragma once

// std c++ headers
#include <cmath>

// mtl headers
#include <boost/numeric/linear_algebra/identity.hpp>
#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/vector/reduction_functors.hpp>
#include <boost/numeric/mtl/operation/dot.hpp>

// AMDiS headers
#include "operations/assign.hpp"
#include "operations/functors.hpp"

namespace AMDiS
{
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
      using result_type = A;

      template <class Value>
      static void init(Value& value)
      {
        using ::math::zero;
        value= zero(value);
      }

      template <class Value, class Element>
      static void update(Value& value, const Element& x)
      {
        using std::abs;
        value+= abs(x);
      }

      template <class Value>
      static void finish(Value& value, const Value& value2)
      {
        value+= value2;
      }

      template <class Value>
      static constexpr Value post_reduction(Value const& value)
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

      template <class Value>
      static void init(Value& value)
      {
        using ::math::zero;
        value= zero(value);
      }

      template <class Value, class Element>
      static void update(Value& value, Element const& x)
      {
        using mtl::squared_abs;
        value+= squared_abs(x);
      }

      template <class Value>
      static void finish(Value& value, Value const& value2)
      {
        value+= value2;
      }

      // After reduction compute square root
      template <class Value>
      static Value post_reduction(Value const& value)
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
      static Value post_reduction(Value const& value)
      {
        return value;
      }
    };

    /// \cond HIDDEN_SYMBOLS
    template <class A, class B, class ConjOp>
    struct dot_functor_aux
    {
      using result_type = decltype( std::declval<ConjOp>()(std::declval<A>()) * std::declval<B>() );

      template <class Value>
      static void init(Value& value)
      {
        using ::math::zero;
        value= zero(value);
      }

      template <class Value, class Element1, class Element2>
      static void update(Value& value, Element1 const& x, Element2 const& y)
      {
        value+= ConjOp()(x) * y;
      }

      template <class Value>
      static void finish(Value& value, Value const& value2, Value const& value3)
      {
        value+= ConjOp()(value2) * value3;
      }

      template <class Value>
      static constexpr Value post_reduction(Value const& value)
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
    using dot_functor
      = dot_functor_aux<A,B, MTL_VEC::detail::with_conj>;


    /// Binary reduction functor (scalar product)
    /** Same as reduction functors, but \ref update has two arguments:
     *  init: result = 0,
     *  update: result += v_i^T * w_i,
     *  post_reduction: result =result
     **/
    template <class A, class B>
    using dot_real_functor
      = dot_functor_aux<A,B, MTL_VEC::detail::without_conj>;


    template <class ResultType, class InitAssign, class UpdateAssign,
              class PostOp = identity<ResultType>, class FinishAssign = UpdateAssign>
    struct general_unary_reduction_functor
    {
      using result_type = ResultType;

      template <class Value>
      static void init(Value& value)
      {
        InitAssign()(value);
      }

      template <class Value, class Element>
      static void update(Value& value, Element const& x)
      {
        UpdateAssign()(value, x);
      }

      template <class Value>
      static void finish(Value& value, Value const& value2)
      {
        FinishAssign()(value, value2);
      }

      // After reduction compute square root
      template <class Value>
      static constexpr Value post_reduction(Value const& value)
      {
        return PostOp()(value);
      }
    };

    // max(v0, v1, v2, v3,...)
    template <class T>
    using max_reduction_functor
      = general_unary_reduction_functor<T,
        AMDiS::assign::min_value<T>, AMDiS::assign::max<T>>;

    // max(|v0|,|v1|,|v2|,...)
    template <class T>
    using abs_max_reduction_functor
      = general_unary_reduction_functor<T,
        AMDiS::assign::ct_value<T, int, 0>,
        AMDiS::assign::compose<AMDiS::assign::max<T>, 2, abs<T>>>;

    // min(v0, v1, v2, v3, ...)
    template <class T>
    using min_reduction_functor
      = general_unary_reduction_functor<T,
        AMDiS::assign::max_value<T>, AMDiS::assign::min<T>>;

    // min(|v0|,|v1|,|v2|,...)
    template <class T>
    using abs_min_reduction_functor
      = general_unary_reduction_functor<T,
        AMDiS::assign::max_value<T>,
        AMDiS::assign::compose<AMDiS::assign::min<T>, 2, abs<T>>>;

    // v0+v1+v2+v3+...
    template <class T>
    using sum_reduction_functor
      = general_unary_reduction_functor<T,
        AMDiS::assign::ct_value<T, int, 0>, AMDiS::assign::plus<T>>;

    // v0*v1*v2*v3*...
    template <class T>
    using prod_reduction_functor
      = general_unary_reduction_functor<T,
        AMDiS::assign::ct_value<T, int, 1>, AMDiS::assign::multiplies<T>>;

  } // end namespace functors
} // end namespace AMDiS
