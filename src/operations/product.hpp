/** \file product.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "Traits.h"
#include "traits/num_rows.hpp"
#include "traits/num_cols.hpp"
#include "traits/size.hpp"
#include "operations/functors.hpp"
#include "operations/meta.hpp"

#include "boost/numeric/mtl/operation/dot.hpp"
#include "boost/numeric/mtl/operation/unary_dot.hpp"
#include "boost/numeric/mtl/operation/cross.hpp"

namespace AMDiS
{
#if 0
  // (v1' * v2)
  // ___________________________________________________________________________

  namespace result_of
  {
    template <typename T1, typename T2, typename Tag1, typename Tag2> struct dot {};

    template <typename T1, typename T2>
    struct dot<T1, T2, tag::vector, tag::vector>
      : boost::mpl::identity<typename traits::mult_type<T1, T2>::type> {};
  }

  template <typename T1, typename T2>
  typename result_of::dot<T1, T2, tag::vector, tag::vector>::type
  inline dot_dispatch(const T1& v1, const T2& v2,
                      tag::vector, tag::vector)
  {
    typedef typename traits::category<T1>::value_type value_type1;
    typedef typename traits::category<T2>::value_type value_type2;
    typedef typename traits::category<T1>::size_type size_type1;
    typedef typename traits::category<T2>::size_type size_type2;

    using namespace mtl::vector;
    typedef parameters<mtl::col_major, non_fixed::dimension, false, size_type1> param1;
    typedef parameters<mtl::col_major, non_fixed::dimension, false, size_type2> param2;
    mtl::dense_vector<value_type1, param1> tmp1(AMDiS::size(v1), const_cast<value_type1*>(v1.begin()));
    mtl::dense_vector<value_type2, param2> tmp2(AMDiS::size(v2), const_cast<value_type2*>(v2.begin()));
    return mtl::dot(tmp1, tmp2);
  }

  template <typename T1, typename T2>
  typename boost::lazy_disable_if<typename boost::mpl::or_
  <
  typename traits::is_mtl<T1>::type,
           typename traits::is_mtl<T2>::type
           >::type,
           typename result_of::dot<T1, T2,
           typename traits::category<T1>::tag,
           typename traits::category<T2>::tag
           >
           >::type
           inline dot(const T1& v1, const T2& v2)
  {
    typedef typename traits::category<T1>::tag tag1_;
    typedef typename traits::category<T2>::tag tag2_;
    return dot_dispatch(v1, v2, tag1_(), tag2_());
  }

  template <typename T1, typename T2>
  typename result_of::dot<T1, T2,
           typename traits::category<T1>::tag,
           typename traits::category<T2>::tag
           >::type
           inline inner(const T1& v1, const T2& v2)
  {
    return dot(v1, v2);
  }
#endif

  namespace functors
  {
    /// inner(v1, v2) == v1' * v2
    template<typename T1, typename T2>
    struct Inner : FunctorBase
    {
      typedef typename traits::mult_type<T1, T2>::type result_type;
      int getDegree(int d0, int d1) const
      {
        return d0+d1;
      }

      static result_type eval(const T1& v1, const T2& v2)
      {
        return dot(v1, v2);
      }
      result_type operator()(const T1& v1, const T2& v2) const
      {
        return eval(v1, v2);
      }
    };
  }
#if 0

  // (v' * v)
  // ___________________________________________________________________________

  namespace result_of
  {
    template <typename T, typename Tag> struct unary_dot {};

    template <typename T>
    struct unary_dot<T, tag::vector>
      : boost::mpl::identity<typename traits::mult_type<T,T>::type> {};
  }

  template <typename T>
  typename result_of::unary_dot<T, tag::vector>::type
  inline unary_dot_dispatch(const T& v, tag::vector)
  {
    return v*v;
  }

  template <typename T>
  typename boost::disable_if<typename traits::is_mtl<T>::type,
           typename result_of::unary_dot<T, typename traits::category<T>::tag>::type
           >::type
           inline unary_dot(const T& t)
  {
    return unary_dot_dispatch(t, typename traits::category<T>::tag());
  }
#endif

  namespace functors
  {
    /// unary_dot(v) = v' * v
    template<typename T>
    struct UnaryDot : FunctorBase
    {
      typedef typename traits::mult_type<T, T>::type result_type;
      int getDegree(int d0) const
      {
        return 2*d0;
      }

      static result_type eval(const T& v)
      {
        return unary_dot(v);
      }
      result_type operator()(const T& v) const
      {
        return eval(v);
      }
    };
  }

#if 0
  // (v1 x v2)
  // ___________________________________________________________________________

  namespace result_of
  {
    template <typename T1, typename T2, typename Tag1, typename Tag2, typename Enabled = void> struct cross {};

    template <typename T1, typename T2>
    struct cross<T1, T2, tag::vector, tag::vector,
             typename boost::enable_if<typename traits::is_mtl<T1>::type>::type>
             : mtl::vector::detail::cross_result<T1, T2> {};

    template <template<typename> class Vector, typename T>
    struct cross<Vector<T>, Vector<T>, tag::vector, tag::vector>
      : boost::mpl::identity<Vector<typename traits::mult_type<T, T>::type>> {};
  }

  template <template<typename> class Vector, typename T>
  typename result_of::cross<Vector<T>, Vector<T>, tag::vector, tag::vector>::type
  inline cross_dispatch(const Vector<T>& v1, const Vector<T>& v2, tag::vector, tag::vector)
  {
    typedef typename traits::mult_type<T, T>::type value_type;
    typedef typename traits::category<Vector<T>>::size_type size_type;
    Vector<value_type> result;

    if (size(v1) != 3 || size(v1) != size(v2))
      throw std::runtime_error("cross: inkompatible sizes!");

    for (size_type i = 0; i < size(v1); i++)
    {
      size_type k = (i+1) % 3, l = (i+2) % 3;
      result[i] = v1[k] * v2[l] - v1[l] * v2[k];
    }
    return result;
  }

  template <typename T1, typename T2>
  typename boost::lazy_disable_if<typename boost::mpl::or_
  <
  typename traits::is_mtl<T1>::type,
           typename traits::is_mtl<T2>::type
           >::type,
           typename result_of::cross<T1, T2,
           typename traits::category<T1>::tag,
           typename traits::category<T2>::tag
           >
           >::type
           inline cross(const T1& v1, const T2& v2)
  {
    typedef typename traits::category<T1>::tag tag1_;
    typedef typename traits::category<T2>::tag tag2_;
    return cross_dispatch(v1, v2, tag1_(), tag2_());
  }
  // #endif

  namespace functors
  {
    /// cross(v1, v2) = v1 x v2
    template<typename T1, typename T2>
    struct Cross : FunctorBase
    {
      typedef typename result_of::cross<T1, T2, typename traits::category<T1>::tag,
              typename traits::category<T2>::tag>::type result_type;
      int getDegree(int d0, int d1) const
      {
        return d0+d1;
      }

      static result_type eval(const T1& v1, const T2& v2)
      {
        return cross(v1, v2);
      }
      result_type operator()(const T1& v1, const T2& v2) const
      {
        return eval(v1, v2);
      }
    };
  }

  // #if 0

  // (v1 * v2')    (only for WorldVector)
  // ___________________________________________________________________________
#endif
  namespace result_of
  {
    template <typename T1, typename T2, typename Tag1, typename Tag2> struct outer {};

    template <typename T>
    struct outer<WorldVector<T>, WorldVector<T>, tag::vector, tag::vector>
      : boost::mpl::identity<WorldMatrix<T>> {};

    template <typename Vector, typename Scalar>
    struct outer<Vector, Scalar, tag::vector, tag::scalar>
      : boost::mpl::identity<Vector> {};

    template <typename Scalar, typename Vector>
    struct outer<Scalar, Vector, tag::scalar, tag::vector>
      : boost::mpl::identity<Vector> {};
  }

#if 0
  template <typename T>
  WorldMatrix<T>
  inline outer_dispatch(const WorldVector<T>& v1, const WorldVector<T>& v2, tag::vector, tag::vector)
  {
    WorldMatrix<T> result;
    result.vecProduct(v1, v2);
    return result;
  }

  template <typename Vector, typename Scalar>
  Vector inline outer_dispatch(const Vector& v, const Scalar& s, tag::vector, tag::scalar)
  {
    Vector result(v);
    result *= s;
    return result;
  }

  template <typename Scalar, typename Vector>
  Vector inline outer_dispatch(const Scalar& s, const Vector& v, tag::scalar, tag::vector)
  {
    return outer_dispatch(v, s, tag::vector(), tag::scalar());
  }

  template <typename T1, typename T2>
  typename result_of::outer<T1, T2, typename traits::category<T1>::tag, typename traits::category<T2>::tag>::type
  inline outer(const T1& v1, const T2& v2)
  {
    typedef typename traits::category<T1>::tag tag1_;
    typedef typename traits::category<T2>::tag tag2_;
    return outer_dispatch(v1, v2, tag1_(), tag2_());
  }
#endif

  namespace functors
  {
    /// outer(v1, v2) = v1 * v2'
    template<typename T1, typename T2>
    struct Outer : FunctorBase
    {
      typedef typename result_of::outer<T1, T2, typename traits::category<T1>::tag, typename traits::category<T2>::tag>::type result_type;
      typedef typename traits::category<T1>::value_type Value1;
      typedef typename traits::category<T2>::value_type Value2;
      int getDegree(int d0, int d1) const
      {
        return d0+d1;
      }

      static result_type eval(const WorldVector<Value1>& v1, const WorldVector<Value2>& v2)
      {
        return outer(v1, v2);
      }
      result_type operator()(const WorldVector<Value1>& v1, const WorldVector<Value2>& v2) const
      {
        return eval(v1, v2);
      }
    };
  }
} // end namespace AMDiS
