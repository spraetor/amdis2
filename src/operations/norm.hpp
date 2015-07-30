/** \file norm.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "Traits.h"
#include "traits/num_rows.hpp"
#include "traits/num_cols.hpp"
#include "traits/size.hpp"
#include "boost/numeric/mtl/operation/two_norm.hpp"
#include "boost/numeric/mtl/operation/one_norm.hpp"
#include "boost/numeric/mtl/vector/dense_vector.hpp"
#include "boost/numeric/mtl/vector/parameter.hpp"

#include "operations/functors.hpp"
#include "operations/generic_loops.hpp"

namespace AMDiS 
{  
  // 2-norm
  // ___________________________________________________________________________
  
  namespace result_of
  {
    template <typename T, typename Tag> struct two_norm {};
    
    template <typename T>
    struct two_norm<T, tag::vector> : boost::mpl::identity<typename traits::category<T>::value_type> {};
  }
  
  template <typename T>
  typename result_of::two_norm<T, tag::vector>::type
  inline two_norm_dispatch(const T& v, tag::vector)
  {
    typedef typename traits::category<T>::value_type value_type;
    typedef typename traits::category<T>::size_type size_type;
    
    using namespace mtl::vector;
    typedef parameters<mtl::col_major, non_fixed::dimension, false, size_type> param;
    mtl::dense_vector<value_type, param> tmp(AMDiS::size(v), const_cast<value_type*>(v.begin()));
    return mtl::two_norm(tmp);
  }
  
  template <typename T>
  typename disable_if< traits::is_mtl<T>,
    typename result_of::two_norm<T, typename traits::category<T>::tag>::type
  >::type
  inline two_norm(const T& t)
  {
    return two_norm_dispatch(t, typename traits::category<T>::tag());
  }
  
  
  namespace functors 
  {
    /// two_norm(v) == ||v||_2 == sqrt(v' * v)
    template<typename T>
    struct TwoNorm : FunctorBase
    {
      typedef typename result_of::two_norm<T, typename traits::category<T>::tag>::type result_type;      
      int getDegree(int d0) const { return 2*d0; }
      
      static result_type eval(const T& v) { return AMDiS::two_norm(v); }
      result_type operator()(const T& v) const { return eval(v); }
    };
  }
  
    
  // 1-norm
  // ___________________________________________________________________________
    
  namespace result_of
  {
    template <typename T, typename Tag> struct one_norm {};
    
    template <typename T>
    struct one_norm<T, tag::vector> : boost::mpl::identity<typename traits::category<T>::value_type> {};
    template <typename T>
    struct one_norm<T, tag::matrix> : boost::mpl::identity<typename traits::category<T>::value_type> {};
  }
  
  template <typename T>
  typename result_of::one_norm<T, tag::vector>::type
  inline one_norm_dispatch(const T& v, tag::vector) // for vector types
  {
    typedef typename traits::category<T>::value_type value_type;
    typedef typename traits::category<T>::size_type size_type;
    
    using namespace mtl::vector;
    typedef parameters<mtl::col_major, non_fixed::dimension, false, size_type> param;
    mtl::dense_vector<value_type, param> tmp(AMDiS::size(v), const_cast<value_type*>(v.begin()));
    return mtl::vector::one_norm(tmp);
  }
  
  template <typename T>
  typename result_of::one_norm<T, tag::matrix>::type
  inline one_norm_dispatch(const T& m, tag::matrix) // for matrix types
  {
    typedef typename traits::category<T>::value_type value_type;
    value_type result; 
    nullify(result);
    
    typedef typename traits::category<T>::size_type size_type;
    for (size_type j = 0; j < num_cols(m); ++j) {
      value_type asum; nullify(asum);
      for (size_type i = 0; i < num_rows(m); ++i)
        asum += std::abs(at(m,i,j));
      result = std::max(result, asum);
    }
    return result;
  }
  
  template <typename T>
  typename disable_if< traits::is_mtl<T>,
    typename result_of::one_norm<T, typename traits::category<T>::tag>::type
  >::type
  inline one_norm(const T& t)
  {
    return one_norm_dispatch(t, typename traits::category<T>::tag());
  }
  
  
  namespace functors 
  {
    /// one_norm(v) == ||v||_1 == max_i(|v_i|)
    template<typename T>
    struct OneNorm : FunctorBase
    {
      typedef typename result_of::one_norm<T, typename traits::category<T>::tag>::type result_type;      
      int getDegree(int d0) const { return d0; }
      
      static result_type eval(const T& v) { return AMDiS::one_norm(v); }
      result_type operator()(const T& v) const { return eval(v); }
    };
  }
    
    
  // p-norm
  // ___________________________________________________________________________
  
  namespace result_of
  {
    template <int P, typename T, typename Tag> struct p_norm {};
    
    template <int P, typename T>
    struct p_norm<P, T, tag::vector> : boost::mpl::identity<typename traits::category<T>::value_type> {};
  }
  
  template <int P, typename T>
  typename result_of::p_norm<P, T, tag::vector>::type
  inline p_norm_dispatch(const T& v, tag::vector)
  {
    typedef typename result_of::p_norm<P, T, tag::vector>::type value_type;
    typedef typename traits::category<T>::size_type size_type;
    static const int BLOCK_SIZE = 8;
    
    value_type result;
    nullify(result);
    
    size_type i = 0;
    for (; i+BLOCK_SIZE < size(v); i+=BLOCK_SIZE)
      result = meta::FOR<0, BLOCK_SIZE>::accumulate(v, result, functors::pow<P, value_type>(), std::plus<value_type>(), i);
    for (; i < size(v); i++)
      result += functors::pow<P, value_type>::eval(v[i]);
    
    return functors::root<P, value_type>::eval(result);
  }
  
  template <int P, typename T>
  typename result_of::p_norm<P, T, typename traits::category<T>::tag>::type
  inline p_norm(const T& t)
  {
    return p_norm_dispatch<P>(t, typename traits::category<T>::tag());
  }
  
  
  namespace functors 
  {
    ///
    template<int p, typename T>
    struct PNorm : FunctorBase
    {
      typedef typename result_of::p_norm<p, T, typename traits::category<T>::tag>::type result_type;      
      int getDegree(int d0) const { return p*d0; }
      
      static result_type eval(const T& v) { return AMDiS::p_norm<p>(v); }
      result_type operator()(const T& v) const { return eval(v); }
    };
  }

} // end namespace AMDiS
