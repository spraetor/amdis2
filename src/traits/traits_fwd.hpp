/** \file traits_fwd.hpp */

#pragma once

#include <traits/tag.hpp>

#include <boost/numeric/mtl/operation/num_cols.hpp>
#include <boost/numeric/mtl/operation/num_rows.hpp>
#include <boost/numeric/mtl/operation/size.hpp>

namespace AMDiS 
{
  namespace traits 
  {
    /// categories
    template <class T, class Enable = void>
    struct category 
    {
      typedef tag::unknown  tag;
      typedef T             value_type;
      typedef size_t        size_type;
    };
    
    /// General declaration, used to disable unsupported types
    template <class T, class Enable = void>
    struct at {};
    
      
    template <class T, class Enable = void>
    struct num_cols 
      : ::mtl::traits::num_cols<T> {}; // import mtl-operation by default
      
  
    template <class T, class Enable = void>
    struct num_rows 
      : ::mtl::traits::num_rows<T> {}; // import mtl-operation by default
      
    /// General declaration, used to disable unsupported types
    template <class T, class Enable = void>
    struct resize {};
    
    template <class T, class Enable = void>
    struct size 
      : ::mtl::traits::size<T> {}; // import mtl-operation by default
      
    // operations on tags
    // _________________________________________________________________________
    
    template <class T, class Tag>
    using has_tag = std::is_same< typename category<T>::tag, Tag >;
    
    template <class T>
    struct is_scalar : has_tag<T, tag::scalar> {};

    template <class T>
    struct is_vector : has_tag<T, tag::vector> {};
      
    template <class T>
    struct is_matrix : has_tag<T, tag::matrix> {};
    
  } // end namespace traits
  
  
  
  
  template <class T>
  size_t inline num_cols(const T& t)
  {
      return traits::num_cols<T>()(t);
  }
  
  template <class T>
  size_t inline num_rows(const T& t)
  {
      return traits::num_rows<T>()(t);
  }
  
  template <class T>
  size_t inline size(const T& t)
  {
      return traits::size<T>()(t);
  }
  
} // end namespace AMDiS
