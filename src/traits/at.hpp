/** \file at.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "MatrixVector_fwd.h"
#include "category.hpp"

namespace AMDiS 
{
  namespace traits 
  {
    /// General declaration, used to disable unsupported types
    template <typename Collection, class Enable = void>
    struct at {};
      
      
    template <typename T>
    struct at<T, typename boost::enable_if< typename has_tag<T, tag::scalar>::type >::type >
    {
      typedef T   value_type;
      value_type& operator()(T& v, size_t r) { 
	TEST_EXIT_DBG(r == 0)
	  ("Skalars have only 1 component!\n");
	return v;
      }
      value_type const& operator()(const T& v, size_t r) const { 
	TEST_EXIT_DBG(r == 0)
	  ("Skalars have only 1 component!\n");
	return v;
      }
      value_type& operator()(T& v, size_t r, size_t c) { 
	TEST_EXIT_DBG(r == 0 && c == 0)
	  ("Skalars have only 1 row/column!\n");
	return v;
      }
      value_type const& operator()(const T& v, size_t r, size_t c) const { 
	TEST_EXIT_DBG(r == 0 && c == 0)
	  ("Skalars have only 1 row/column!\n");
	return v;
      }
    };
      
    
// == vectors ===
    
    
    template <typename T>
    struct at<T, typename boost::enable_if
    < typename boost::mpl::and_< typename is_vector<T>::type, typename is_mtl<T>::type >::type >::type > 
    {
      typedef typename T::value_type   value_type;
      value_type& operator()(T& v, size_t r) { 
	return v[r];
      }
      value_type const& operator()(const T& v, size_t r) const { 
	return v[r];
      }
    };
    
    /// AMDiS::WorldVector
    template <typename Value>
    struct at< WorldVector<Value> > 
    {
      typedef Value   value_type;
      value_type& operator()(WorldVector<Value>& v, size_t r) { 
	return v[r];
      }
      value_type const& operator()(const WorldVector<Value>& v, size_t r) const { 
	return v[r];
      }
    };

    
// === matrices ===
    
    
    template <typename T>
    struct at<T, typename boost::enable_if
    < typename boost::mpl::and_< typename is_matrix<T>::type, typename is_mtl<T>::type >::type >::type > 
    {
      typedef typename T::value_type   value_type;
      value_type& operator()(T& v, size_t r, size_t c) { 
	return v(r,c);
      }
      value_type const& operator()(const T& v, size_t r, size_t c) { 
	return v(r,c);
      }
    };

    
    template <typename Value>
    struct at< WorldMatrix<Value> > 
    {
      typedef Value   value_type;
      value_type& operator()(WorldMatrix<Value>& v, size_t r, size_t c) { 
	return v[r][c];
      }
      value_type const& operator()(const WorldMatrix<Value>& v, size_t r, size_t c) { 
	return v[r][c];
      }
    };
    
  } // end namespace traits
  
  
  /// at function for vectors
  template <typename Collection>
  typename boost::enable_if< typename traits::is_vector<Collection>::type,
    typename traits::at<Collection>::value_type const&
  >::type
  at(const Collection& c, size_t rows)
  {  
    return traits::at<Collection>()(c, rows);
  }
  
  /// at function for vectors
  template <typename Collection>
  typename boost::enable_if< typename traits::is_vector<Collection>::type,
    typename traits::at<Collection>::value_type&
  >::type
  at(Collection& c, size_t rows)
  {  
    return traits::at<Collection>()(c, rows);
  }
  
  
  /// at function for matrices
  template <typename Collection>
  typename boost::enable_if< typename traits::is_matrix<Collection>::type,
    typename traits::at<Collection>::value_type const&
  >::type
  at(const Collection& c, size_t rows, size_t cols)
  {  
    return traits::at<Collection>()(c, rows, cols);
  }
  
  /// at function for matrices
  template <typename Collection>
  typename boost::enable_if< typename traits::is_matrix<Collection>::type,
    typename traits::at<Collection>::value_type&
  >::type
  at(Collection& c, size_t rows, size_t cols)
  {  
    return traits::at<Collection>()(c, rows, cols);
  }
    
} // end namespace AMDiS

#endif // AMDIS_TYPE_TRAITS_AT_HPP
