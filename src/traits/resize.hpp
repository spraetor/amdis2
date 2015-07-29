/** \file resize.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "category.hpp"

namespace AMDiS 
{
  namespace traits 
  {
  
    /// General declaration, used to disable unsupported types
    template <typename Collection, class Enable = void>
    struct resize {};
      
      
    template <typename T>
    struct resize<T, typename boost::enable_if< typename has_tag<T, tag::scalar>::type >::type >
    {
	void operator()(T& v, size_t r) const { 
	TEST_EXIT_DBG(r == 1)
	  ("Skalars can not be resized!\n");
	}
    };
      
// == vectors ===
    
    template <typename T>
    struct resize<T, typename boost::enable_if
    < typename boost::mpl::and_< typename is_vector<T>::type, typename is_mtl<T>::type >::type >::type > 
    {
	void operator()(T& v, size_t r) const { v.change_dim(r); }
    };

    
    /// change_dim implementation for STL vectors
    template <typename Value>
    struct resize< std::vector<Value> > 
    {
      void operator()(std::vector<Value>& v, size_t r) { 
	v.resize(r);
      }
    };
        
  
// === matrices ===
    
    template <typename T>
    struct resize<T, typename boost::enable_if
    < typename boost::mpl::and_< typename is_matrix<T>::type, typename is_mtl<T>::type >::type >::type > 
    {
	void operator()(T& m, size_t r, size_t c) { return m.change_dim(r,c); }
    };
        
  } // end namespace traits
  
    
  /// resize function for vectors
  template <typename Collection>
  typename boost::enable_if< typename traits::is_vector<Collection>::type >::type
  resize(Collection& c, size_t rows)
  {  
    traits::resize<Collection>()(c, rows);
  }
  
  /// resize function for matrices
  template <typename Collection>
  typename boost::enable_if< typename traits::is_matrix<Collection>::type >::type
  resize(Collection& c, size_t rows, size_t cols)
  {  
    traits::resize<Collection>()(c, rows, cols);
  }
    
} // end namespace AMDiS

#endif // AMDIS_TYPE_TRAITS_RESIZE_HPP
