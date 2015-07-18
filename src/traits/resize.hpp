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



/** \file resize.hpp */

#ifndef AMDIS_TYPE_TRAITS_RESIZE_HPP
#define AMDIS_TYPE_TRAITS_RESIZE_HPP

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

    
    /// change_dim implementation for AMDiS WorldVectors
    template <typename Value>
    struct resize< WorldVector<Value> > 
    {
      void operator()(WorldVector<Value>& v, size_t r) { 
	TEST_EXIT_DBG(Global::getGeo(WORLD) == r)
	  ("WorldVectors can not be resized!\n");
      }
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
        
    /// change_dim implementation for AMDiS WorldMatrices
    template <typename Value>
    struct resize< WorldMatrix<Value> > 
    {
      void operator()(WorldMatrix<Value>& v, size_t r, size_t c) { 
#ifndef NDEBUG
	size_t dow = static_cast<size_t>(Global::getGeo(WORLD));
	TEST_EXIT_DBG(dow == r && dow == c)
	  ("WorldMatrices can not be resized!\n");
#endif
      }
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
