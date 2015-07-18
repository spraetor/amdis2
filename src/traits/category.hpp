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



/** \file category.hpp */

#ifndef AMDIS_TYPE_TRAITS_CATEGORY_HPP
#define AMDIS_TYPE_TRAITS_CATEGORY_HPP

#include "AMDiS_fwd.h"
#include "tag.hpp"
#include "types.hpp"

#include "DOFMatrix.h"
#include "solver/BlockMTLMatrix.h"

#include <boost/mpl/if.hpp>
#include <boost/mpl/or.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_signed.hpp>
#include <boost/type_traits/is_unsigned.hpp>

#include "boost/numeric/mtl/vector/dense_vector.hpp"
#include "boost/numeric/mtl/matrix/dense2D.hpp"
#include "boost/numeric/mtl/matrix/compressed2D.hpp"
#include "boost/numeric/mtl/matrix/coordinate2D.hpp"
#include "boost/numeric/mtl/matrix/morton_dense.hpp"

namespace AMDiS 
{
  namespace traits 
  {
    // categories
    // _________________________________________________________________________
    
    template<typename T, typename Enabled = void>
    struct category 
    {
      typedef tag::unknown  tag;
      typedef T             value_type;
      typedef size_t        size_type;
    };
    
    // scalars
    template < typename T >
    struct category<T, typename boost::enable_if< typename is_numeric<T>::type >::type >
    {
      typedef tag::scalar tag;
      typedef T           value_type;
//       typedef size_t      size_type;
    };
    
    
    // vectors
    template<typename T>
    struct category<Vector<T> > 
    {
      typedef tag::vector  tag;
      typedef T           value_type;
      typedef int         size_type;
    };
    
    template<typename T, GeoIndex d>
    struct category<FixVec<T, d> > 
    {
      typedef tag::vector  tag;
      typedef T           value_type;
      typedef int         size_type;
    };
    
    template<typename T>
    struct category<DimVec<T> > 
    {
      typedef tag::vector  tag;
      typedef T           value_type;
      typedef int         size_type;
    };
    
    template<typename T>
    struct category<WorldVector<T> > 
    {
      typedef tag::vector  tag;
      typedef T           value_type;
      typedef int         size_type;
    };
    
    template<typename T, typename Allocator>
    struct category<std::vector<T, Allocator> > 
    {
      typedef tag::vector  tag;
      typedef T           value_type;
      typedef size_t      size_type;
    };
    
    template<typename T, typename Parameters>
    struct category<mtl::dense_vector<T, Parameters> > 
    {
      typedef tag::vector  tag;
      typedef typename mtl::dense_vector<T, Parameters>::value_type  value_type;
      typedef typename mtl::dense_vector<T, Parameters>::size_type   size_type;
    };
    
    template <typename T>
    struct category< AMDiS::DOFVector<T> > 
    {
      typedef tag::vector      tag;
      typedef T                value_type;
      typedef DegreeOfFreedom  size_type;
    };
    
    template <typename T>
    struct category< AMDiS::DOFIndexed<T> > 
    {
      typedef tag::vector      tag;
      typedef T                value_type;
      typedef DegreeOfFreedom  size_type;
    };
    
    template <>
    struct category< AMDiS::SystemVector > 
    {
      typedef tag::vector      tag;
      typedef double           value_type;
      typedef int              size_type;
    };
    
    
    // matrices    
    template<typename T>
    struct category<Matrix<T> > 
    {
      typedef tag::matrix  tag;
      typedef T            value_type;
    };
    
    template<typename T>
    struct category<DimMat<T> > 
    {
      typedef tag::matrix  tag;
      typedef T            value_type;
      typedef int         size_type;
    };
    
    template<typename T>
    struct category<WorldMatrix<T> > 
    {
      typedef tag::matrix  tag;
      typedef T            value_type;
      typedef int         size_type;
    };
    
    template<typename T>
    struct category<WorldVector<WorldVector<T> > > 
    {
      typedef tag::matrix  tag;
      typedef T            value_type;
      typedef int          size_type;
    };
    
    template<typename T, typename Parameters>
    struct category<mtl::dense2D<T, Parameters> > 
    {
      typedef tag::matrix  tag;
      typedef typename mtl::dense2D<T, Parameters>::value_type   value_type;
      typedef typename mtl::dense2D<T, Parameters>::size_type    size_type;
    };
    
    template<typename T, typename Parameters>
    struct category<mtl::compressed2D<T, Parameters> > 
    {
      typedef tag::matrix  tag;
      typedef typename mtl::compressed2D<T, Parameters>::value_type   value_type;
      typedef typename mtl::compressed2D<T, Parameters>::size_type    size_type;
    };
    
    template<typename T, typename Parameters>
    struct category<mtl::matrix::coordinate2D<T, Parameters> > 
    {
      typedef tag::matrix  tag;
      typedef typename mtl::coordinate2D<T, Parameters>::value_type   value_type;
      typedef typename mtl::coordinate2D<T, Parameters>::size_type    size_type;
    };
    
    template <typename T, std::size_t BitMask, typename Parameters>
    struct category<mtl::matrix::morton_dense<T, BitMask, Parameters> > 
    {
      typedef tag::matrix  tag;
      typedef typename mtl::morton_dense<T, BitMask, Parameters>::value_type value_type;
      typedef typename mtl::morton_dense<T, BitMask, Parameters>::size_type  size_type;
    };
  
    template <>
    struct category< AMDiS::DOFMatrix > 
    {
      typedef tag::matrix  tag;
      typedef AMDiS::DOFMatrix::value_type  value_type;
      typedef AMDiS::DOFMatrix::size_type   size_type;
    };
  
    template <>
    struct category< AMDiS::BlockMTLMatrix > 
    {
      typedef tag::matrix  tag;
      typedef AMDiS::DOFMatrix::value_type  value_type;
      typedef AMDiS::DOFMatrix::size_type   size_type;
    };
  
    template < class Matrix >
    struct category< AMDiS::SolverMatrix<Matrix> > 
    {
      typedef tag::matrix  tag;
    };
    
  
    template < class T >
    struct category< const T > : category< T > {};
    
  
    // operations on tags
    // _________________________________________________________________________
    
    template< typename T, typename Tag >
    struct has_tag : boost::is_same< typename category<T>::tag, Tag > {};
    
    template<typename T>
    struct is_scalar : has_tag<T, tag::scalar> {};

    template<typename T>
    struct is_vector : has_tag<T, tag::vector> {};
      
    template<typename T>
    struct is_matrix : has_tag<T, tag::matrix> {};
    
    // -------------------------------------------------------------------------
    template<typename T, typename S, typename Enable = void>
    struct is_convertible : boost::is_convertible<T, S> {};
    
    template<typename T, typename S>
    struct is_convertible<T, S, typename boost::enable_if< 
				  typename boost::mpl::and_<
				    typename is_vector<T>::type,
				    typename is_vector<S>::type
				  >::type
				>::type > 
      : boost::is_convertible<typename T::value_type, typename S::value_type> {};
      
    template<typename T, typename S>
    struct is_convertible<T, S, typename boost::enable_if< 
				  typename boost::mpl::and_<
				    typename is_matrix<T>::type,
				    typename is_matrix<S>::type
				  >::type
				>::type > 
      : boost::is_convertible<typename T::value_type, typename S::value_type> {};
    
    
    
  } // end namespace traits

} // end namespace AMDiS

#endif // AMDIS_TYPE_TRAITS_CATEGORY_HPP
