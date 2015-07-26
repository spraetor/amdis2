/** \file category.hpp */

#pragma once

#include <type_traits>

#include "AMDiS_fwd.h"
#include "tag.hpp"
#include "types.hpp"

#include "DOFMatrix.h"
#include "solver/BlockMTLMatrix.h"

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
    
    template <class T, class Enabled = void>
    struct category 
    {
      typedef tag::unknown  tag;
      typedef T             value_type;
      typedef size_t        size_type;
    };
    
    // scalars
    template <class T>
    struct category<T, typename enable_if< std::is_arithmetic<T> >::type >
    {
      typedef tag::scalar tag;
      typedef T           value_type;
//       typedef size_t      size_type;
    };
    
    
    // vectors
    template <class T>
    struct category<Vector<T> > 
    {
      typedef tag::vector  tag;
      typedef T           value_type;
      typedef int         size_type;
    };
    
    template <class T, GeoIndex d>
    struct category<FixVec<T, d> > 
    {
      typedef tag::vector  tag;
      typedef T           value_type;
      typedef int         size_type;
    };
    
    template <class T>
    struct category<DimVec<T> > 
    {
      typedef tag::vector  tag;
      typedef T           value_type;
      typedef int         size_type;
    };
    
    template <class T>
    struct category<WorldVector<T> > 
    {
      typedef tag::vector  tag;
      typedef T           value_type;
      typedef int         size_type;
    };
    
    template <class T, class Allocator>
    struct category<std::vector<T, Allocator> > 
    {
      typedef tag::vector  tag;
      typedef T           value_type;
      typedef size_t      size_type;
    };
    
    template <class T, class Parameters>
    struct category<mtl::dense_vector<T, Parameters> > 
    {
      typedef tag::vector  tag;
      typedef typename mtl::dense_vector<T, Parameters>::value_type  value_type;
      typedef typename mtl::dense_vector<T, Parameters>::size_type   size_type;
    };
    
    template <class T>
    struct category< AMDiS::DOFVector<T> > 
    {
      typedef tag::vector      tag;
      typedef T                value_type;
      typedef DegreeOfFreedom  size_type;
    };
    
    template <class T>
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
    template <class T>
    struct category<Matrix<T> > 
    {
      typedef tag::matrix  tag;
      typedef T            value_type;
    };
    
    template <class T>
    struct category<DimMat<T> > 
    {
      typedef tag::matrix  tag;
      typedef T            value_type;
      typedef int         size_type;
    };
    
    template <class T>
    struct category<WorldMatrix<T> > 
    {
      typedef tag::matrix  tag;
      typedef T            value_type;
      typedef int         size_type;
    };
    
    template <class T>
    struct category<WorldVector<WorldVector<T> > > 
    {
      typedef tag::matrix  tag;
      typedef T            value_type;
      typedef int          size_type;
    };
    
    template <class T, class Parameters>
    struct category<mtl::dense2D<T, Parameters> > 
    {
      typedef tag::matrix  tag;
      typedef typename mtl::dense2D<T, Parameters>::value_type   value_type;
      typedef typename mtl::dense2D<T, Parameters>::size_type    size_type;
    };
    
    template <class T, class Parameters>
    struct category<mtl::compressed2D<T, Parameters> > 
    {
      typedef tag::matrix  tag;
      typedef typename mtl::compressed2D<T, Parameters>::value_type   value_type;
      typedef typename mtl::compressed2D<T, Parameters>::size_type    size_type;
    };
    
    template <class T, class Parameters>
    struct category<mtl::matrix::coordinate2D<T, Parameters> > 
    {
      typedef tag::matrix  tag;
      typedef typename mtl::coordinate2D<T, Parameters>::value_type   value_type;
      typedef typename mtl::coordinate2D<T, Parameters>::size_type    size_type;
    };
    
    template <class T, std::size_t BitMask, class Parameters>
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
  
    template <class Matrix>
    struct category< AMDiS::SolverMatrix<Matrix> > 
    {
      typedef tag::matrix  tag;
    };
    
  
    template <class T>
    struct category< const T > : category< T > {};
    
  
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
    
    // -------------------------------------------------------------------------
    template <class T, class S, class Enable = void>
    struct is_convertible : std::is_convertible<T, S> {};
    
    template <class T, class S>
    struct is_convertible<T, S, typename enable_if< 
				  and_< is_vector<T>, is_vector<S> >
				>::type > 
      : std::is_convertible<typename T::value_type, typename S::value_type> {};
      
    template <class T, class S>
    struct is_convertible<T, S, typename enable_if< 
				  and_< is_matrix<T>, is_matrix<S> >
				>::type > 
      : std::is_convertible<typename T::value_type, typename S::value_type> {};
    
  } // end namespace traits

} // end namespace AMDiS
