/** \file size.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "category.hpp"
#include "boost/numeric/mtl/operation/size.hpp"

namespace AMDiS 
{
  namespace traits 
  {
    template <typename T, typename Enabled = void>
    struct size 
      : ::mtl::traits::size<T> {}; // import mtl-operation by default
      
      
    template <typename T>
    struct size<T, typename boost::enable_if< typename has_tag<T, tag::scalar>::type >::type >
    {
	typedef size_t type;
	type operator()(const T& v) const { return 1; }
    };
      
// == vectors ===
    
    template <typename T>
    struct size<T, typename boost::enable_if
    < typename boost::mpl::and_< typename is_vector<T>::type, typename is_mtl<T>::type >::type >::type > 
    {
	typedef typename T::size_type type;
	type operator()(const T& v) const { return ::mtl::vector::size(v); }
    };
      
    /// size implementation for AMDiS::DOFVector
    template <typename Value>
    struct size< AMDiS::DOFVector<Value> > 
    {
	typedef typename category< AMDiS::DOFVector<Value> >::size_type   type;
	type operator()(const AMDiS::DOFVector<Value>& v) const { return v.getUsedSize(); }
    };
  
// === matrices ===
    
    template <typename T>
    struct size<T, typename boost::enable_if
    < typename boost::mpl::and_< typename is_matrix<T>::type, typename is_mtl<T>::type >::type >::type > 
    {
	typedef typename T::size_type type;
	type operator()(const T& v) const { return ::mtl::matrix::size(v); }
    };
      
    /// size implementation for AMDiS::DOFMatrix
    template <>
    struct size< AMDiS::DOFMatrix > 
    {
	typedef category< AMDiS::DOFMatrix >::size_type   type;
	type operator()(const AMDiS::DOFMatrix& v) const { return mtl::matrix::size(v.getBaseMatrix()); }
    };
      
    /// size implementation for AMDiS::BlockMTLMatrix
    template <>
    struct size< AMDiS::BlockMTLMatrix > 
    {
	typedef category< AMDiS::BlockMTLMatrix >::size_type   type;
	type operator()(const AMDiS::BlockMTLMatrix& M) const { return (M.n) * (M.m); }
    };
    
// === multi-vector or multi-matrix ===
      
    /// size implementation for AMDiS::VectorOfFixVecs
    template <typename FixVecType>
    struct size< AMDiS::VectorOfFixVecs<FixVecType> > 
    {
	typedef typename category< AMDiS::VectorOfFixVecs<FixVecType> >::size_type   type;
	type operator()(const AMDiS::VectorOfFixVecs<FixVecType>& v) const 
	{ return v.getSize() * v.getSizeOfFixVec(); }
    };
    
    /// AMDiS::MatrixOfFixVecs
    template <typename FixVecType>
    struct size< AMDiS::MatrixOfFixVecs<FixVecType> >
    {
	typedef typename category< AMDiS::MatrixOfFixVecs<FixVecType> >::size_type   type;
	type operator()(const AMDiS::MatrixOfFixVecs<FixVecType>& v) const 
	{ 
	  return v.getNumberOfRows() == 0 ? 0 : v.getNumberOfRows() * v.getNumberOfColumns() * v[0].getSizeOfFixVec(); 
	}
    };
    
  } // end namespace traits
  
  template <typename T>
  typename traits::size<T>::type 
  inline size(const T& t)
  {
      return traits::size<T>()(t);
  }
    
} // end namespace AMDiS
