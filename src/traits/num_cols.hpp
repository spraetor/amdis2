/** \file num_cols.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include "category.hpp"
#include "boost/numeric/mtl/operation/num_cols.hpp"

namespace AMDiS 
{
  namespace traits 
  {
      
    template <typename T, typename Enabled = void>
    struct num_cols 
      : ::mtl::traits::num_cols<T> {}; // import mtl-operation by default
      
    template <typename T>
    struct num_cols<T, typename boost::enable_if< typename has_tag<T, tag::scalar>::type >::type >
    {
	typedef size_t type;
	type operator()(const T& v) const { return 1; }
    };
    
// == vectors ===
    
    template <typename T>
    struct num_cols<T, typename boost::enable_if
    < typename boost::mpl::and_< typename is_vector<T>::type, typename is_mtl<T>::type >::type >::type > 
    {
	typedef typename T::size_type type;
	type operator()(const T& v) const { return ::mtl::vector::num_cols(v); }
    };

    /// size implementation for AMDiS::DOFVector
    template <typename Value>
    struct num_cols< AMDiS::DOFVector<Value> > 
    {
	typedef typename category< AMDiS::DOFVector<Value> >::size_type   type;
	type operator()(const AMDiS::DOFVector<Value>& v) const { return 1; }
    };
  
// === matrices ===
    
    template <typename T>
    struct num_cols<T, typename boost::enable_if
    < typename boost::mpl::and_< typename is_matrix<T>::type, typename is_mtl<T>::type >::type >::type > 
    {
	typedef typename T::size_type type;
	type operator()(const T& v) const { return ::mtl::matrix::num_cols(v); }
    };
      
    /// size implementation for AMDiS::DOFMatrix
    template <>
    struct num_cols< AMDiS::DOFMatrix > 
    {
	typedef category< AMDiS::DOFMatrix >::size_type   type;
	type operator()(const AMDiS::DOFMatrix& v) const { return mtl::matrix::num_cols(v.getBaseMatrix()); }
    };
      
    /// size implementation for AMDiS::BlockMTLMatrix
    template <>
    struct num_cols< AMDiS::BlockMTLMatrix > 
    {
	typedef category< AMDiS::BlockMTLMatrix >::size_type   type;
	type operator()(const AMDiS::BlockMTLMatrix& M) const { return M.m; }
    };
    
// === multi-vector or multi-multi-vector ===
      
    /// size implementation for AMDiS::VectorOfFixVecs
    template <typename FixVecType>
    struct num_cols< AMDiS::VectorOfFixVecs<FixVecType> > 
    {
	typedef typename category< AMDiS::VectorOfFixVecs<FixVecType> >::size_type   type;
	type operator()(const AMDiS::VectorOfFixVecs<FixVecType>& v) const
	{ return v.getSize() == 0 ? 0 : v[0].getSize(); }
    };
    
    /// AMDiS::MatrixOfFixVecs
    template <typename FixVecType>
    struct num_cols< AMDiS::MatrixOfFixVecs<FixVecType> >
    {
	typedef typename category< AMDiS::MatrixOfFixVecs<FixVecType> >::size_type   type;
	type operator()(const AMDiS::MatrixOfFixVecs<FixVecType>& v) const
	{ return v.getNumberOfColumns(); }
    };
    
  } // end namespace traits
  
  template <typename T>
  typename traits::num_cols<T>::type 
  inline num_cols(const T& t)
  {
      return traits::num_cols<T>()(t);
  }
    
} // end namespace AMDiS

#endif // AMDIS_TYPE_TRAITS_NUM_COLS_HPP
