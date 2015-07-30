/** \file num_cols.hpp */

#pragma once

#include "AMDiS_fwd.h"
#include <traits/traits_fwd.hpp>

namespace AMDiS 
{
  namespace traits 
  {
    template <typename T>
    struct num_cols<T, typename enable_if< is_scalar<T> >::type >
    {
	typedef size_t type;
	type operator()(const T& v) const { return 1; }
    };
    
// == vectors ===
    
    template <typename T>
    struct num_cols<T, typename enable_if< and_< is_vector<T>, is_mtl<T> > >::type > 
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
    struct num_cols<T, typename enable_if< and_< is_matrix<T>, is_mtl<T> > >::type > 
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
    
} // end namespace AMDiS
