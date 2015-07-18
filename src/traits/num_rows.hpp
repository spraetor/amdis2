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



/** \file num_rows.hpp */

#ifndef AMDIS_TYPE_TRAITS_NUM_ROWS_HPP
#define AMDIS_TYPE_TRAITS_NUM_ROWS_HPP

#include "AMDiS_fwd.h"
#include "category.hpp"
#include "boost/numeric/mtl/operation/num_rows.hpp"

namespace AMDiS 
{
  namespace traits 
  {
  
    template <typename T, typename Enabled = void>
    struct num_rows 
      : ::mtl::traits::num_rows<T> {}; // import mtl-operation by default
      
    template <typename T>
    struct num_rows<T, typename boost::enable_if< typename has_tag<T, tag::scalar>::type >::type >
    {
	typedef size_t type;
	type operator()(const T& v) const { return 1; }
    };
      
// == vectors ===
    
    template <typename T>
    struct num_rows<T, typename boost::enable_if
    < typename boost::mpl::and_< typename is_vector<T>::type, typename is_mtl<T>::type >::type >::type > 
    {
	typedef typename T::size_type type;
	type operator()(const T& v) const { return ::mtl::vector::num_rows(v); }
    };

    /// size implementation for AMDiS::Vector
    template <typename Value>
    struct num_rows< AMDiS::Vector<Value> > 
    {
	typedef typename category< AMDiS::Vector<Value> >::size_type   type;
	type operator()(const AMDiS::Vector<Value>& v) const { return v.getSize(); }
    };
      
    /// size implementation for AMDiS::FixVec
    template <typename Value, AMDiS::GeoIndex d>
    struct num_rows< AMDiS::FixVec<Value, d> >
    {
	typedef typename category< AMDiS::FixVec<Value, d> >::size_type   type;
	type operator()(const AMDiS::FixVec<Value, d>& v) const { return v.getSize(); }
    };
      
    /// size implementation for AMDiS::DimVec
    template <typename Value>
    struct num_rows< AMDiS::DimVec<Value> > 
    {
	typedef typename category< AMDiS::DimVec<Value> >::size_type   type;
	type operator()(const AMDiS::DimVec<Value>& v) const { return v.getSize(); }
    };
      
    /// size implementation for AMDiS::WorldVector
    template <typename Value>
    struct num_rows< AMDiS::WorldVector<Value> > 
    {
	typedef typename category< AMDiS::WorldVector<Value> >::size_type   type;
	type operator()(const AMDiS::WorldVector<Value>& v) const { return v.getSize(); }
    };
      
    /// size implementation for AMDiS::DOFVector
    template <typename Value>
    struct num_rows< AMDiS::DOFVector<Value> > 
    {
	typedef typename category< AMDiS::DOFVector<Value> >::size_type   type;
	type operator()(const AMDiS::DOFVector<Value>& v) const { return v.getUsedSize(); }
    };
  
// === matrices ===
    
    template <typename T>
    struct num_rows<T, typename boost::enable_if
    < typename boost::mpl::and_< typename is_matrix<T>::type, typename is_mtl<T>::type >::type >::type > 
    {
	typedef typename T::size_type type;
	type operator()(const T& v) const { return ::mtl::matrix::num_rows(v); }
    };
      
    /// AMDiS::Matrix
    template <typename Value>
    struct num_rows< AMDiS::Matrix<Value> >
    {
	typedef typename category< AMDiS::Matrix<Value> >::size_type   type;
	type operator()(const AMDiS::Matrix<Value>& v) const 
	{ return v.getNumRows(); }
    };
    
    /// AMDiS::DimMat
    template <typename Value>
    struct num_rows< AMDiS::DimMat<Value> >
    {
	typedef typename category< AMDiS::DimMat<Value> >::size_type   type;
	type operator()(const AMDiS::DimMat<Value>& v) const 
	{ return v.getNumRows(); }
    };
    
    /// AMDiS::WorldMatrix
    template <typename Value>
    struct num_rows< AMDiS::WorldMatrix<Value> >
    {
	typedef typename category< AMDiS::WorldMatrix<Value> >::size_type   type;
	type operator()(const AMDiS::WorldMatrix<Value>& v) const 
	{ return v.getNumRows(); }
    };
      
    /// size implementation for AMDiS::DOFMatrix
    template <>
    struct num_rows< AMDiS::DOFMatrix > 
    {
	typedef category< AMDiS::DOFMatrix >::size_type   type;
	type operator()(const AMDiS::DOFMatrix& v) const { return mtl::matrix::num_rows(v.getBaseMatrix()); }
    };
      
    /// size implementation for AMDiS::BlockMTLMatrix
    template <>
    struct num_rows< AMDiS::BlockMTLMatrix > 
    {
	typedef category< AMDiS::BlockMTLMatrix >::size_type   type;
	type operator()(const AMDiS::BlockMTLMatrix& M) const { return M.n; }
    };
    
// === multi-vector or multi-multi-vector ===
      
    /// size implementation for AMDiS::VectorOfFixVecs
    template <typename FixVecType>
    struct num_rows< AMDiS::VectorOfFixVecs<FixVecType> > 
    {
	typedef typename category< AMDiS::VectorOfFixVecs<FixVecType> >::size_type   type;
	type operator()(const AMDiS::VectorOfFixVecs<FixVecType>& v) const { return v.getSize(); }
    };
    
    /// AMDiS::MatrixOfFixVecs
    template <typename FixVecType>
    struct num_rows< AMDiS::MatrixOfFixVecs<FixVecType> >
    {
	typedef typename category< AMDiS::MatrixOfFixVecs<FixVecType> >::size_type   type;
	type operator()(const AMDiS::MatrixOfFixVecs<FixVecType>& v) const 
	{ return v.getNumberOfRows(); }
    };
    
  } // end namespace traits
  
  template <typename T>
  typename traits::num_rows<T>::type 
  inline num_rows(const T& t)
  {
      return traits::num_rows<T>()(t);
  }
    
} // end namespace AMDiS

#endif // AMDIS_TYPE_TRAITS_NUM_ROWS_HPP
