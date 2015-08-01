/** \file num_rows.hpp */

#pragma once

#include <AMDiS_fwd.h>
#include <traits/traits_fwd.hpp>

namespace AMDiS 
{
  namespace traits 
  {
    template <class T>
    struct num_rows<T, typename enable_if< is_scalar<T> >::type >
    {
      typedef size_t type;
      type operator()(const T& v) const { return 1; }
    };
      
// == vectors ===
    
    template <class T>
    struct num_rows<T, typename enable_if< and_< is_vector<T>, is_mtl<T> > >::type > 
    {
      typedef typename T::size_type type;
      type operator()(const T& v) const { return ::mtl::vector::num_rows(v); }
    };
      
    /// size implementation for AMDiS::DOFVector
    template <class Value>
    struct num_rows< AMDiS::DOFVector<Value> > 
    {
      typedef typename category< AMDiS::DOFVector<Value> >::size_type   type;
      type operator()(const AMDiS::DOFVector<Value>& v) const { return v.getUsedSize(); }
    };
  
// === matrices ===
    
    template <class T>
    struct num_rows<T, typename enable_if< and_< is_matrix<T>, is_mtl<T> > >::type > 
    {
      typedef typename T::size_type type;
      type operator()(const T& v) const { return ::mtl::matrix::num_rows(v); }
    };
    
// === multi-vector or multi-multi-vector ===
      
    /// size implementation for AMDiS::VectorOfFixVecs
    template <class FixVecType>
    struct num_rows< AMDiS::VectorOfFixVecs<FixVecType> > 
    {
      typedef typename category< AMDiS::VectorOfFixVecs<FixVecType> >::size_type   type;
      type operator()(const AMDiS::VectorOfFixVecs<FixVecType>& v) const { return v.getSize(); }
    };
    
    /// AMDiS::MatrixOfFixVecs
    template <class FixVecType>
    struct num_rows< AMDiS::MatrixOfFixVecs<FixVecType> >
    {
      typedef typename category< AMDiS::MatrixOfFixVecs<FixVecType> >::size_type   type;
      type operator()(const AMDiS::MatrixOfFixVecs<FixVecType>& v) const 
      { return v.getNumberOfRows(); }
    };
    
  } // end namespace traits
    
} // end namespace AMDiS
