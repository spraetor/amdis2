/** \file types.hpp */

#pragma once

// mtl4 headers
#include <boost/numeric/mtl/mtl_fwd.hpp>

// AMDiS headers
#include "AMDiS_fwd.h"
// #include "FixVec.h"
#include "traits/meta_basic.hpp"
#include "traits/scalar_types.hpp"

namespace AMDiS 
{
  namespace traits 
  {
    // dummy type
    class no_valid_type {};
      
    // test for mtl4 types
    // _________________________________________________________________________      
    template<typename T>
    struct is_mtl : false_ {};
    
    template<typename T, typename Param>
    struct is_mtl<mtl::dense_vector<T, Param> > : true_ {};
    
    template<typename T, typename Param>
    struct is_mtl<mtl::dense2D<T, Param> > : true_ {};
    
    template<typename T, typename Param>
    struct is_mtl<mtl::compressed2D<T, Param> > : true_ {};
    
    template<typename T, typename Param>
    struct is_mtl<mtl::matrix::coordinate2D<T, Param> > : true_ {};
    
    template<typename T, std::size_t BitMask, typename Param>
    struct is_mtl<mtl::matrix::morton_dense<T, BitMask, Param> > : true_ {};
    
      
    // test for AMDiS types
    // _________________________________________________________________________     
    template<typename T>
    struct is_amdis : false_ {};
    
    template<typename T>
    struct is_amdis<AMDiS::Vector<T> > : true_ {};
    
    template<typename T>
    struct is_amdis<AMDiS::Matrix<T> > : true_ {};
    
    template<typename T,GeoIndex d>
    struct is_amdis<AMDiS::FixVec<T,d> > : true_ {};
    
    template<typename T>
    struct is_amdis<AMDiS::VectorOfFixVecs<T> > : true_ {};
    
    template<typename T>
    struct is_amdis<AMDiS::MatrixOfFixVecs<T> > : true_ {};
    
    template<typename T>
    struct is_amdis<AMDiS::DimVec<T> > : true_ {};
    
    template<typename T>
    struct is_amdis<AMDiS::DimMat<T> > : true_ {};
    
    template<typename T>
    struct is_amdis<AMDiS::WorldVector<T> > : true_ {};
    
    template<typename T>
    struct is_amdis<AMDiS::WorldMatrix<T> > : true_ {};
    
    template<typename T>
    struct is_amdis<AMDiS::DOFVectorBase<T> > : true_ {};
    
    template<typename T>
    struct is_amdis<AMDiS::DOFVector<T> > : true_ {};
    
    template<>
    struct is_amdis<AMDiS::DOFVectorDOF > : true_ {};
    
    template<>
    struct is_amdis<AMDiS::DOFMatrix > : true_ {};
    
    template<>
    struct is_amdis<AMDiS::DOFContainer > : true_ {};
    
    template<>
    struct is_amdis<AMDiS::SystemVector > : true_ {};
    
  } // end namespace traits
  
} // end namespace AMDiS
