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



/** \file types.hpp */

#ifndef AMDIS_TYPE_TRAITS_TYPES_HPP
#define AMDIS_TYPE_TRAITS_TYPES_HPP

#include "AMDiS_fwd.h"
#include "FixVec.h"

// mtl4 types
#include "boost/numeric/mtl/vector/dense_vector.hpp"
#include "boost/numeric/mtl/matrix/dense2D.hpp"
#include "boost/numeric/mtl/matrix/compressed2D.hpp"
#include "boost/numeric/mtl/matrix/coordinate2D.hpp"
#include "boost/numeric/mtl/matrix/morton_dense.hpp"

#include "boost/numeric/ublas/detail/returntype_deduction.hpp"

#include "traits/meta_basic.hpp"
#include "traits/scalar_types.hpp"

namespace AMDiS 
{
  namespace traits 
  {
  
    // dummy type
    typedef boost::numeric::ublas::error_cant_deduce_type no_valid_type;
      
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

#endif // AMDIS_TYPE_TRAITS_TYPES_HPP
