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


#include "DOFIndexed.h"
#include "DOFMatrix.h"

// #include <boost/numeric/mtl/mtl.hpp>
// #include <boost/numeric/mtl/utility/tag.hpp>
// #include <boost/numeric/mtl/utility/category.hpp>
// #include <boost/numeric/linear_algebra/identity.hpp>


// Defining the interface for MTL4
// namespace mtl 
// {
//   // Let MTL4 know that DOFIndexed it is a column vector
//   namespace traits 
//   {
//     template <class T>
//     struct category< AMDiS::DOFIndexed<T> > 
//     {
//       typedef tag::dense_col_vector type;
//     };
//   } // end namespace traits
//   
// 
//   namespace ashape 
//   {
//     template <class T>
//     struct ashape< AMDiS::DOFIndexed<T> > 
//     {
//       typedef cvec<typename ashape<T>::type> type;
//     };
//   } // end namespace ashape
//   
// 
//   // Modelling Collection and MutableCollection
//   template <class T>
//   struct Collection< AMDiS::DOFIndexed<T> > 
//   {
//     typedef typename AMDiS::DOFIndexed<T>::value_type       value_type;
//     typedef typename AMDiS::DOFIndexed<T>::const_reference  const_reference;
//     typedef typename AMDiS::DOFIndexed<T>::size_type        size_type;
//   };
//   
// 
//   template <class T>
//   struct MutableCollection< AMDiS::DOFIndexed<T> > 
//     : public Collection< AMDiS::DOFIndexed<T> > 
//   {
//     typedef typename AMDiS::DOFIndexed<T>::reference       reference;
//   };
// 
// 
// } // namespace mtl
// 
// 
// 
// namespace AMDiS {
// 
//   // Some free functions used in MTL4
// 
//   template <typename T>
//   inline std::size_t size(const AMDiS::DOFIndexed<T>& v)
//   {
//     return v.getSize();
//   }
// 
//   template <typename T>
//   inline std::size_t num_rows(const AMDiS::DOFIndexed<T>& v)
//   {
//     return v.getSize();
//   }
// 
//   template <typename T>
//   inline std::size_t num_cols(const AMDiS::DOFIndexed<T>& v)
//   {
//     return 1;
//   }
//  
// 
//   template <class T>
//   inline void set_to_zero(AMDiS::DOFIndexed<T>& v)
//   {
//     T my_zero; nullify(my_zero);
// 
//     std::fill(v.begin(), v.end(), my_zero);
//   }
// }
