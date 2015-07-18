// Software License for MTL
// 
// Copyright (c) 2007 The Trustees of Indiana University.
//               2008 Dresden University of Technology and the Trustees of Indiana University.
//               2010 SimuNova UG (haftungsbeschränkt), www.simunova.com.
// All rights reserved.
// Authors: Peter Gottschling and Andrew Lumsdaine
// 
// This file is part of the Matrix Template Library
// 
// See also license.mtl.txt in the distribution.

#ifndef MTL_IS_ROW_MAJOR_INCLUDE
#define MTL_IS_ROW_MAJOR_INCLUDE

#include <boost/mpl/bool.hpp>
#include <boost/mpl/not.hpp>
#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/matrix/parameter.hpp>
#include <boost/numeric/mtl/vector/parameter.hpp>

namespace mtl { namespace traits {

    /// Meta-function whether a tag is row_major or col_major
    /** For convenience, directly applicable to matrix::parameter and vector::parameter. 
	Refined from boost::mpl::true_ or boost::mpl::false_ if defined.
    **/
    template <typename Parameter>
    struct is_row_major {};

    template <>
    struct is_row_major<row_major>
      : boost::mpl::true_ {};

    template <>
    struct is_row_major<col_major>
      : boost::mpl::false_ {};

    template <typename Dimension, bool OnStack, typename SizeType>
    struct is_row_major<vector::parameters<row_major, Dimension, OnStack, SizeType> >
      : boost::mpl::true_ {};

    template <typename T>
    struct is_row_major<const T>
      : is_row_major<T> {};

    template <typename Dimension, bool OnStack, typename SizeType>
    struct is_row_major<vector::parameters<col_major, Dimension, OnStack, SizeType> >
      : boost::mpl::false_ {};

    template <typename Index, typename Dimension, bool OnStack, typename SizeType>
    struct is_row_major<matrix::parameters<row_major, Index, Dimension, OnStack, SizeType> >
      : boost::mpl::true_ {};

    template <typename Index, typename Dimension, bool OnStack, typename SizeType>
    struct is_row_major<matrix::parameters<col_major, Index, Dimension, OnStack, SizeType> >
      : boost::mpl::false_ {};

    template <typename Value, typename Parameters>
    struct is_row_major<vector::dense_vector<Value, Parameters> >
      : is_row_major<Parameters> {};

    template <typename Value, typename Parameters>
    struct is_row_major<vector::strided_vector_ref<Value, Parameters> >
      : is_row_major<Parameters> {};

    template <typename Value, typename Parameters>
    struct is_row_major<mtl::matrix::compressed2D<Value, Parameters> >
      : is_row_major<Parameters> {};

    template <typename Value, typename Parameters>
    struct is_row_major<mtl::matrix::ell_matrix<Value, Parameters> >
      : is_row_major<Parameters> {};

    template <typename Value, typename Parameters>
    struct is_row_major<mtl::matrix::dense2D<Value, Parameters> >
      : is_row_major<Parameters> {};

    template <typename Value, std::size_t Mask, typename Parameters>
    struct is_row_major<mtl::matrix::morton_dense<Value, Mask, Parameters> >
      : is_row_major<Parameters> {};

    template <typename Matrix>
    struct is_row_major<matrix::banded_view<Matrix> >
      : is_row_major<Matrix> {};

    template <typename Matrix> 
    struct is_row_major<matrix::transposed_view<Matrix> >
      : boost::mpl::not_<is_row_major<Matrix> > {};

    template <typename Matrix> 
    struct is_row_major<matrix::hermitian_view<Matrix> >
      : boost::mpl::not_<is_row_major<Matrix> > {};

}} // namespace mtl::traits

#endif // MTL_IS_ROW_MAJOR_INCLUDE
