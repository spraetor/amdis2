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

#ifndef MTL_TRAITS_TRANSPOSED_MATRIX_TYPE_INCLUDE
#define MTL_TRAITS_TRANSPOSED_MATRIX_TYPE_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/transposed_orientation.hpp>

namespace mtl { namespace traits {

template <class T> struct transposed_matrix_parameter {};

template <typename O, typename I, typename D, bool S, typename ST>
struct transposed_matrix_parameter<matrix::parameters<O, I, D, S, ST> >
{
    typedef matrix::parameters<typename transposed_orientation<O>::type, I, D, S, ST>  type;
};

template <class T> struct transposed_matrix_type {};

template <typename Value, typename Parameters>
struct transposed_matrix_type<matrix::dense2D<Value, Parameters> >
{
    typedef matrix::dense2D<Value, typename transposed_matrix_parameter<Parameters>::type> type;
};

template <typename Value, typename Parameters>
struct transposed_matrix_type<matrix::compressed2D<Value, Parameters> >
{
    typedef matrix::compressed2D<Value, typename transposed_matrix_parameter<Parameters>::type> type;
};

template <typename Value, std::size_t Mask, typename Parameters>
struct transposed_matrix_type<matrix::morton_dense<Value, Mask, Parameters> >
{
    typedef matrix::morton_dense<Value, Mask, typename transposed_matrix_parameter<Parameters>::type> type;
};




template <class T> struct transposed_sparse_matrix_type {};

template <typename Value, typename Parameters>
struct transposed_sparse_matrix_type<matrix::compressed2D<Value, Parameters> >
{
    typedef matrix::compressed2D<Value, typename transposed_matrix_parameter<Parameters>::type> type;
};

template <typename Matrix>
struct transposed_sparse_matrix_type<matrix::banded_view<Matrix> >
{
    typedef typename transposed_sparse_matrix_type<Matrix>::type type;
};


template <typename Value, typename Parameters>
struct transposed_sparse_matrix_type<matrix::transposed_view<matrix::compressed2D<Value, Parameters> > >
{
    typedef matrix::compressed2D<Value, Parameters> type;
};

template <typename Value, typename Parameters>
struct transposed_sparse_matrix_type<matrix::transposed_view<const matrix::compressed2D<Value, Parameters> > >
{
    typedef matrix::compressed2D<Value, Parameters> type;
};

}} // namespace mtl::traits

#endif // MTL_TRAITS_TRANSPOSED_MATRIX_TYPE_INCLUDE
