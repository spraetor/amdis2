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

#ifndef MTL_PROPERTY_MAP_INCLUDE
#define MTL_PROPERTY_MAP_INCLUDE

#include <boost/numeric/mtl/mtl_fwd.hpp>
#include <boost/numeric/mtl/utility/property_map_impl.hpp>

namespace mtl { namespace traits {    

template <class Matrix> struct row {};
template <class Matrix> struct col {};
template <class Matrix> struct const_value {};
template <class Matrix> struct value {};
template <class Matrix> struct offset {};

// For vectors
template <class Vector> struct index {};

// ===========
// For dense2D
// ===========

template <typename Value, class Parameters>
struct row<mtl::matrix::dense2D<Value, Parameters> >
{
    typedef mtl::detail::indexer_row_ref<mtl::matrix::dense2D<Value, Parameters> > type;
};

template <typename Value, class Parameters>
struct col<mtl::matrix::dense2D<Value, Parameters> >
{
    typedef mtl::detail::indexer_col_ref<mtl::matrix::dense2D<Value, Parameters> > type;
};

template <typename Value, class Parameters>
struct const_value<mtl::matrix::dense2D<Value, Parameters> >
{
    typedef mtl::detail::direct_const_value<mtl::matrix::dense2D<Value, Parameters> > type;
};

template <typename Value, class Parameters>
struct value<mtl::matrix::dense2D<Value, Parameters> >
{
    typedef mtl::detail::direct_value<mtl::matrix::dense2D<Value, Parameters> > type;
};


// ================
// For morton_dense
// ================


template <class Elt, std::size_t BitMask, class Parameters>
struct row<mtl::matrix::morton_dense<Elt, BitMask, Parameters> >
{
    typedef mtl::detail::row_in_key<mtl::matrix::morton_dense<Elt, BitMask, Parameters> > type;
};

template <class Elt, std::size_t BitMask, class Parameters>
struct col<mtl::matrix::morton_dense<Elt, BitMask, Parameters> >
{
    typedef mtl::detail::col_in_key<mtl::matrix::morton_dense<Elt, BitMask, Parameters> > type;
};

template <class Elt, std::size_t BitMask, class Parameters>
struct const_value<mtl::matrix::morton_dense<Elt, BitMask, Parameters> >
{
    typedef mtl::detail::matrix_const_value_ref<mtl::matrix::morton_dense<Elt, BitMask, Parameters> > type;
};

template <class Elt, std::size_t BitMask, class Parameters>
struct value<mtl::matrix::morton_dense<Elt, BitMask, Parameters> >
{
    typedef mtl::detail::matrix_value_ref<mtl::matrix::morton_dense<Elt, BitMask, Parameters> > type;
};


// ================
// For compressed2D
// ================

template <class Elt, class Parameters>
struct row<mtl::matrix::compressed2D<Elt, Parameters> >
{
    typedef typename boost::mpl::if_<
	boost::is_same<typename Parameters::orientation, row_major>
      , mtl::detail::major_in_key<mtl::matrix::compressed2D<Elt, Parameters> >
      , mtl::detail::indexer_minor_ref<mtl::matrix::compressed2D<Elt, Parameters> >
    >::type type;  
};

template <class Elt, class Parameters>
struct col<mtl::matrix::compressed2D<Elt, Parameters> >
{
    typedef typename boost::mpl::if_<
	boost::is_same<typename Parameters::orientation, row_major>
      , mtl::detail::indexer_minor_ref<mtl::matrix::compressed2D<Elt, Parameters> >
      , mtl::detail::major_in_key<mtl::matrix::compressed2D<Elt, Parameters> >
    >::type type;  
};

template <class Elt, class Parameters>
struct const_value<mtl::matrix::compressed2D<Elt, Parameters> >
{
    typedef mtl::detail::matrix_offset_const_value<mtl::matrix::compressed2D<Elt, Parameters> > type;
};

template <class Elt, class Parameters>
struct value<mtl::matrix::compressed2D<Elt, Parameters> >
{
    typedef mtl::detail::matrix_offset_value<mtl::matrix::compressed2D<Elt, Parameters> > type;
};
 
// Offset that corresponds to cursor, e.g. to set values in a matrix with same pattern 
// needed in ILU_0, so far only for compressed2D, could be useful for algos on sparse and dense
template <class Elt, class Parameters>
struct offset<mtl::matrix::compressed2D<Elt, Parameters> >
{
    typedef mtl::detail::offset_from_key<mtl::matrix::compressed2D<Elt, Parameters> > type;
};
  
  
// ================
// For coordinate2D
// ================

template <class Value, class Parameters>
struct row<mtl::matrix::coordinate2D<Value, Parameters> >
{
    typedef mtl::detail::coordinate2D_row<Value, Parameters>   type; 
};

template <class Value, class Parameters>
struct col<mtl::matrix::coordinate2D<Value, Parameters> >
{
    typedef mtl::detail::coordinate2D_col<Value, Parameters>   type; 
};

template <class Value, class Parameters>
struct const_value<mtl::matrix::coordinate2D<Value, Parameters> >
{
    typedef mtl::detail::coordinate2D_const_value<Value, Parameters>   type; 
};

// =================
// For sparse_banded
// =================

template <class Value, class Parameters>
struct row<mtl::matrix::sparse_banded<Value, Parameters> >
{
    typedef mtl::detail::sparse_banded_row<Value, Parameters>   type; 
};

template <class Value, class Parameters>
struct col<mtl::matrix::sparse_banded<Value, Parameters> >
{
    typedef mtl::detail::sparse_banded_col<Value, Parameters>   type; 
};

template <class Value, class Parameters>
struct const_value<mtl::matrix::sparse_banded<Value, Parameters> >
{
    typedef mtl::detail::sparse_banded_const_value<Value, Parameters>   type; 
};

// ==================
// For implicit_dense
// ==================

template <typename Functor>
struct row<mtl::matrix::implicit_dense<Functor> >
{
    typedef mtl::detail::row_in_element_key<mtl::matrix::implicit_dense<Functor> > type;
};

template <typename Functor>
struct col<mtl::matrix::implicit_dense<Functor> >
{
    typedef mtl::detail::col_in_element_key<mtl::matrix::implicit_dense<Functor> > type;
};

template <typename Functor>
struct const_value<mtl::matrix::implicit_dense<Functor> >
{
    typedef mtl::detail::const_value_in_element_key<mtl::matrix::implicit_dense<Functor> > type;
};


// ===============
// For ones_matrix
// ===============

template <typename Value>
struct row<mtl::matrix::ones_matrix<Value> >
  : public row<mtl::matrix::implicit_dense<mtl::matrix::ones_functor<Value> > > 
{};

template <typename Value>
struct col<mtl::matrix::ones_matrix<Value> >
  : public col<mtl::matrix::implicit_dense<mtl::matrix::ones_functor<Value> > > 
{};

template <typename Value>
struct const_value<mtl::matrix::ones_matrix<Value> >
  : public const_value<mtl::matrix::implicit_dense<mtl::matrix::ones_functor<Value> > > 
{};


// ===============
// For hilbert_matrix
// ===============

template <typename Value>
struct row<mtl::matrix::hilbert_matrix<Value> >
  : public row<mtl::matrix::implicit_dense<mtl::matrix::hilbert_functor<Value> > > 
{};

template <typename Value>
struct col<mtl::matrix::hilbert_matrix<Value> >
  : public col<mtl::matrix::implicit_dense<mtl::matrix::hilbert_functor<Value> > > 
{};

template <typename Value>
struct const_value<mtl::matrix::hilbert_matrix<Value> >
  : public const_value<mtl::matrix::implicit_dense<mtl::matrix::hilbert_functor<Value> > > 
{};


// ========================
// For outer_product_matrix
// ========================

template <typename Vector1, typename Vector2>
struct row<mtl::matrix::outer_product_matrix<Vector1, Vector2> >
  : public row<mtl::matrix::implicit_dense<mtl::matrix::outer_product_functor<Vector1, Vector2> > > 
{};

template <typename Vector1, typename Vector2>
struct col<mtl::matrix::outer_product_matrix<Vector1, Vector2> >
  : public col<mtl::matrix::implicit_dense<mtl::matrix::outer_product_functor<Vector1, Vector2> > > 
{};

template <typename Vector1, typename Vector2>
struct const_value<mtl::matrix::outer_product_matrix<Vector1, Vector2> >
  : public const_value<mtl::matrix::implicit_dense<mtl::matrix::outer_product_functor<Vector1, Vector2> > > 
{};


// ====================
// For matrix::indirect
// ====================

template <typename Matrix>
struct row<mtl::matrix::indirect<Matrix> >
{
    typedef mtl::detail::row_in_element_key<mtl::matrix::indirect<Matrix> > type;
};

template <typename Matrix>
struct col<mtl::matrix::indirect<Matrix> >
{
    typedef mtl::detail::col_in_element_key<mtl::matrix::indirect<Matrix> > type;
};

template <typename Matrix>
struct const_value<mtl::matrix::indirect<Matrix> >
{
    typedef mtl::detail::const_value_in_element_key<mtl::matrix::indirect<Matrix> > type;
};


// ================
// For dense_vector
// ================

template <class Elt, class Parameters>
struct index<mtl::vector::dense_vector<Elt, Parameters> >
{
    typedef mtl::detail::index_from_offset< mtl::vector::dense_vector<Elt, Parameters> > type;
};

template <typename Value, class Parameters>
struct const_value<mtl::vector::dense_vector<Value, Parameters> >
{
    typedef mtl::detail::direct_const_value<mtl::vector::dense_vector<Value, Parameters> > type;
};

template <typename Value, class Parameters>
struct value<mtl::vector::dense_vector<Value, Parameters> >
{
    typedef mtl::detail::direct_value<mtl::vector::dense_vector<Value, Parameters> > type;
};
// ================
// For vector::strided_vector_ref
// ================

template <class Elt, class Parameters>
struct index<vector::strided_vector_ref<Elt, Parameters> >
{
    typedef mtl::detail::index_from_offset< vector::strided_vector_ref<Elt, Parameters> > type;
};

template <typename Value, class Parameters>
struct const_value<vector::strided_vector_ref<Value, Parameters> >
{
    typedef mtl::detail::direct_const_value<vector::strided_vector_ref<Value, Parameters> > type;
};

template <typename Value, class Parameters>
struct value<vector::strided_vector_ref<Value, Parameters> >
{
    typedef mtl::detail::direct_value<vector::strided_vector_ref<Value, Parameters> > type;
};




}} // namespace mtl::traits

namespace mtl { namespace matrix {

// Helpers

/// Row map of matrix A
template <typename Matrix>
typename mtl::traits::row<Matrix>::type
inline row_map(const Matrix& A)
{
    return typename mtl::traits::row<Matrix>::type(A);
}

/// Column map of matrix A
template <typename Matrix>
typename mtl::traits::col<Matrix>::type
inline col_map(const Matrix& A)
{
    return typename mtl::traits::col<Matrix>::type(A);
}

/// Constant value map of matrix A
template <typename Matrix>
typename mtl::traits::const_value<Matrix>::type
inline const_value_map(const Matrix& A)
{
    return typename mtl::traits::const_value<Matrix>::type(A);
}

/// Value map of matrix A
template <typename Matrix>
typename mtl::traits::value<Matrix>::type
inline value_map(Matrix& A)
{
    return typename mtl::traits::value<Matrix>::type(A);
}

/// Offset map of matrix A
template <typename Matrix>
typename mtl::traits::offset<Matrix>::type
inline offset_map(const Matrix& A)
{
    return typename mtl::traits::offset<Matrix>::type(A);
}

}} // namespace typename mtl::matrix

namespace mtl { namespace vector {

/// Index map of vector A
template <typename Vector>
typename mtl::traits::index<Vector>::type
inline index_map(const Vector& A)
{
    return typename mtl::traits::index<Vector>::type(A);
}

/// Constant value map of vector A
template <typename Vector>
typename mtl::traits::const_value<Vector>::type
inline const_value_map(const Vector& A)
{
    return typename mtl::traits::const_value<Vector>::type(A);
}

/// Value map of vector A
template <typename Vector>
typename mtl::traits::value<Vector>::type
inline value_map(Vector& A)
{
    return typename mtl::traits::value<Vector>::type(A);
}

}} // namespace typename mtl::vector



#endif // MTL_PROPERTY_MAP_INCLUDE


