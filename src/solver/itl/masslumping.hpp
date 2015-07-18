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

#ifndef ITL_PC_MASSLUMPING_INCLUDE
#define ITL_PC_MASSLUMPING_INCLUDE

#include <boost/numeric/linear_algebra/inverse.hpp>

#include <boost/numeric/mtl/vector/dense_vector.hpp>
#include <boost/numeric/mtl/utility/exception.hpp>
#include <boost/numeric/mtl/utility/tag.hpp>
#include <boost/numeric/mtl/detail/base_cursor.hpp>
#include <boost/numeric/mtl/concept/collection.hpp>
#include <boost/numeric/mtl/operation/resource.hpp>
#include <boost/numeric/mtl/interface/vpt.hpp>
#include <boost/numeric/itl/pc/solver.hpp>

namespace itl { namespace pc {

/// Diagonal Preconditioner
template <typename Matrix, typename Value= typename mtl::Collection<Matrix>::value_type>
class masslumping
{
  public:
    typedef Value                                         value_type;
    typedef typename mtl::Collection<Matrix>::size_type   size_type;
    typedef masslumping                                   self;

    /// Constructor takes matrix reference
    explicit masslumping(const Matrix& A) : inv_diag(num_rows(A))
    {
// 	mtl::vampir_trace<5050> tracer;
	MTL_THROW_IF(num_rows(A) != num_cols(A), mtl::matrix_not_square());
	using namespace mtl; using namespace mtl::tag; 
	using mtl::traits::range_generator; 
	using math::reciprocal; 
	
	typedef typename range_generator<tag::row, Matrix>::type c_type;
	typedef typename range_generator<tag::nz, c_type>::type  ic_type;

	typename mtl::traits::row<Matrix>::type row(A);
	typename mtl::traits::const_value<Matrix>::type value(A);
	
	value_type row_sum;
	c_type cursor(begin<tag::row>(A));
	for (c_type cend(end<tag::row>(A)); cursor != cend; ++cursor) {
	  row_sum = math::zero(value_type());
	  ic_type icursor(begin<tag::nz>(cursor));
	  size_type r = row(icursor);
	
	  for (ic_type icend(end<tag::nz>(cursor)); icursor != icend; ++icursor)
	    row_sum += value(*icursor);
	  
	  inv_diag[r]= reciprocal(row_sum);
	}
    }

    /// Member function solve, better use free function solve
    template <typename Vector>
    Vector solve(const Vector& x) const
    {
	Vector y(resource(x));
	solve(x, y);
	return y;
    }

    template <typename VectorIn, typename VectorOut>
    void solve(const VectorIn& x, VectorOut& y) const
    {
	mtl::vampir_trace<5051> tracer;
	y.checked_change_resource(x);
	MTL_THROW_IF(size(x) != size(inv_diag), mtl::incompatible_size());
	for (size_type i= 0; i < size(inv_diag); ++i)
	    y[i]= inv_diag[i] * x[i];
    }

    /// Member function for solving adjoint problem, better use free function adjoint_solve
    template <typename Vector>
    Vector adjoint_solve(const Vector& x) const
    {
	Vector y(resource(x));
	adjoint_solve(x, y);
	return y;
    }

    template <typename VectorIn, typename VectorOut>
    void adjoint_solve(const VectorIn& x, VectorOut& y) const
    {
	using mtl::conj;
	y.checked_change_resource(x);
	MTL_THROW_IF(size(x) != size(inv_diag), mtl::incompatible_size());
	for (size_type i= 0; i < size(inv_diag); ++i)
	    y[i]= conj(inv_diag[i]) * x[i];
    }

 protected:
    mtl::vector::dense_vector<value_type>    inv_diag;
}; 

/// Solve approximately a sparse system in terms of inverse lumped diagonal
template <typename Matrix, typename Vector>
solver<masslumping<Matrix>, Vector, false>
inline solve(const masslumping<Matrix>& P, const Vector& x)
{
    return solver<masslumping<Matrix>, Vector, false>(P, x);
}

/// Solve approximately the adjoint of a sparse system in terms of inverse lumped diagonal
template <typename Matrix, typename Vector>
solver<masslumping<Matrix>, Vector, true>
inline adjoint_solve(const masslumping<Matrix>& P, const Vector& x)
{
    return solver<masslumping<Matrix>, Vector, true>(P, x);
}


}} // namespace itl::pc

#endif // ITL_PC_MASSLUMPING_INCLUDE
