/** \file MatVecExpr.hpp */

#pragma once

#include <boost/numeric/linear_algebra/identity.hpp>	// mtl::math::zero

#include <traits/basic.hpp>
#include <traits/traits_fwd.hpp>

#include "matrix_vector/expr/base_expr.hpp"
#include <operations/meta.hpp>
#include <matrix_vector/Vector.hpp>
#include <Math.hpp>

namespace AMDiS
{
  // helper class for matrix element access
  template <class MatE>
  struct ColumnIndex
  {
    using value_type = Value_t<MatE>;
    using size_type  = Size_t<MatE>;

    ColumnIndex(MatE const& matrix_, size_type row_)
      : matrix(matrix_), row(row_) { }

    value_type const& operator()(size_type col) const
    {
      return matrix(row,col);
    }

  private:
    MatE const& matrix;
    size_type row;
  };


  // if necessary assign vector_expr to buffer-vector, otherwise store an expr.
  template <class VecE, bool use_buffer>
  using BufferType
    = if_then_else< use_buffer,
	if_then_else< (VecE::_SIZE> 0),
	  VectorBase<MemoryBaseStatic<Value_t<VecE>, math::max(0,VecE::_SIZE), 1>,
		     StaticSizePolicy<math::max(0,VecE::_SIZE)> >,
	  VectorBase<MemoryBaseDynamic<Value_t<VecE>, false>> >,
	VecE const& >;


  /// Expression with two arguments, that multiplied a matrix_expr with a vector_expr
  template <class MatE, class VecE, bool use_buffer>
  struct MatVecExpr
    : public VectorExpr<MatVecExpr<MatE, VecE, use_buffer>>
  {
    using Self      = MatVecExpr;
    using expr_base = VectorExpr<Self>;

    using value_type = Value_t<VecE>;
    using size_type  = traits::MaxSizeType<MatE, VecE>;

    using matrix_type = MatE;
    using vector_type = BufferType<VecE, use_buffer>;

    // sizes of the resulting expr.
    static constexpr int _SIZE = MatE::_ROWS;
    static constexpr int _ROWS = MatE::_ROWS;
    static constexpr int _COLS = math::max(VecE::_COLS, 1);

  private:
    static constexpr int ARG_COLS = math::max(VecE::_ROWS, MatE::_COLS);

  public:
    /// constructor takes a matrix expression \p mat and a
    /// vector expression \p vec for the matrix-vector product.
    MatVecExpr(matrix_type const& mat, vector_type const& vec)
      : matrix(mat), vector(vec)
    {
      TEST_EXIT_DBG( num_cols(mat) == num_rows(vec) )("Sizes do not match!\n");
    }

    /// access the elements of an expr.
    value_type operator()(size_type i) const
    {
      return reduce(i, int_<ARG_COLS>());
    }

    matrix_type const& get_matrix() const
    {
      return matrix;
    }
    vector_type const& get_vector() const
    {
      return vector;
    }

  protected:
    // reduce the expression, if length is known at compiletime
    template <int N>
    value_type reduce(size_type row, int_<N>) const
    {
      using meta::FOR;
      using ::math::zero;
      value_type erg = zero(value_type());
      FOR<0,N>::inner_product(ColumnIndex<MatE>(matrix,row), vector, erg,
                              functors::dot_functor<value_type, value_type>());
      return erg;
    }

    // reduce the expression, if length is known only at runtime
    value_type reduce(size_type r, int_<-1>) const
    {
      using ::math::zero;
      value_type erg = zero(value_type());
      for (size_type c = 0; c < num_cols(matrix); ++c)
        erg += matrix(r,c) * vector(c);
      return erg;
    }

  private:
    matrix_type const&  matrix;
    vector_type vector;
  };

  // TODO: generalize MatVec expression to allow a functor passed to the expr.
  //       Example: diagonal(matrix) -> vector_expr


  /// Size of MatVecExpr
  template <class E1, class E2, bool b>
  inline size_t size(MatVecExpr<E1,E2,b> const& expr)
  {
    return num_rows(expr.get_matrix());
  }

  /// number of rows of MatVecExpr
  template <class E1, class E2, bool b>
  inline size_t num_rows(MatVecExpr<E1,E2,b> const& expr)
  {
    return num_rows(expr.get_matrix());
  }

  /// number of columns of MatVecExpr
  template <class E1, class E2, bool b>
  inline size_t num_cols(MatVecExpr<E1,E2,b> const& /*expr*/)
  {
    return 1;
  }

  namespace traits
  {
    /// \cond HIDDEN_SYMBOLS
    template <class E1, class E2, bool b>
    struct category<MatVecExpr<E1,E2,b>>
    {
      using tag = typename category<E2>::tag;
      using value_type = Value_t<MatVecExpr<E1,E2,b>>;
      using size_type  = Size_t<MatVecExpr<E1,E2,b>>;
    };
    /// \endcond

  } // end namespace traits

} // end namespace AMDiS
