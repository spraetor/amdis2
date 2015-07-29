/** \file MatVecExpr.hpp */

#pragma once

#include <boost/numeric/linear_algebra/identity.hpp>	// mtl::math::zero

#include "traits/basic.hpp"
#include "traits/num_rows.hpp"
#include "traits/num_cols.hpp"
#include "base_expr.hpp"
#include "operations/meta.hpp"

#include "Vector.hpp"

namespace AMDiS {

  template <class MatE>
  struct ColumnIndex
  {
    typedef typename MatE::value_type  value_type;
    typedef typename MatE::size_type    size_type;
    
    ColumnIndex(MatE const& matrix_, size_type row_) 
	: matrix(matrix_), row(row_) { }
  
    inline value_type const& operator()(size_type col) const 
    {
      return matrix(row,col); 
    }
      
  private:
    MatE const& matrix;
    size_type row;
  };
  
  
  // if necessary assign vector_expr to buffer-vector, otherwise store an expr.
  template <class VecE, bool use_buffer>
  struct BufferType {
    typedef typename if_c< use_buffer, 
      typename if_c< (VecE::_SIZE > 0), 
        VectorBase<MemoryBaseStatic<typename VecE::value_type, MAX(0,VecE::_SIZE), 1>, 
				   StaticSizePolicy<MAX(0,VecE::_SIZE)> >,
        VectorBase<MemoryBaseDynamic<typename VecE::value_type, false> > >::type, 
      VecE const& >::type type;
  };

      
  /// Expression with two arguments, that multiplied a matrix_expr with a vector_expr
  template <class MatE, class VecE, bool use_buffer>
  struct MatVecExpr 
      : public VectorExpr< MatVecExpr<MatE, VecE, use_buffer> >
  {
    typedef MatVecExpr                                             self;
    typedef VectorExpr<self>                                  expr_base;
    
    typedef typename VecE::value_type                        value_type;
    typedef typename traits::max_size_type<MatE, VecE>::type  size_type;
	
    typedef MatE                                            matrix_type;
    typedef typename BufferType<VecE, use_buffer>::type     vector_type;
    
    // sizes of the resulting expr.
    static constexpr int _SIZE = MatE::_ROWS;
    static constexpr int _ROWS = MatE::_ROWS;
    static constexpr int _COLS = MAX(VecE::_COLS, 1);
    
  private:
    static constexpr int ARG_COLS = MAX(VecE::_ROWS, MatE::_COLS);
    
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
    
    matrix_type const& get_matrix() const { return matrix; }
    vector_type const& get_vector() const { return vector; }
    
  protected:  
    template <int N>
    value_type reduce(size_type row, int_<N>) const
    {
      using meta::FOR;
      value_type erg = math::zero(value_type());
      FOR<0,N>::inner_product(ColumnIndex<MatE>(matrix,row), vector, erg, 
			      functors::dot_functor<value_type, value_type>());
      return erg;
    }
    
    value_type reduce(size_type r, int_<-1>) const
    {
      value_type erg = math::zero(value_type());	  
      for (size_type c = 0; c < num_cols(matrix); ++c)
	erg += matrix(r,c) * vector(c);
      return erg;
    }
    
  private:
    matrix_type const&  matrix;
    vector_type vector;
  };
  
  
  /// Size of MatVecExpr
  template <class E1, class E2, bool b>
  size_t size(MatVecExpr<E1,E2,b> const& expr)
  {
    return num_rows(expr.get_matrix());
  }
  
  /// number of rows of MatVecExpr
  template <class E1, class E2, bool b>
  size_t num_rows(MatVecExpr<E1,E2,b> const& expr)
  {
    return num_rows(expr.get_matrix());
  }
  
  /// number of columns of MatVecExpr
  template <class E1, class E2, bool b>
  size_t num_cols(MatVecExpr<E1,E2,b> const& expr)
  {
    return 1;
  }
  
  namespace traits {
    
    /// \cond HIDDEN_SYMBOLS
    template <class E1, class E2, bool b>
    struct category<MatVecExpr<E1,E2,b> > 
    {
      typedef typename category<E2>::tag                       tag;
      typedef typename MatVecExpr<E1,E2,b>::value_type  value_type;
      typedef typename MatVecExpr<E1,E2,b>::size_type    size_type;
    };
    /// \endcond
  }
  
} // end namespace AMDiS
